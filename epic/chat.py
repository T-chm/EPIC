"""ChatEngine: central orchestrator used by both the API and CLI."""

from __future__ import annotations

import asyncio
import logging
import re
from collections.abc import AsyncIterator, Iterator
from concurrent.futures import ThreadPoolExecutor

from epic.models import CodeResult, Message, Role, StreamEvent
from epic.prompts import SYSTEM_PROMPT
from epic.providers.base import LLMProvider
from epic.tools.interpreter import CodeInterpreter

logger = logging.getLogger(__name__)

# Regex to extract ```python ... ``` code blocks
_CODE_BLOCK_RE = re.compile(r"```python\n(.*?)\n```", re.DOTALL)

_INTERPRETATION_INSTRUCTION = (
    "Above is the actual code execution output. Give a brief interpretation (3-5 sentences max):\n"
    "- State whether the design passes or fails each criterion\n"
    "- Flag specific values that are outside range\n"
    "- If there are issues, suggest ONE concrete next step\n"
    "Do NOT repeat the data in tables. Do NOT generate new code. Be concise."
)

# Sentinel to signal end of stream
_DONE = object()

# Shared thread pool for blocking provider calls
_executor = ThreadPoolExecutor(max_workers=8)


class ChatEngine:
    """Per-session orchestrator that streams LLM responses and executes code.

    Used identically by the FastAPI WebSocket handler and the CLI.
    """

    def __init__(
        self,
        provider: LLMProvider,
        system_prompt: str = SYSTEM_PROMPT,
        auto_execute: bool = True,
        temperature: float = 0.2,
        max_tokens: int = 4096,
        compression_threshold: int = 20,
        compression_keep_recent: int = 4,
    ) -> None:
        self._provider = provider
        self._system_prompt = system_prompt
        self._auto_execute = auto_execute
        self._temperature = temperature
        self._max_tokens = max_tokens
        self._compression_threshold = compression_threshold
        self._compression_keep_recent = compression_keep_recent
        self._history: list[Message] = []
        self._interpreter = CodeInterpreter()
        self._setup_interpreter()

    def _setup_interpreter(self) -> None:
        """Pre-load bioinformatics libraries into the interpreter namespace."""
        try:
            import primer3
            import Bio
            from Bio import SeqUtils
            from Bio.Seq import Seq
            from Bio.SeqUtils import MeltingTemp as mt
            from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_GC, Tm_Wallace
            from Bio.SeqUtils import gc_fraction

            from epic.tools.mafft import run_mafft

            self._interpreter.inject_globals(
                {
                    "Bio": Bio,
                    "Seq": Seq,
                    "mt": mt,
                    "MeltingTemp": mt,
                    "Tm_NN": Tm_NN,
                    "Tm_GC": Tm_GC,
                    "Tm_Wallace": Tm_Wallace,
                    "gc_fraction": gc_fraction,
                    "SeqUtils": SeqUtils,
                    "primer3": primer3,
                    "design_primers": primer3.design_primers,
                    "run_mafft": run_mafft,
                }
            )
        except ImportError as e:
            logger.warning("Could not preload some bioinformatics libraries: %s", e)

    @property
    def provider(self) -> LLMProvider:
        return self._provider

    @property
    def history(self) -> list[Message]:
        return list(self._history)

    def switch_provider(self, provider: LLMProvider) -> None:
        """Switch LLM provider mid-session, preserving conversation history."""
        logger.info(
            "Switching provider from %s/%s to %s/%s",
            self._provider.name,
            self._provider.model_name,
            provider.name,
            provider.model_name,
        )
        self._provider = provider

    def clear_history(self) -> None:
        """Clear conversation history."""
        self._history.clear()

    def set_temperature(self, temperature: float) -> None:
        self._temperature = temperature

    def set_auto_execute(self, auto_execute: bool) -> None:
        self._auto_execute = auto_execute

    def set_max_tokens(self, max_tokens: int) -> None:
        self._max_tokens = max_tokens

    def _maybe_compress_history(self) -> None:
        """Compress old conversation history into a summary when threshold is exceeded.

        Keeps the most recent messages intact and summarizes older ones into a
        single system-role message, preventing context window overflow.
        """
        if len(self._history) < self._compression_threshold:
            return

        keep = self._compression_keep_recent
        to_compress = self._history[:-keep]
        to_keep = self._history[-keep:]

        # Build a text representation of old messages for summarization
        lines = []
        for m in to_compress:
            lines.append(f"{m.role.value}: {m.content[:500]}")
        old_text = "\n".join(lines)

        # Use a non-streaming LLM call to generate the summary
        summary_prompt = (
            "Summarize the key information from this conversation history in 3-5 sentences. "
            "Focus on: what was designed, what parameters were used, what results were obtained.\n\n"
            f"{old_text}"
        )
        summary_messages = [Message(role=Role.USER, content=summary_prompt)]

        try:
            summary_chunks = list(self._provider.stream_chat(
                messages=summary_messages,
                system_prompt="You are a concise summarizer. Output only the summary.",
                temperature=0.1,
                max_tokens=500,
            ))
            summary = "".join(summary_chunks)
        except Exception as e:
            logger.warning("History compression failed: %s", e)
            return

        # Replace history with summary + recent messages
        self._history = [
            Message(role=Role.ASSISTANT, content=f"[Previous conversation summary]\n{summary}"),
            *to_keep,
        ]
        logger.info(
            "Compressed %d messages into summary, keeping %d recent",
            len(to_compress), keep,
        )

    def execute_code(self, code_str: str) -> CodeResult:
        """Execute Python code directly and return the result."""
        stdout, stderr, success = self._interpreter.execute(code_str)
        return CodeResult(code=code_str, stdout=stdout, stderr=stderr, success=success)

    def _stream_llm(self, messages: list[Message]) -> Iterator[tuple[str, str]]:
        """Stream one LLM turn. Yields (event_type, data) tuples.

        Detects <think>...</think> blocks (from qwen3.5 and similar models)
        and yields them as thinking_start/thinking_token/thinking_end events,
        keeping them out of the main response text.
        """
        full_response = ""
        in_thinking = False
        buffer = ""

        try:
            for chunk in self._provider.stream_chat(
                messages=messages,
                system_prompt=self._system_prompt,
                temperature=self._temperature,
                max_tokens=self._max_tokens,
            ):
                buffer += chunk

                while buffer:
                    if not in_thinking:
                        think_idx = buffer.find("<think>")
                        if think_idx == -1:
                            # No think tag — emit all as tokens
                            full_response += buffer
                            yield ("token", buffer)
                            buffer = ""
                        else:
                            # Emit text before <think>
                            if think_idx > 0:
                                pre = buffer[:think_idx]
                                full_response += pre
                                yield ("token", pre)
                            yield ("thinking_start", "")
                            in_thinking = True
                            buffer = buffer[think_idx + len("<think>"):]
                    else:
                        end_idx = buffer.find("</think>")
                        if end_idx == -1:
                            # Still in thinking — emit what we have
                            if buffer:
                                yield ("thinking_token", buffer)
                            buffer = ""
                        else:
                            # Emit thinking content before </think>
                            if end_idx > 0:
                                yield ("thinking_token", buffer[:end_idx])
                            yield ("thinking_end", "")
                            in_thinking = False
                            buffer = buffer[end_idx + len("</think>"):]

        except Exception as e:
            logger.error("Provider error: %s", e)
            yield ("error", str(e))
            return

        # Flush any remaining buffer
        if buffer:
            if in_thinking:
                yield ("thinking_token", buffer)
                yield ("thinking_end", "")
            else:
                full_response += buffer
                yield ("token", buffer)

        yield ("_full_response", full_response)

    def stream_response(self, user_input: str) -> Iterator[StreamEvent]:
        """Synchronous streaming generator. Used by CLI directly.

        Two-turn architecture:
        Turn 1: LLM generates plan + code.
        Code is extracted and executed.
        Turn 2: Code output is fed back to the LLM for interpretation.

        Yields StreamEvent objects:
        - token: text chunks from the LLM
        - code_detected: when a Python code block is found
        - code_executing: before running code
        - code_result: after running code
        - done: final event
        - error: if something goes wrong
        """
        # Compress old history if needed before adding new message
        self._maybe_compress_history()

        self._history.append(Message(role=Role.USER, content=user_input))

        # ── Direct dispatch: bypass LLM for obvious QDB/PCR requests ──
        from epic.tools.direct_dispatch import try_direct_dispatch

        direct_code = try_direct_dispatch(user_input)
        if direct_code and self._auto_execute:
            logger.info("Direct dispatch: bypassing LLM for QDB design")
            preamble = "Executing QDB design directly...\n"
            yield StreamEvent(event_type="token", data=preamble)
            full_response = preamble
            self._history.append(Message(role=Role.ASSISTANT, content=full_response))

            # Execute the generated code
            code_outputs: list[str] = []
            yield StreamEvent(event_type="code_detected", data=direct_code)
            yield StreamEvent(event_type="code_executing", data=direct_code)
            result = self.execute_code(direct_code)
            output = result.stdout if result.success else f"Error: {result.stderr}"
            yield StreamEvent(event_type="code_result", data=output)
            code_outputs.append(output)

            # Skip to Turn 2 for interpretation
            results_text = f"Code execution output:\n```\n{output}\n```"
            followup = f"{results_text}\n\n{_INTERPRETATION_INSTRUCTION}"
            self._history.append(Message(role=Role.USER, content=followup))

            interpretation = ""
            for event_type, data in self._stream_llm(self._history):
                if event_type in ("token", "thinking_start", "thinking_token", "thinking_end"):
                    yield StreamEvent(event_type=event_type, data=data)
                elif event_type == "error":
                    yield StreamEvent(event_type="error", data=data)
                    return
                elif event_type == "_full_response":
                    interpretation = data

            self._history.append(Message(role=Role.ASSISTANT, content=interpretation))
            yield StreamEvent(event_type="done", data=full_response)
            return

        # ── Turn 1: LLM generates plan + code ──
        full_response = ""
        for event_type, data in self._stream_llm(self._history):
            if event_type in ("token", "thinking_start", "thinking_token", "thinking_end"):
                yield StreamEvent(event_type=event_type, data=data)
            elif event_type == "error":
                yield StreamEvent(event_type="error", data=data)
                return
            elif event_type == "_full_response":
                full_response = data

        self._history.append(Message(role=Role.ASSISTANT, content=full_response))

        # ── Execute code blocks ──
        code_blocks = _CODE_BLOCK_RE.findall(full_response)
        code_outputs: list[str] = []

        for code_str in code_blocks:
            code_str = code_str.strip()
            if not code_str:
                continue
            yield StreamEvent(event_type="code_detected", data=code_str)
            if self._auto_execute:
                yield StreamEvent(event_type="code_executing", data=code_str)
                result = self.execute_code(code_str)
                output = result.stdout if result.success else f"Error: {result.stderr}"
                yield StreamEvent(event_type="code_result", data=output)
                code_outputs.append(output)

        # ── Turn 2: Feed code output back for interpretation ──
        if code_outputs and self._auto_execute:
            results_text = "\n\n".join(
                f"Code execution output:\n```\n{out}\n```" for out in code_outputs
            )
            followup = f"{results_text}\n\n{_INTERPRETATION_INSTRUCTION}"
            self._history.append(Message(role=Role.USER, content=followup))

            interpretation = ""
            for event_type, data in self._stream_llm(self._history):
                if event_type in ("token", "thinking_start", "thinking_token", "thinking_end"):
                    yield StreamEvent(event_type=event_type, data=data)
                elif event_type == "error":
                    yield StreamEvent(event_type="error", data=data)
                    return
                elif event_type == "_full_response":
                    interpretation = data

            self._history.append(Message(role=Role.ASSISTANT, content=interpretation))

        yield StreamEvent(event_type="done", data=full_response)

    async def astream_response(self, user_input: str) -> AsyncIterator[StreamEvent]:
        """Async wrapper for FastAPI. Runs the sync generator in a background thread
        and forwards events via an asyncio.Queue.

        Token events are batched: the producer pushes raw events, and the consumer
        coalesces consecutive token events into a single yield every ~50ms. This
        dramatically reduces WebSocket messages and React re-renders.
        """
        queue: asyncio.Queue[StreamEvent | object] = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def _produce() -> None:
            """Run in thread pool — consumes sync generator, pushes to queue."""
            try:
                for event in self.stream_response(user_input):
                    loop.call_soon_threadsafe(queue.put_nowait, event)
            except Exception as e:
                loop.call_soon_threadsafe(
                    queue.put_nowait,
                    StreamEvent(event_type="error", data=str(e)),
                )
            finally:
                loop.call_soon_threadsafe(queue.put_nowait, _DONE)

        # Start producer in thread pool
        fut = loop.run_in_executor(_executor, _produce)

        # Consume from queue, batching consecutive token events
        BATCH_INTERVAL = 0.05  # 50ms batching window
        token_buffer = ""

        try:
            while True:
                # Wait for at least one item
                try:
                    item = await asyncio.wait_for(queue.get(), timeout=BATCH_INTERVAL)
                except asyncio.TimeoutError:
                    # Flush any buffered tokens on timeout
                    if token_buffer:
                        yield StreamEvent(event_type="token", data=token_buffer)
                        token_buffer = ""
                    continue

                if item is _DONE:
                    # Flush remaining buffer before finishing
                    if token_buffer:
                        yield StreamEvent(event_type="token", data=token_buffer)
                        token_buffer = ""
                    break

                event: StreamEvent = item  # type: ignore[assignment]

                if event.event_type == "token":
                    token_buffer += event.data
                    # Drain any additional token events already in the queue
                    while not queue.empty():
                        try:
                            peek = queue.get_nowait()
                        except asyncio.QueueEmpty:
                            break
                        if peek is _DONE:
                            if token_buffer:
                                yield StreamEvent(event_type="token", data=token_buffer)
                                token_buffer = ""
                            # Re-signal done for the outer loop
                            await queue.put(_DONE)
                            break
                        if peek.event_type == "token":  # type: ignore[union-attr]
                            token_buffer += peek.data  # type: ignore[union-attr]
                        else:
                            # Non-token event: flush buffer first, then put this back
                            if token_buffer:
                                yield StreamEvent(event_type="token", data=token_buffer)
                                token_buffer = ""
                            yield peek  # type: ignore[misc]
                            break
                else:
                    # Non-token event: flush any buffered tokens first
                    if token_buffer:
                        yield StreamEvent(event_type="token", data=token_buffer)
                        token_buffer = ""
                    yield event
        finally:
            await fut
