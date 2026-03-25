"""Anthropic LLM provider using the native anthropic SDK."""

from __future__ import annotations

import logging
from collections.abc import Iterator

from epic.models import Message, Role
from epic.providers.base import LLMProvider

logger = logging.getLogger(__name__)


class AnthropicProvider(LLMProvider):
    """Anthropic provider using the official anthropic Python SDK."""

    KNOWN_MODELS = [
        "claude-sonnet-4-20250514",
        "claude-opus-4-20250514",
        "claude-haiku-4-20250414",
    ]

    def __init__(self, api_key: str, model: str = "claude-sonnet-4-20250514") -> None:
        import anthropic

        self._client = anthropic.Anthropic(api_key=api_key)
        self._model = model

    @property
    def name(self) -> str:
        return "anthropic"

    @property
    def model_name(self) -> str:
        return self._model

    def stream_chat(
        self,
        messages: list[Message],
        system_prompt: str,
        temperature: float = 0.2,
        max_tokens: int = 4096,
    ) -> Iterator[str]:
        # Anthropic: system prompt is a top-level parameter, not in messages.
        # Filter out any system messages from the list.
        sdk_messages = [
            {"role": m.role.value, "content": m.content}
            for m in messages
            if m.role != Role.SYSTEM
        ]

        logger.debug("Anthropic request: model=%s, messages=%d", self._model, len(sdk_messages))

        with self._client.messages.stream(
            model=self._model,
            system=system_prompt,
            messages=sdk_messages,
            temperature=temperature,
            max_tokens=max_tokens,
        ) as stream:
            for text in stream.text_stream:
                yield text

    def list_models(self) -> list[str]:
        return list(self.KNOWN_MODELS)
