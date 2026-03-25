"""Ollama LLM provider using the native ollama SDK."""

from __future__ import annotations

import logging
from collections.abc import Iterator

from epic.models import Message
from epic.providers.base import LLMProvider

logger = logging.getLogger(__name__)

# Singleton clients keyed by host URL — connection reuse across sessions
_clients: dict[str, object] = {}


def _get_client(host: str):
    """Get or create a shared Ollama client for the given host."""
    if host not in _clients:
        from ollama import Client as OllamaClient

        _clients[host] = OllamaClient(host=host)
        logger.info("Created Ollama client for %s", host)
    return _clients[host]


class OllamaProvider(LLMProvider):
    """Ollama provider using the official ollama Python SDK."""

    def __init__(
        self,
        host: str = "http://localhost:11434",
        model: str = "qwen3.5:4b",
        keep_alive: str = "10m",
    ) -> None:
        self._client = _get_client(host)
        self._model = model
        self._host = host
        self._keep_alive = keep_alive

    @property
    def name(self) -> str:
        return "ollama"

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
        sdk_messages: list[dict] = [{"role": "system", "content": system_prompt}]
        for m in messages:
            sdk_messages.append({"role": m.role.value, "content": m.content})

        logger.debug("Ollama request: model=%s, messages=%d", self._model, len(sdk_messages))

        stream = self._client.chat(
            model=self._model,
            messages=sdk_messages,
            stream=True,
            keep_alive=self._keep_alive,
            options={
                "temperature": temperature,
                "num_predict": max_tokens,
            },
        )
        for chunk in stream:
            content = chunk.message.content
            if content:
                yield content

    def list_models(self) -> list[str]:
        try:
            response = self._client.list()
            return [m.model for m in response.models]
        except Exception as e:
            logger.warning("Failed to list Ollama models: %s", e)
            return [self._model]
