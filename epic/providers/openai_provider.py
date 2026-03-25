"""OpenAI LLM provider using the native openai SDK."""

from __future__ import annotations

import logging
from collections.abc import Iterator

from epic.models import Message
from epic.providers.base import LLMProvider

logger = logging.getLogger(__name__)


class OpenAIProvider(LLMProvider):
    """OpenAI provider using the official openai Python SDK."""

    def __init__(self, api_key: str, model: str = "gpt-4o") -> None:
        import openai

        self._client = openai.OpenAI(api_key=api_key)
        self._model = model

    @property
    def name(self) -> str:
        return "openai"

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

        logger.debug("OpenAI request: model=%s, messages=%d", self._model, len(sdk_messages))

        stream = self._client.chat.completions.create(
            model=self._model,
            messages=sdk_messages,
            temperature=temperature,
            max_tokens=max_tokens,
            stream=True,
        )
        for chunk in stream:
            delta = chunk.choices[0].delta
            if delta and delta.content:
                yield delta.content

    def list_models(self) -> list[str]:
        try:
            models = self._client.models.list()
            return sorted(m.id for m in models.data if "gpt" in m.id)
        except Exception as e:
            logger.warning("Failed to list OpenAI models: %s", e)
            return ["gpt-4o", "gpt-4o-mini", "gpt-4-turbo"]
