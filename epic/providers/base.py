"""Abstract base class for LLM providers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator

from epic.models import Message


class LLMProvider(ABC):
    """Provider-neutral base class. No SDK concepts leak into this interface."""

    @abstractmethod
    def stream_chat(
        self,
        messages: list[Message],
        system_prompt: str,
        temperature: float = 0.2,
        max_tokens: int = 4096,
    ) -> Iterator[str]:
        """Yield text chunks from the model.

        Args:
            messages: Conversation history (user/assistant only).
            system_prompt: System instructions, handled per-provider.
            temperature: Sampling temperature.
            max_tokens: Maximum tokens to generate.

        Yields:
            Text chunks as they arrive from the model.
        """
        ...

    @abstractmethod
    def list_models(self) -> list[str]:
        """Return available model identifiers for this provider."""
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Provider identifier string (e.g. 'openai', 'anthropic', 'ollama')."""
        ...

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Currently configured model identifier."""
        ...
