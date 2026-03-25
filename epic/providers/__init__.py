"""LLM provider registry with lazy imports."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from epic.config import Settings
    from epic.providers.base import LLMProvider


def get_provider(
    name: str,
    settings: "Settings",
    model: str | None = None,
) -> "LLMProvider":
    """Factory function. Lazily imports provider SDKs so only the needed one must be installed."""
    if name == "openai":
        from epic.providers.openai_provider import OpenAIProvider

        api_key = settings.openai_api_key
        if not api_key:
            raise ValueError("OpenAI API key not configured. Set OPENAI_API_KEY in .env")
        return OpenAIProvider(api_key=api_key, model=model or "gpt-4o")

    elif name == "anthropic":
        from epic.providers.anthropic_provider import AnthropicProvider

        api_key = settings.anthropic_api_key
        if not api_key:
            raise ValueError("Anthropic API key not configured. Set ANTHROPIC_API_KEY in .env")
        return AnthropicProvider(api_key=api_key, model=model or "claude-sonnet-4-20250514")

    elif name == "ollama":
        from epic.providers.ollama_provider import OllamaProvider

        return OllamaProvider(
            host=settings.ollama_host,
            model=model or settings.ollama_default_model,
            keep_alive=settings.ollama_keep_alive,
        )

    else:
        raise ValueError(f"Unknown provider: {name!r}. Available: openai, anthropic, ollama")
