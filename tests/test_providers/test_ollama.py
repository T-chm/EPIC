"""Tests for Ollama provider — verify it uses native ollama SDK, NOT openai."""

from __future__ import annotations

import importlib
import inspect
from unittest.mock import MagicMock, patch

from epic.models import Message, Role


def test_ollama_does_not_import_openai():
    """The ollama provider module must NOT import from openai."""
    import epic.providers.ollama_provider as mod

    source = inspect.getsource(mod)
    assert "import openai" not in source
    assert "from openai" not in source


@patch("epic.providers.ollama_provider.OllamaClient", create=True)
def test_ollama_provider_uses_native_sdk(MockOllamaClient):
    """Ollama provider should call ollama.Client.chat(), not openai."""
    # Need to patch at import time
    with patch.dict("sys.modules", {}):
        mock_client = MagicMock()
        MockOllamaClient.return_value = mock_client

        # Ollama SDK returns ChatResponse objects with attribute access
        chunk1 = MagicMock()
        chunk1.message.content = "Hello"
        chunk2 = MagicMock()
        chunk2.message.content = " world"
        mock_client.chat.return_value = [chunk1, chunk2]

        from epic.providers.ollama_provider import OllamaProvider

        with patch.object(OllamaProvider, "__init__", lambda self, **kw: None):
            provider = OllamaProvider.__new__(OllamaProvider)
            provider._client = mock_client
            provider._model = "qwen3.5:4b"
            provider._host = "http://localhost:11434"
            provider._keep_alive = "10m"

            messages = [Message(role=Role.USER, content="Hi")]
            result = list(
                provider.stream_chat(messages, system_prompt="Be helpful")
            )

            assert result == ["Hello", " world"]
            mock_client.chat.assert_called_once()

            call_kwargs = mock_client.chat.call_args.kwargs
            assert call_kwargs["model"] == "qwen3.5:4b"
            assert call_kwargs["stream"] is True
            assert call_kwargs["messages"][0]["role"] == "system"


def test_ollama_provider_name():
    with patch("epic.providers.ollama_provider.OllamaClient", create=True):
        from epic.providers.ollama_provider import OllamaProvider

        with patch.object(OllamaProvider, "__init__", lambda self, **kw: None):
            p = OllamaProvider.__new__(OllamaProvider)
            p._model = "qwen3.5:4b"
            p._host = "http://localhost:11434"
            p._keep_alive = "10m"
            p._client = MagicMock()
            assert p.name == "ollama"
            assert p.model_name == "qwen3.5:4b"
