"""Tests for OpenAI provider."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from epic.models import Message, Role


def test_openai_provider_formats_system_prompt():
    """System prompt should become the first message with role='system'."""
    mock_openai_module = MagicMock()
    mock_client = MagicMock()
    mock_openai_module.OpenAI.return_value = mock_client

    chunk = MagicMock()
    chunk.choices = [MagicMock()]
    chunk.choices[0].delta = MagicMock()
    chunk.choices[0].delta.content = "Hello"
    mock_client.chat.completions.create.return_value = [chunk]

    with patch.dict("sys.modules", {"openai": mock_openai_module}):
        from importlib import reload
        import epic.providers.openai_provider as mod

        reload(mod)

        provider = mod.OpenAIProvider(api_key="test-key", model="gpt-4o")
        messages = [Message(role=Role.USER, content="Hi")]
        result = list(provider.stream_chat(messages, system_prompt="Be helpful"))

        assert result == ["Hello"]
        call_kwargs = mock_client.chat.completions.create.call_args.kwargs
        assert call_kwargs["messages"][0] == {"role": "system", "content": "Be helpful"}
        assert call_kwargs["messages"][1] == {"role": "user", "content": "Hi"}
        assert call_kwargs["stream"] is True


def test_openai_provider_name():
    """Provider name should be 'openai'."""
    mock_openai_module = MagicMock()
    with patch.dict("sys.modules", {"openai": mock_openai_module}):
        from importlib import reload
        import epic.providers.openai_provider as mod

        reload(mod)
        p = mod.OpenAIProvider(api_key="test", model="gpt-4o")
        assert p.name == "openai"
        assert p.model_name == "gpt-4o"
