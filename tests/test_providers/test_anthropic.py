"""Tests for Anthropic provider."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from epic.models import Message, Role


def test_anthropic_uses_system_parameter():
    """System prompt must be passed as 'system' kwarg, NOT in messages."""
    mock_anthropic_module = MagicMock()
    mock_client = MagicMock()
    mock_anthropic_module.Anthropic.return_value = mock_client

    mock_stream_ctx = MagicMock()
    mock_stream_ctx.__enter__ = MagicMock(return_value=mock_stream_ctx)
    mock_stream_ctx.__exit__ = MagicMock(return_value=False)
    mock_stream_ctx.text_stream = iter(["Hello", " world"])
    mock_client.messages.stream.return_value = mock_stream_ctx

    with patch.dict("sys.modules", {"anthropic": mock_anthropic_module}):
        from importlib import reload
        import epic.providers.anthropic_provider as mod

        reload(mod)

        provider = mod.AnthropicProvider(api_key="test-key")
        messages = [Message(role=Role.USER, content="Hi")]
        result = list(provider.stream_chat(messages, system_prompt="Be helpful"))

        assert result == ["Hello", " world"]
        call_kwargs = mock_client.messages.stream.call_args.kwargs
        assert call_kwargs["system"] == "Be helpful"
        for msg in call_kwargs["messages"]:
            assert msg["role"] != "system"


def test_anthropic_filters_system_messages():
    """System messages in the list should be filtered out."""
    mock_anthropic_module = MagicMock()
    mock_client = MagicMock()
    mock_anthropic_module.Anthropic.return_value = mock_client

    mock_stream_ctx = MagicMock()
    mock_stream_ctx.__enter__ = MagicMock(return_value=mock_stream_ctx)
    mock_stream_ctx.__exit__ = MagicMock(return_value=False)
    mock_stream_ctx.text_stream = iter(["ok"])
    mock_client.messages.stream.return_value = mock_stream_ctx

    with patch.dict("sys.modules", {"anthropic": mock_anthropic_module}):
        from importlib import reload
        import epic.providers.anthropic_provider as mod

        reload(mod)

        provider = mod.AnthropicProvider(api_key="test-key")
        messages = [
            Message(role=Role.SYSTEM, content="old system msg"),
            Message(role=Role.USER, content="Hi"),
        ]
        list(provider.stream_chat(messages, system_prompt="New system"))

        call_kwargs = mock_client.messages.stream.call_args.kwargs
        assert len(call_kwargs["messages"]) == 1
        assert call_kwargs["messages"][0]["role"] == "user"


def test_anthropic_provider_name():
    mock_anthropic_module = MagicMock()
    with patch.dict("sys.modules", {"anthropic": mock_anthropic_module}):
        from importlib import reload
        import epic.providers.anthropic_provider as mod

        reload(mod)
        p = mod.AnthropicProvider(api_key="test")
        assert p.name == "anthropic"
