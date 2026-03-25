"""Tests for epic.chat.ChatEngine."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from epic.chat import ChatEngine
from epic.models import Message, Role


def _make_mock_provider(responses: list[str], second_turn: list[str] | None = None) -> MagicMock:
    """Create a mock provider that yields the given tokens.

    Args:
        responses: Tokens for the first LLM call.
        second_turn: If provided, tokens for the second LLM call (interpretation turn).
    """
    provider = MagicMock()
    provider.name = "mock"
    provider.model_name = "mock-model"
    if second_turn is not None:
        provider.stream_chat.side_effect = [iter(responses), iter(second_turn)]
    else:
        provider.stream_chat.return_value = iter(responses)
    return provider


def test_stream_response_yields_tokens():
    """ChatEngine should yield token events for each chunk."""
    provider = _make_mock_provider(["Hello", " world"])
    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)

    events = list(engine.stream_response("Hi"))
    token_events = [e for e in events if e.event_type == "token"]
    done_events = [e for e in events if e.event_type == "done"]

    assert len(token_events) == 2
    assert token_events[0].data == "Hello"
    assert token_events[1].data == " world"
    assert len(done_events) == 1
    assert done_events[0].data == "Hello world"


def test_history_management():
    """Messages should be appended to history."""
    provider = _make_mock_provider(["Response"])
    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)

    list(engine.stream_response("Hello"))

    assert len(engine.history) == 2
    assert engine.history[0].role == Role.USER
    assert engine.history[0].content == "Hello"
    assert engine.history[1].role == Role.ASSISTANT
    assert engine.history[1].content == "Response"


def test_clear_history():
    provider = _make_mock_provider(["Response"])
    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)
    list(engine.stream_response("Hello"))
    assert len(engine.history) == 2

    engine.clear_history()
    assert len(engine.history) == 0


def test_switch_provider():
    provider1 = _make_mock_provider(["From provider 1"])
    provider2 = _make_mock_provider(["From provider 2"])

    engine = ChatEngine(provider=provider1, system_prompt="test", auto_execute=False)
    assert engine.provider.name == "mock"

    engine.switch_provider(provider2)
    assert engine.provider is provider2


def test_code_extraction_and_execution():
    """ChatEngine should detect and execute code blocks when auto_execute=True."""
    code_response = "Here is code:\n```python\nprint('hello from code')\n```"
    provider = _make_mock_provider(
        [code_response],
        second_turn=["The code printed hello."],
    )
    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)

    events = list(engine.stream_response("Write code"))
    event_types = [e.event_type for e in events]

    assert "code_detected" in event_types
    assert "code_result" in event_types

    code_result = [e for e in events if e.event_type == "code_result"][0]
    assert "hello from code" in code_result.data


def test_code_not_executed_when_disabled():
    """With auto_execute=False, code should be detected but NOT executed."""
    code_response = "```python\nprint('should not run')\n```"
    provider = _make_mock_provider([code_response])
    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)

    events = list(engine.stream_response("Write code"))
    event_types = [e.event_type for e in events]

    # Code blocks are still detected, but not executed
    assert "code_detected" in event_types
    assert "code_executing" not in event_types
    assert "code_result" not in event_types


def test_error_handling():
    """Provider errors should yield error events."""
    provider = MagicMock()
    provider.name = "mock"
    provider.model_name = "mock-model"
    provider.stream_chat.side_effect = RuntimeError("API down")

    engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)
    events = list(engine.stream_response("Hi"))

    error_events = [e for e in events if e.event_type == "error"]
    assert len(error_events) == 1
    assert "API down" in error_events[0].data
