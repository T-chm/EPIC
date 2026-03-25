"""Tests for epic.config."""

from __future__ import annotations

import os

import pytest


def test_settings_defaults():
    """Settings should have sensible defaults."""
    from epic.config import Settings

    s = Settings()
    assert s.default_provider == "ollama"
    assert s.ollama_host == "http://localhost:11434"
    assert s.default_temperature == 0.2
    assert s.server_port == 8000


def test_settings_openai_key_compat(monkeypatch):
    """Plain OPENAI_API_KEY env var should be picked up."""
    monkeypatch.setenv("OPENAI_API_KEY", "sk-test-123")
    from epic.config import Settings

    s = Settings()
    assert s.openai_api_key == "sk-test-123"


def test_settings_strips_quotes(monkeypatch):
    """API keys with surrounding quotes should be stripped."""
    monkeypatch.setenv("OPENAI_API_KEY", '"sk-quoted"')
    from epic.config import Settings

    s = Settings()
    assert s.openai_api_key == "sk-quoted"
