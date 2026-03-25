"""Shared test fixtures."""

from __future__ import annotations

import pytest

from epic.models import Message, Role


@pytest.fixture
def sample_messages() -> list[Message]:
    return [
        Message(role=Role.USER, content="Design QDB probes for ATGCATGCATGC"),
        Message(role=Role.ASSISTANT, content="I'll design probes for that sequence."),
        Message(role=Role.USER, content="What is the GC content?"),
    ]


@pytest.fixture
def sample_system_prompt() -> str:
    return "You are a helpful bioinformatics assistant."
