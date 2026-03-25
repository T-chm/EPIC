"""Tests for epic.tools.interpreter."""

from __future__ import annotations

from epic.tools.interpreter import CodeInterpreter


def test_execute_captures_stdout():
    interp = CodeInterpreter()
    stdout, stderr, success = interp.execute("print('hello world')")
    assert stdout == "hello world"
    assert success is True


def test_execute_captures_error():
    interp = CodeInterpreter()
    stdout, stderr, success = interp.execute("1/0")
    assert success is False
    assert "ZeroDivision" in stderr


def test_execute_empty_code():
    interp = CodeInterpreter()
    stdout, stderr, success = interp.execute("")
    assert "No code provided" in stdout
    assert success is False


def test_execute_whitespace_only():
    interp = CodeInterpreter()
    stdout, stderr, success = interp.execute("   \n  ")
    assert "No code provided" in stdout


def test_inject_globals():
    interp = CodeInterpreter()
    interp.inject_globals({"my_var": 42})
    stdout, stderr, success = interp.execute("print(my_var)")
    assert stdout == "42"
    assert success is True


def test_state_persists_across_executions():
    interp = CodeInterpreter()
    interp.execute("x = 10")
    stdout, _, success = interp.execute("print(x + 5)")
    assert stdout == "15"
    assert success is True
