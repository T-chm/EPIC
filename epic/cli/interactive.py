"""Interactive chat mode for the CLI."""

from __future__ import annotations

import logging

from epic.chat import ChatEngine
from epic.cli.display import (
    console,
    finish_assistant_response,
    print_assistant_token,
    print_banner,
    print_code_block,
    print_code_result,
    print_error,
    print_help,
    print_provider_info,
    print_user_message,
)
from epic.config import get_settings
from epic.models import StreamEvent
from epic.providers import get_provider

logger = logging.getLogger(__name__)


def _handle_slash_command(
    command: str,
    engine: ChatEngine,
) -> bool:
    """Handle slash commands. Returns True if command was handled."""
    parts = command.strip().split(maxsplit=1)
    cmd = parts[0].lower()
    arg = parts[1] if len(parts) > 1 else ""

    settings = get_settings()

    if cmd == "/help":
        print_help()
        return True

    elif cmd == "/exit":
        console.print("[dim]Goodbye![/]")
        raise SystemExit(0)

    elif cmd == "/clear":
        engine.clear_history()
        console.print("[dim]Conversation history cleared.[/]")
        return True

    elif cmd == "/history":
        history = engine.history
        if not history:
            console.print("[dim]No conversation history.[/]")
        else:
            for msg in history:
                role_color = "blue" if msg.role.value == "user" else "green"
                preview = msg.content[:120] + ("..." if len(msg.content) > 120 else "")
                console.print(f"[bold {role_color}]{msg.role.value}:[/] {preview}")
        console.print()
        return True

    elif cmd == "/provider":
        if not arg:
            print_error("Usage: /provider <name>  (openai, anthropic, ollama)")
            return True
        try:
            new_provider = get_provider(arg.strip(), settings)
            engine.switch_provider(new_provider)
            print_provider_info(new_provider.name, new_provider.model_name)
        except (ValueError, ImportError) as e:
            print_error(str(e))
        return True

    elif cmd == "/model":
        if not arg:
            print_error("Usage: /model <name>")
            return True
        try:
            new_provider = get_provider(
                engine.provider.name, settings, model=arg.strip()
            )
            engine.switch_provider(new_provider)
            print_provider_info(new_provider.name, new_provider.model_name)
        except (ValueError, ImportError) as e:
            print_error(str(e))
        return True

    return False


def _process_event(event: StreamEvent, no_execute: bool) -> None:
    """Process a single stream event for display."""
    if event.event_type == "token":
        print_assistant_token(event.data)
    elif event.event_type == "code_detected":
        if no_execute:
            print_code_block(event.data)
    elif event.event_type == "code_executing":
        pass  # Displayed when result comes
    elif event.event_type == "code_result":
        is_error = event.data.startswith("Error:")
        print_code_result(event.data, is_error=is_error)
    elif event.event_type == "done":
        finish_assistant_response()
    elif event.event_type == "error":
        print_error(event.data)


def run_interactive(
    provider_name: str,
    model: str | None,
    temperature: float,
    no_execute: bool,
) -> None:
    """Run the interactive chat loop."""
    settings = get_settings()

    try:
        provider = get_provider(provider_name, settings, model=model)
    except (ValueError, ImportError) as e:
        print_error(str(e))
        raise SystemExit(1)

    engine = ChatEngine(
        provider=provider,
        auto_execute=not no_execute,
        temperature=temperature,
    )

    print_banner()
    print_provider_info(provider.name, provider.model_name)
    print_help()

    while True:
        try:
            user_input = console.input("[bold blue]You:[/] ").strip()
        except (EOFError, KeyboardInterrupt):
            console.print("\n[dim]Goodbye![/]")
            break

        if not user_input:
            continue

        # Slash commands
        if user_input.startswith("/"):
            if _handle_slash_command(user_input, engine):
                continue

        # Direct code execution
        if user_input.startswith("!code"):
            code_str = user_input[len("!code"):].strip()
            if code_str:
                result = engine.execute_code(code_str)
                print_code_block(code_str)
                print_code_result(
                    result.stdout if result.success else result.stderr,
                    is_error=not result.success,
                )
            else:
                print_error("No code provided. Usage: !code <python code>")
            continue

        # Normal chat with streaming
        print_user_message("")  # Blank line before assistant response
        console.print("[bold green]Assistant:[/] ", end="")

        # If auto-execute is off, prompt user for each code block
        if no_execute:
            for event in engine.stream_response(user_input):
                _process_event(event, no_execute=True)
        else:
            for event in engine.stream_response(user_input):
                _process_event(event, no_execute=False)
