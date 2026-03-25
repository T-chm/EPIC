"""Rich console rendering helpers for the CLI."""

from __future__ import annotations

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table

console = Console()


def print_banner() -> None:
    """Print the EPIC startup banner."""
    console.print(
        Panel(
            "[bold cyan]EPIC[/] — Engineering Probes via Instructing a Chatbot\n"
            "[dim]Nucleic acid diagnostic assay design assistant[/]",
            border_style="cyan",
        )
    )


def print_provider_info(provider_name: str, model: str) -> None:
    """Print current provider/model info."""
    console.print(f"  Provider: [bold green]{provider_name}[/]  Model: [bold green]{model}[/]")
    console.print()


def print_user_message(message: str) -> None:
    """Print a user message."""
    console.print(f"[bold blue]You:[/] {message}")
    console.print()


def print_assistant_token(token: str) -> None:
    """Print a streaming token without newline."""
    console.print(token, end="", highlight=False)


def finish_assistant_response() -> None:
    """End the assistant response line."""
    console.print()
    console.print()


def print_code_block(code: str) -> None:
    """Print a detected code block with syntax highlighting."""
    console.print()
    console.print(
        Panel(
            Syntax(code, "python", theme="monokai", line_numbers=True),
            title="[bold yellow]Detected Code[/]",
            border_style="yellow",
        )
    )


def print_code_result(result: str, is_error: bool = False) -> None:
    """Print code execution result."""
    style = "red" if is_error else "green"
    title = "Execution Error" if is_error else "Execution Result"
    console.print(
        Panel(result, title=f"[bold {style}]{title}[/]", border_style=style)
    )
    console.print()


def print_error(message: str) -> None:
    """Print an error message."""
    console.print(f"[bold red]Error:[/] {message}")


def print_help() -> None:
    """Print interactive mode help."""
    table = Table(title="Commands", show_header=True)
    table.add_column("Command", style="cyan")
    table.add_column("Description")
    table.add_row("/provider <name>", "Switch LLM provider (openai, anthropic, ollama)")
    table.add_row("/model <name>", "Switch model")
    table.add_row("/history", "Show conversation history")
    table.add_row("/clear", "Clear conversation history")
    table.add_row("/help", "Show this help")
    table.add_row("/exit", "Quit")
    table.add_row("!code <python>", "Execute Python code directly")
    console.print(table)
    console.print()


def print_providers_table(providers: list[dict]) -> None:
    """Print a table of available providers and models."""
    table = Table(title="Available Providers", show_header=True)
    table.add_column("Provider", style="cyan")
    table.add_column("Models")
    for p in providers:
        table.add_row(p["name"], ", ".join(p["models"]))
    console.print(table)
