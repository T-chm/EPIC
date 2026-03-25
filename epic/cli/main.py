"""EPIC CLI — main entry point using Typer."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(
    name="epic",
    help="EPIC: Engineering Probes via Instructing a Chatbot — CLI tool for nucleic acid assay design.",
    no_args_is_help=True,
)


@app.command()
def chat(
    provider: str = typer.Option("ollama", "--provider", "-p", help="LLM provider (openai, anthropic, ollama)"),
    model: Optional[str] = typer.Option(None, "--model", "-m", help="Model name (uses provider default if omitted)"),
    temperature: float = typer.Option(0.2, "--temperature", "-t", help="Sampling temperature"),
    no_execute: bool = typer.Option(False, "--no-execute", help="Disable automatic code execution"),
    log_level: str = typer.Option("WARNING", "--log-level", help="Logging level"),
) -> None:
    """Interactive chat mode with streaming LLM responses."""
    logging.basicConfig(level=getattr(logging, log_level.upper(), logging.WARNING))
    from epic.cli.interactive import run_interactive

    run_interactive(provider, model, temperature, no_execute)


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Input FASTA file", exists=True),
    provider: str = typer.Option("ollama", "--provider", "-p", help="LLM provider"),
    model: Optional[str] = typer.Option(None, "--model", "-m", help="Model name"),
    pipeline: str = typer.Option("qdb", "--pipeline", help="Pipeline: qdb, pcr, or full"),
    output_dir: Path = typer.Option(Path("./results"), "--output", "-o", help="Output directory"),
    log_level: str = typer.Option("WARNING", "--log-level", help="Logging level"),
) -> None:
    """Batch mode: process FASTA sequences through automated design pipelines."""
    logging.basicConfig(level=getattr(logging, log_level.upper(), logging.WARNING))
    from epic.cli.batch import run_batch

    run_batch(input_file, provider, model, pipeline, output_dir)


@app.command()
def serve(
    host: str = typer.Option("0.0.0.0", "--host", "-h", help="Server host"),
    port: int = typer.Option(8000, "--port", help="Server port"),
    log_level: str = typer.Option("INFO", "--log-level", help="Logging level"),
) -> None:
    """Start the FastAPI web server."""
    logging.basicConfig(level=getattr(logging, log_level.upper(), logging.INFO))
    import uvicorn

    uvicorn.run("epic.api.app:create_app", host=host, port=port, factory=True)


@app.command()
def providers(
    log_level: str = typer.Option("WARNING", "--log-level", help="Logging level"),
) -> None:
    """List available LLM providers and their models."""
    logging.basicConfig(level=getattr(logging, log_level.upper(), logging.WARNING))
    from epic.cli.display import print_providers_table
    from epic.config import get_settings

    settings = get_settings()
    results = []

    # Try each provider
    for name in ("openai", "anthropic", "ollama"):
        try:
            from epic.providers import get_provider

            p = get_provider(name, settings)
            models = p.list_models()
            results.append({"name": name, "models": models})
        except (ValueError, ImportError, Exception) as e:
            results.append({"name": f"{name} (unavailable)", "models": [str(e)]})

    print_providers_table(results)


if __name__ == "__main__":
    app()
