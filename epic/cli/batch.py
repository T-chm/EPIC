"""Batch mode: process FASTA sequences through automated pipelines."""

from __future__ import annotations

import logging
from pathlib import Path

from rich.progress import Progress, SpinnerColumn, TextColumn

from epic.chat import ChatEngine
from epic.cli.display import console, print_error
from epic.config import get_settings
from epic.providers import get_provider

logger = logging.getLogger(__name__)

PIPELINE_PROMPTS = {
    "qdb": (
        "Design quantum dot barcode (QDB) probes for the following target sequence. "
        "Follow the default QDB design criteria. Generate and execute Python code to "
        "calculate GC content and melting temperature for the probes.\n\n"
        "Sequence name: {name}\n"
        "Sequence: {sequence}"
    ),
    "pcr": (
        "Design PCR primers for the following target sequence. "
        "Follow the default PCR primer design criteria. Use primer3-py to design "
        "the primers and report the results.\n\n"
        "Sequence name: {name}\n"
        "Sequence: {sequence}"
    ),
    "full": (
        "For the following target sequence, perform a complete assay design:\n"
        "1. Design QDB probes following default criteria\n"
        "2. Design PCR primers using primer3-py\n"
        "3. Report all results including GC content, Tm, and probe/primer sequences\n\n"
        "Sequence name: {name}\n"
        "Sequence: {sequence}"
    ),
}


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into (name, sequence) pairs."""
    try:
        from Bio import SeqIO

        records = list(SeqIO.parse(str(path), "fasta"))
        return [(r.id, str(r.seq)) for r in records]
    except ImportError:
        # Fallback: manual FASTA parsing
        sequences = []
        current_name = ""
        current_seq: list[str] = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_name:
                        sequences.append((current_name, "".join(current_seq)))
                    current_name = line[1:].split()[0]
                    current_seq = []
                elif line:
                    current_seq.append(line)
        if current_name:
            sequences.append((current_name, "".join(current_seq)))
        return sequences


def run_batch(
    input_path: Path,
    provider_name: str,
    model: str | None,
    pipeline: str,
    output_dir: Path,
) -> None:
    """Run batch processing on a FASTA file."""
    if pipeline not in PIPELINE_PROMPTS:
        print_error(f"Unknown pipeline: {pipeline!r}. Available: {', '.join(PIPELINE_PROMPTS)}")
        raise SystemExit(1)

    if not input_path.exists():
        print_error(f"Input file not found: {input_path}")
        raise SystemExit(1)

    settings = get_settings()

    try:
        provider = get_provider(provider_name, settings, model=model)
    except (ValueError, ImportError) as e:
        print_error(str(e))
        raise SystemExit(1)

    sequences = _parse_fasta(input_path)
    if not sequences:
        print_error(f"No sequences found in {input_path}")
        raise SystemExit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    prompt_template = PIPELINE_PROMPTS[pipeline]

    console.print(f"[bold cyan]EPIC Batch Mode[/]")
    console.print(f"  Provider: [green]{provider.name}[/]  Model: [green]{provider.model_name}[/]")
    console.print(f"  Pipeline: [green]{pipeline}[/]")
    console.print(f"  Sequences: [green]{len(sequences)}[/]")
    console.print(f"  Output: [green]{output_dir}[/]")
    console.print()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Processing sequences...", total=len(sequences))

        for i, (name, sequence) in enumerate(sequences):
            progress.update(task, description=f"Processing {name} ({i+1}/{len(sequences)})")

            engine = ChatEngine(
                provider=provider,
                auto_execute=True,
            )

            prompt = prompt_template.format(name=name, sequence=sequence)

            full_response = ""
            code_results: list[str] = []

            for event in engine.stream_response(prompt):
                if event.event_type == "token":
                    full_response += event.data
                elif event.event_type == "code_result":
                    code_results.append(event.data)

            # Save results
            result_file = output_dir / f"{name}_{pipeline}_result.txt"
            with open(result_file, "w") as f:
                f.write(f"# {name} - {pipeline.upper()} Design Results\n\n")
                f.write("## LLM Response\n\n")
                f.write(full_response)
                if code_results:
                    f.write("\n\n## Code Execution Results\n\n")
                    for j, result in enumerate(code_results, 1):
                        f.write(f"### Execution {j}\n{result}\n\n")

            progress.advance(task)

    console.print(f"\n[bold green]Done![/] Results saved to {output_dir}")
