"""MAFFT integration for multiple sequence alignment."""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile

logger = logging.getLogger(__name__)


def run_mafft(
    sequences: dict[str, str],
    output_format: str = "clustal",
) -> str:
    """Run MAFFT on a set of sequences for multiple sequence alignment.

    Args:
        sequences: Dictionary with sequence names as keys and sequences as values.
        output_format: 'clustal' for CLUSTAL format or 'default' for FASTA format.

    Returns:
        The alignment result as a string.

    Raises:
        FileNotFoundError: If MAFFT is not installed.
        RuntimeError: If MAFFT execution fails.
    """
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".fasta", delete=False
    ) as temp_in:
        for name, seq in sequences.items():
            temp_in.write(f">{name}\n{seq}\n")
        temp_in_name = temp_in.name

    temp_out_handle, temp_out_name = tempfile.mkstemp(suffix=".aln")
    os.close(temp_out_handle)

    try:
        cmd = ["mafft", "--quiet"]
        if output_format == "clustal":
            cmd.append("--clustalout")
        cmd.append(temp_in_name)

        logger.info("Running MAFFT command: %s", " ".join(cmd))

        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return process.stdout.strip()

    except subprocess.CalledProcessError as e:
        msg = f"MAFFT execution failed: {e.stderr}"
        logger.error(msg)
        raise RuntimeError(msg) from e
    except FileNotFoundError:
        msg = (
            "MAFFT not found. Please install MAFFT and ensure it is in your PATH. "
            "See: https://mafft.cbrc.jp/alignment/software/"
        )
        logger.error(msg)
        raise FileNotFoundError(msg)
    finally:
        for path in (temp_in_name, temp_out_name):
            try:
                os.remove(path)
            except OSError:
                pass
