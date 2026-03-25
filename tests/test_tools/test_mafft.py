"""Tests for epic.tools.mafft."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from epic.tools.mafft import run_mafft


@patch("epic.tools.mafft.subprocess.run")
def test_run_mafft_clustal(mock_run):
    """run_mafft should call MAFFT with --clustalout for clustal format."""
    mock_run.return_value = MagicMock(stdout="CLUSTAL alignment\n\nseq1  ATGC\nseq2  ATGC")

    result = run_mafft({"seq1": "ATGC", "seq2": "ATGC"}, output_format="clustal")

    assert "CLUSTAL" in result
    call_args = mock_run.call_args
    cmd = call_args[0][0]
    assert "--clustalout" in cmd
    assert cmd[0] == "mafft"


@patch("epic.tools.mafft.subprocess.run")
def test_run_mafft_fasta(mock_run):
    """run_mafft with default format should not add --clustalout."""
    mock_run.return_value = MagicMock(stdout=">seq1\nATGC\n>seq2\nATGC")

    result = run_mafft({"seq1": "ATGC", "seq2": "ATGC"}, output_format="default")

    call_args = mock_run.call_args
    cmd = call_args[0][0]
    assert "--clustalout" not in cmd


@patch("epic.tools.mafft.subprocess.run", side_effect=FileNotFoundError)
def test_run_mafft_not_installed(mock_run):
    """Should raise FileNotFoundError if MAFFT is not installed."""
    with pytest.raises(FileNotFoundError, match="MAFFT not found"):
        run_mafft({"seq1": "ATGC"})
