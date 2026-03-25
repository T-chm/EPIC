"""Unit and integration tests for combined PCR + QDB design.

Test case:
  "Design PCR primers then the QDB assay for the following template sequence:
   5'-CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG
   CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG-3'"

This tests the full pipeline: both PCR primer design and QDB probe design
from a single user request, using the same target sequence.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_NN

# ═══════════════════════════════════════════════════════════════════════════════
# Ground Truth
# ═══════════════════════════════════════════════════════════════════════════════

TARGET_SEQUENCE = (
    "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG"
    "CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
)

USER_PROMPT = (
    "Design PCR primers then the QDB assay for the following template sequence: "
    f"5'-{TARGET_SEQUENCE}-3'"
)

# QDB ground truth (from test_qdb_design.py)
FIRST_HALF = "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG"
SECOND_HALF = "CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
CAPTURE_PROBE = "CCTTGTGTGGTCTGCATGAGTTTAGGCCTGAGTTGAGTCAGCACTG"
REPORTER_PROBE = "CTCATGGATTGTTGCAATTGTTTGGAGAAATCATCCAAATCTGCAG"
FINAL_CAPTURE = f"/5AmMC6/{CAPTURE_PROBE}"
FINAL_REPORTER = f"{REPORTER_PROBE}/3Cy5Sp/"
CAPTURE_GC = 52.17
REPORTER_GC = 39.13
CAPTURE_TM = 69.72
REPORTER_TM = 64.92

# PCR ground truth (primer3 with default parameters)
PCR_RESULT = primer3.design_primers(
    seq_args={"SEQUENCE_ID": "target", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
    global_args={
        "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 30,
        "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 55.0, "PRIMER_MAX_TM": 72.0,
        "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
        "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 12, "PRIMER_MAX_SELF_END": 8,
        "PRIMER_PAIR_MAX_COMPL_ANY": 12, "PRIMER_PAIR_MAX_COMPL_END": 8,
        "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100], [100, 125], [125, 150],
                                       [150, 175], [175, 200], [200, 225]],
    },
)
PCR_NUM_PAIRS = PCR_RESULT["PRIMER_PAIR_NUM_RETURNED"]
PCR_FWD_0 = PCR_RESULT["PRIMER_LEFT_0_SEQUENCE"]
PCR_REV_0 = PCR_RESULT["PRIMER_RIGHT_0_SEQUENCE"]
PCR_PRODUCT_0 = PCR_RESULT["PRIMER_PAIR_0_PRODUCT_SIZE"]


# ═══════════════════════════════════════════════════════════════════════════════
# Unit Tests — Verify Both Designs Are Correct Independently
# ═══════════════════════════════════════════════════════════════════════════════


class TestCombinedGroundTruth:
    """Verify that PCR and QDB ground truth values are consistent."""

    def test_target_length(self):
        assert len(TARGET_SEQUENCE) == 92

    def test_pcr_finds_primers(self):
        assert PCR_NUM_PAIRS >= 1

    def test_pcr_forward_primer_in_template(self):
        assert PCR_FWD_0 in TARGET_SEQUENCE

    def test_pcr_reverse_primer_rc_in_template(self):
        rev_rc = str(Seq(PCR_REV_0).reverse_complement())
        assert rev_rc in TARGET_SEQUENCE

    def test_pcr_product_size_reasonable(self):
        assert 75 <= PCR_PRODUCT_0 <= 225

    def test_pcr_primer_tm_in_range(self):
        fwd_tm = PCR_RESULT["PRIMER_LEFT_0_TM"]
        rev_tm = PCR_RESULT["PRIMER_RIGHT_0_TM"]
        assert 50 <= fwd_tm <= 72
        assert 50 <= rev_tm <= 72

    def test_qdb_probes_match_expected(self):
        computed_capture = str(Seq(SECOND_HALF).reverse_complement())
        computed_reporter = str(Seq(FIRST_HALF).reverse_complement())
        assert computed_capture == CAPTURE_PROBE
        assert computed_reporter == REPORTER_PROBE

    def test_qdb_gc_values(self):
        assert abs(gc_fraction(Seq(CAPTURE_PROBE)) * 100 - CAPTURE_GC) < 0.01
        assert abs(gc_fraction(Seq(REPORTER_PROBE)) * 100 - REPORTER_GC) < 0.01

    def test_qdb_tm_values(self):
        assert abs(Tm_NN(Seq(CAPTURE_PROBE)) - CAPTURE_TM) < 0.1
        assert abs(Tm_NN(Seq(REPORTER_PROBE)) - REPORTER_TM) < 0.1

    def test_qdb_modifications(self):
        assert FINAL_CAPTURE == "/5AmMC6/" + CAPTURE_PROBE
        assert FINAL_REPORTER == REPORTER_PROBE + "/3Cy5Sp/"


# ═══════════════════════════════════════════════════════════════════════════════
# Direct Dispatch Tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestCombinedDirectDispatch:
    """Verify the combined regex pattern and code template."""

    def test_dispatch_matches_combined_request(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        assert code is not None, "Combined request should match direct dispatch"

    def test_dispatch_code_has_both_designs(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        assert "primer3.design_primers" in code, "Should contain PCR primer design"
        assert "/5AmMC6/" in code, "Should contain QDB capture modification"
        assert "/3Cy5Sp/" in code, "Should contain QDB reporter modification"

    def test_dispatch_code_has_correct_sequence(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        assert TARGET_SEQUENCE in code

    def test_dispatch_alternate_phrasings(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        # "and" instead of "then"
        alt1 = f"Design PCR primers and the QDB assay for the following template sequence: {TARGET_SEQUENCE}"
        assert try_direct_dispatch(alt1) is not None

    def test_dispatch_does_not_match_partial(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        # Just PCR should NOT match the combined pattern
        pcr_only = f"Design PCR primers for {TARGET_SEQUENCE}"
        code = try_direct_dispatch(pcr_only)
        # Should match PCR pattern, not combined
        assert code is not None
        assert "/5AmMC6/" not in code  # No QDB in PCR-only dispatch

    def test_individual_dispatches_still_work(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        qdb = try_direct_dispatch(f"Design QDB assay for {TARGET_SEQUENCE}")
        assert qdb is not None
        assert "/5AmMC6/" in qdb
        assert "primer3" not in qdb

        pcr = try_direct_dispatch(f"Design PCR primers for the sequence: {TARGET_SEQUENCE}")
        assert pcr is not None
        assert "primer3" in pcr
        assert "/5AmMC6/" not in pcr


# ═══════════════════════════════════════════════════════════════════════════════
# Interpreter Execution Tests
# ═══════════════════════════════════════════════════════════════════════════════


def _make_interpreter():
    """Create interpreter with BioPython + primer3 pre-loaded."""
    from Bio import SeqUtils
    from Bio.SeqUtils import MeltingTemp as mt
    from epic.tools.interpreter import CodeInterpreter

    interp = CodeInterpreter()
    interp.inject_globals({
        "Seq": Seq, "gc_fraction": gc_fraction, "Tm_NN": Tm_NN,
        "SeqUtils": SeqUtils, "mt": mt, "MeltingTemp": mt,
        "primer3": primer3, "design_primers": primer3.design_primers,
    })
    return interp


def _parse_output(stdout: str) -> dict[str, str]:
    return dict(
        line.split("=", 1) for line in stdout.strip().split("\n") if "=" in line and not line.startswith("=")
    )


class TestCombinedInterpreterExecution:
    """Run the combined dispatch code in the interpreter and verify results."""

    def test_execution_succeeds(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, stderr, success = interp.execute(code)
        assert success, f"Combined code execution failed:\n{stderr}"

    def test_output_contains_pcr_results(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(code)

        assert "PCR PRIMER DESIGN RESULTS" in stdout
        assert "Primer pairs found:" in stdout
        assert "Fwd" in stdout
        assert "Rev" in stdout

    def test_output_contains_qdb_results(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(code)

        assert "QDB ASSAY DESIGN RESULTS" in stdout
        assert CAPTURE_PROBE in stdout
        assert REPORTER_PROBE in stdout
        assert "/5AmMC6/" in stdout
        assert "/3Cy5Sp/" in stdout

    def test_output_pcr_primer_sequences(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(code)

        # First primer pair should be present
        assert PCR_FWD_0 in stdout
        assert PCR_REV_0 in stdout

    def test_output_qdb_gc_and_tm(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(code)

        assert f"{CAPTURE_GC:.2f}" in stdout
        assert f"{REPORTER_GC:.2f}" in stdout
        assert f"{CAPTURE_TM:.2f}" in stdout
        assert f"{REPORTER_TM:.2f}" in stdout

    def test_output_template_length(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(USER_PROMPT)
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(code)

        assert "Template length: 92 nucleotides" in stdout


# ═══════════════════════════════════════════════════════════════════════════════
# ChatEngine Integration Tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestChatEngineCombinedDesign:
    """Full integration: ChatEngine handles the combined request via direct dispatch."""

    def test_chatengine_dispatches_combined(self):
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        # Turn 2: interpretation
        provider.stream_chat.return_value = iter(["Both designs look good."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        events = list(engine.stream_response(USER_PROMPT))

        event_types = [e.event_type for e in events]
        assert "code_detected" in event_types
        assert "code_result" in event_types
        assert "done" in event_types

    def test_chatengine_code_result_has_both(self):
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["Interpretation done."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        events = list(engine.stream_response(USER_PROMPT))

        code_result = next(e for e in events if e.event_type == "code_result")

        # PCR results
        assert "PCR PRIMER DESIGN RESULTS" in code_result.data
        assert "Primer pairs found:" in code_result.data

        # QDB results
        assert "QDB ASSAY DESIGN RESULTS" in code_result.data
        assert CAPTURE_PROBE in code_result.data
        assert REPORTER_PROBE in code_result.data

    def test_chatengine_sends_results_for_interpretation(self):
        """Turn 2 should receive the actual code output for interpretation."""
        from epic.chat import ChatEngine
        from epic.models import Role

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["The results are valid."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        list(engine.stream_response(USER_PROMPT))

        # History should be: user → assistant(direct) → user(results) → assistant(interp)
        assert len(engine.history) == 4
        assert engine.history[0].role == Role.USER
        assert engine.history[2].role == Role.USER
        assert "Code execution output" in engine.history[2].content
        assert "PCR PRIMER DESIGN RESULTS" in engine.history[2].content
        assert "QDB ASSAY DESIGN RESULTS" in engine.history[2].content

    def test_chatengine_llm_called_once_for_interpretation(self):
        """Direct dispatch skips Turn 1 LLM — only Turn 2 interpretation uses LLM."""
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["Done."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        list(engine.stream_response(USER_PROMPT))

        # Only 1 LLM call (Turn 2 interpretation), not 2
        assert provider.stream_chat.call_count == 1
