"""Unit and integration tests for QDB assay design.

Primary test case:
  "Design QDB probes for this sequence:
   5'-CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG
   CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG-3'"

Ground truth computed independently with BioPython and verified by hand.
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_NN

# ═══════════════════════════════════════════════════════════════════════════════
# Ground Truth Constants
# ═══════════════════════════════════════════════════════════════════════════════

TARGET_SEQUENCE = (
    "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG"
    "CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
)
TARGET_LENGTH = 92

FIRST_HALF = "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG"
SECOND_HALF = "CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
HALF_LENGTH = 46

RC_FIRST_HALF = "CTCATGGATTGTTGCAATTGTTTGGAGAAATCATCCAAATCTGCAG"
RC_SECOND_HALF = "CCTTGTGTGGTCTGCATGAGTTTAGGCCTGAGTTGAGTCAGCACTG"

CAPTURE_PROBE = RC_SECOND_HALF
REPORTER_PROBE = RC_FIRST_HALF

CAPTURE_GC = 52.17
REPORTER_GC = 39.13
CAPTURE_TM = 69.72
REPORTER_TM = 64.92

FINAL_CAPTURE = f"/5AmMC6/{CAPTURE_PROBE}"
FINAL_REPORTER = f"{REPORTER_PROBE}/3Cy5Sp/"


# ═══════════════════════════════════════════════════════════════════════════════
# Unit Tests — Individual Design Steps
# ═══════════════════════════════════════════════════════════════════════════════


class TestSequenceInput:
    """Unit tests for Step 1: target sequence validation."""

    def test_target_length(self):
        assert len(TARGET_SEQUENCE) == 92

    def test_target_is_valid_dna(self):
        assert set(TARGET_SEQUENCE).issubset({"A", "T", "G", "C"})

    def test_target_matches_raw_input(self):
        """The sequence from the user prompt (with 5'- and -3' stripped)."""
        raw = "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAGCAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
        assert TARGET_SEQUENCE == raw
        assert len(raw) == 92


class TestSequenceSplitting:
    """Unit tests for Step 2: splitting into equal halves."""

    def test_even_split(self):
        half = len(TARGET_SEQUENCE) // 2
        first = TARGET_SEQUENCE[:half]
        second = TARGET_SEQUENCE[half:]
        assert first == FIRST_HALF
        assert second == SECOND_HALF

    def test_halves_equal_length(self):
        assert len(FIRST_HALF) == len(SECOND_HALF) == HALF_LENGTH

    def test_halves_reconstruct_target(self):
        assert FIRST_HALF + SECOND_HALF == TARGET_SEQUENCE

    def test_odd_length_sequence_uses_spacer(self):
        """Odd-length sequences should leave the middle nucleotide as spacer."""
        odd_seq = "ATGCATGCATGCA"  # 13 nt
        half = len(odd_seq) // 2  # 6
        first = odd_seq[:half]       # "ATGCAT"
        spacer = odd_seq[half]       # "G" (position 6, the middle)
        second = odd_seq[half + 1:]  # "CATGCA"
        assert len(first) == len(second) == 6
        assert spacer == "G"
        assert first + spacer + second == odd_seq

    def test_split_preserves_no_overlap(self):
        """No nucleotide should appear in both halves."""
        half = HALF_LENGTH
        # Indices [0..45] and [46..91] — no overlap
        assert TARGET_SEQUENCE[:half][-1] == "G"  # pos 45
        assert TARGET_SEQUENCE[half] == "C"  # pos 46


class TestReverseComplement:
    """Unit tests for Step 3: reverse complement computation."""

    def test_rc_first_half(self):
        assert str(Seq(FIRST_HALF).reverse_complement()) == RC_FIRST_HALF

    def test_rc_second_half(self):
        assert str(Seq(SECOND_HALF).reverse_complement()) == RC_SECOND_HALF

    def test_rc_is_involution(self):
        """RC of RC should return the original."""
        assert str(Seq(RC_FIRST_HALF).reverse_complement()) == FIRST_HALF
        assert str(Seq(RC_SECOND_HALF).reverse_complement()) == SECOND_HALF

    def test_rc_preserves_length(self):
        assert len(RC_FIRST_HALF) == len(FIRST_HALF)
        assert len(RC_SECOND_HALF) == len(SECOND_HALF)

    def test_rc_base_pairing(self):
        """Each base in the RC should be the complement of the reversed original."""
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        for orig, rc in [(FIRST_HALF, RC_FIRST_HALF), (SECOND_HALF, RC_SECOND_HALF)]:
            for i, base in enumerate(reversed(orig)):
                assert rc[i] == complement[base], (
                    f"Mismatch at RC position {i}: "
                    f"original(reversed)[{i}]={base}, RC[{i}]={rc[i]}"
                )

    def test_known_short_rc(self):
        """Sanity check: 5'-ATGC-3' → RC = 5'-GCAT-3'."""
        assert str(Seq("ATGC").reverse_complement()) == "GCAT"


class TestProbeAssignment:
    """Unit tests for Step 5: probe assignment."""

    def test_capture_is_rc_of_second_half(self):
        assert CAPTURE_PROBE == RC_SECOND_HALF

    def test_reporter_is_rc_of_first_half(self):
        assert REPORTER_PROBE == RC_FIRST_HALF

    def test_capture_and_reporter_are_different(self):
        assert CAPTURE_PROBE != REPORTER_PROBE

    def test_probes_bind_correct_halves(self):
        """RC of capture should bind second half, RC of reporter should bind first half."""
        assert str(Seq(CAPTURE_PROBE).reverse_complement()) == SECOND_HALF
        assert str(Seq(REPORTER_PROBE).reverse_complement()) == FIRST_HALF

    def test_probes_together_cover_full_target(self):
        """Reporter binds first half + capture binds second half = full target."""
        reporter_target = str(Seq(REPORTER_PROBE).reverse_complement())
        capture_target = str(Seq(CAPTURE_PROBE).reverse_complement())
        assert reporter_target + capture_target == TARGET_SEQUENCE


class TestGCContent:
    """Unit tests for Step 6: GC content validation."""

    def test_capture_gc_value(self):
        gc = gc_fraction(Seq(CAPTURE_PROBE)) * 100
        assert abs(gc - CAPTURE_GC) < 0.01

    def test_reporter_gc_value(self):
        gc = gc_fraction(Seq(REPORTER_PROBE)) * 100
        assert abs(gc - REPORTER_GC) < 0.01

    def test_capture_gc_in_range(self):
        gc = gc_fraction(Seq(CAPTURE_PROBE)) * 100
        assert 35 <= gc <= 60, f"Capture GC {gc:.2f}% outside 35-60%"

    def test_reporter_gc_in_range(self):
        gc = gc_fraction(Seq(REPORTER_PROBE)) * 100
        assert 35 <= gc <= 60, f"Reporter GC {gc:.2f}% outside 35-60%"

    def test_gc_manual_count_matches(self):
        """Cross-check gc_fraction against manual G+C count."""
        for probe, expected_gc in [(CAPTURE_PROBE, CAPTURE_GC), (REPORTER_PROBE, REPORTER_GC)]:
            g_count = probe.count("G")
            c_count = probe.count("C")
            manual_gc = (g_count + c_count) / len(probe) * 100
            assert abs(manual_gc - expected_gc) < 0.01

    def test_known_gc(self):
        """Sanity: 50% GC sequence."""
        assert abs(gc_fraction(Seq("ATGC")) * 100 - 50.0) < 0.01


class TestMeltingTemperature:
    """Unit tests for Step 7: melting temperature validation."""

    def test_capture_tm_value(self):
        tm = Tm_NN(Seq(CAPTURE_PROBE))
        assert abs(tm - CAPTURE_TM) < 0.1

    def test_reporter_tm_value(self):
        tm = Tm_NN(Seq(REPORTER_PROBE))
        assert abs(tm - REPORTER_TM) < 0.1

    def test_capture_tm_in_range(self):
        tm = Tm_NN(Seq(CAPTURE_PROBE))
        assert 55 <= tm <= 72, f"Capture Tm {tm:.2f}°C outside 55-72°C"

    def test_reporter_tm_in_range(self):
        tm = Tm_NN(Seq(REPORTER_PROBE))
        assert 55 <= tm <= 72, f"Reporter Tm {tm:.2f}°C outside 55-72°C"

    def test_tm_increases_with_gc(self):
        """Higher GC probe should generally have higher Tm."""
        # Capture has higher GC (52.17%) than reporter (39.13%)
        assert CAPTURE_TM > REPORTER_TM


class TestModifications:
    """Unit tests for Step 8: chemical modifications."""

    def test_capture_5prime_modification(self):
        assert FINAL_CAPTURE == "/5AmMC6/" + CAPTURE_PROBE

    def test_capture_starts_with_modification(self):
        assert FINAL_CAPTURE.startswith("/5AmMC6/")

    def test_capture_probe_sequence_after_mod(self):
        assert FINAL_CAPTURE[len("/5AmMC6/"):] == CAPTURE_PROBE

    def test_reporter_3prime_modification(self):
        assert FINAL_REPORTER == REPORTER_PROBE + "/3Cy5Sp/"

    def test_reporter_ends_with_modification(self):
        assert FINAL_REPORTER.endswith("/3Cy5Sp/")

    def test_reporter_probe_sequence_before_mod(self):
        assert FINAL_REPORTER[: -len("/3Cy5Sp/")] == REPORTER_PROBE

    def test_no_cross_modification(self):
        """Capture should NOT have /3Cy5Sp/, reporter should NOT have /5AmMC6/."""
        assert "/3Cy5Sp/" not in FINAL_CAPTURE
        assert "/5AmMC6/" not in FINAL_REPORTER


# ═══════════════════════════════════════════════════════════════════════════════
# Integration Tests — Full Pipeline
# ═══════════════════════════════════════════════════════════════════════════════


def _make_interpreter():
    """Create a CodeInterpreter with BioPython pre-loaded (same as ChatEngine)."""
    from Bio import SeqUtils
    from Bio.SeqUtils import MeltingTemp as mt

    from epic.tools.interpreter import CodeInterpreter

    interp = CodeInterpreter()
    interp.inject_globals({
        "Seq": Seq,
        "mt": mt,
        "MeltingTemp": mt,
        "Tm_NN": Tm_NN,
        "gc_fraction": gc_fraction,
        "SeqUtils": SeqUtils,
    })
    return interp


# The canonical QDB design script — what a correct LLM response should produce
QDB_DESIGN_CODE = f'''\
target = "{TARGET_SEQUENCE}"
seq = Seq(target)

# Step 1: Count
length = len(seq)
print(f"LENGTH={{length}}")

# Step 2: Split
half = length // 2
first_half = str(seq[:half])
second_half = str(seq[half:])
print(f"FIRST_HALF={{first_half}}")
print(f"SECOND_HALF={{second_half}}")

# Step 3: Reverse complement
rc_first = str(Seq(first_half).reverse_complement())
rc_second = str(Seq(second_half).reverse_complement())
print(f"RC_FIRST={{rc_first}}")
print(f"RC_SECOND={{rc_second}}")

# Step 5: Assign probes
capture = rc_second
reporter = rc_first
print(f"CAPTURE={{capture}}")
print(f"REPORTER={{reporter}}")

# Step 6: GC content
cap_gc = gc_fraction(Seq(capture)) * 100
rep_gc = gc_fraction(Seq(reporter)) * 100
print(f"CAP_GC={{cap_gc:.2f}}")
print(f"REP_GC={{rep_gc:.2f}}")

# Step 7: Melting temperature
cap_tm = Tm_NN(Seq(capture))
rep_tm = Tm_NN(Seq(reporter))
print(f"CAP_TM={{cap_tm:.2f}}")
print(f"REP_TM={{rep_tm:.2f}}")

# Step 8: Modifications
final_capture = "/5AmMC6/" + capture
final_reporter = reporter + "/3Cy5Sp/"
print(f"FINAL_CAP={{final_capture}}")
print(f"FINAL_REP={{final_reporter}}")
'''


def _parse_output(stdout: str) -> dict[str, str]:
    """Parse KEY=VALUE output lines into a dict."""
    return dict(line.split("=", 1) for line in stdout.strip().split("\n") if "=" in line)


class TestInterpreterFullPipeline:
    """Integration: run the full QDB design code in the interpreter and verify every value."""

    def test_execution_succeeds(self):
        interp = _make_interpreter()
        stdout, stderr, success = interp.execute(QDB_DESIGN_CODE)
        assert success, f"Execution failed:\nstderr: {stderr}"

    def test_length(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["LENGTH"] == str(TARGET_LENGTH)

    def test_split(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["FIRST_HALF"] == FIRST_HALF
        assert r["SECOND_HALF"] == SECOND_HALF

    def test_reverse_complements(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["RC_FIRST"] == RC_FIRST_HALF
        assert r["RC_SECOND"] == RC_SECOND_HALF

    def test_probe_assignment(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["CAPTURE"] == CAPTURE_PROBE
        assert r["REPORTER"] == REPORTER_PROBE

    def test_gc_content(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["CAP_GC"] == f"{CAPTURE_GC:.2f}"
        assert r["REP_GC"] == f"{REPORTER_GC:.2f}"

    def test_melting_temperature(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["CAP_TM"] == f"{CAPTURE_TM:.2f}"
        assert r["REP_TM"] == f"{REPORTER_TM:.2f}"

    def test_final_probes(self):
        interp = _make_interpreter()
        stdout, _, _ = interp.execute(QDB_DESIGN_CODE)
        r = _parse_output(stdout)
        assert r["FINAL_CAP"] == FINAL_CAPTURE
        assert r["FINAL_REP"] == FINAL_REPORTER


class TestInterpreterImportSanitization:
    """Integration: verify the interpreter handles bad LLM imports gracefully."""

    def test_wrong_import_path_still_works(self):
        """LLM generates 'from Bio.SeqUtils import Tm_NN' (wrong) — should still work."""
        interp = _make_interpreter()
        code = (
            "from Bio.SeqUtils import Tm_NN\n"
            "from Bio.SeqUtils.MeltingTemp import mt\n"
            "from Bio.Seq import Seq\n"
            "import primer3\n"
            "\n"
            f'seq = Seq("{TARGET_SEQUENCE}")\n'
            'print(f"LEN={len(seq)}")\n'
            'print(f"TM={Tm_NN(seq[:46]):.2f}")\n'
        )
        stdout, stderr, success = interp.execute(code)
        assert success, f"Failed with sanitized imports: {stderr}"
        r = _parse_output(stdout)
        assert r["LEN"] == "92"

    def test_no_import_at_all_still_works(self):
        """Code that uses pre-loaded symbols without any imports."""
        interp = _make_interpreter()
        code = f"""\
seq = Seq("{FIRST_HALF}")
gc = gc_fraction(seq) * 100
tm = Tm_NN(seq)
print(f"GC={{gc:.2f}}")
print(f"TM={{tm:.2f}}")
"""
        stdout, stderr, success = interp.execute(code)
        assert success, f"Failed without imports: {stderr}"
        r = _parse_output(stdout)
        assert float(r["GC"]) == pytest.approx(REPORTER_GC, abs=0.01)

    def test_non_bio_imports_preserved(self):
        """Imports for os, math, etc. should NOT be stripped."""
        interp = _make_interpreter()
        code = """\
import math
print(f"PI={math.pi:.2f}")
"""
        stdout, stderr, success = interp.execute(code)
        assert success, f"Non-bio import was incorrectly stripped: {stderr}"
        r = _parse_output(stdout)
        assert r["PI"] == "3.14"


class TestChatEngineQDBDesign:
    """Integration: verify ChatEngine correctly streams, executes, and triggers interpretation."""

    def test_chatengine_extracts_and_executes_code(self):
        """Mock provider returns a response containing a python code block.
        ChatEngine should extract it, execute it, and yield code_result events."""
        from epic.chat import ChatEngine

        llm_response = f"Here is the QDB design:\n```python\n{QDB_DESIGN_CODE}\n```"

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        # Turn 1: code response, Turn 2: interpretation
        provider.stream_chat.side_effect = [
            iter([llm_response]),
            iter(["The capture probe GC is 52.17% which is within range."]),
        ]

        engine = ChatEngine(
            provider=provider,
            system_prompt="test",
            auto_execute=True,
        )

        events = list(engine.stream_response(
            "Design QDB probes for: " + TARGET_SEQUENCE
        ))

        event_types = [e.event_type for e in events]
        assert "token" in event_types
        assert "code_detected" in event_types
        assert "code_result" in event_types
        assert "done" in event_types

        # Verify the code result contains correct values
        code_result = next(e for e in events if e.event_type == "code_result")
        assert f"CAPTURE={CAPTURE_PROBE}" in code_result.data
        assert f"REPORTER={REPORTER_PROBE}" in code_result.data
        assert f"FINAL_CAP={FINAL_CAPTURE}" in code_result.data
        assert f"FINAL_REP={FINAL_REPORTER}" in code_result.data

    def test_chatengine_two_turn_interpretation(self):
        """After code executes, ChatEngine should call the LLM a second time
        with the actual code output for interpretation."""
        from epic.chat import ChatEngine
        from epic.models import Role

        llm_response = f"Plan:\n```python\n{QDB_DESIGN_CODE}\n```"

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.side_effect = [
            iter([llm_response]),
            iter(["Interpretation: all probes meet criteria."]),
        ]

        engine = ChatEngine(
            provider=provider,
            system_prompt="test",
            auto_execute=True,
        )

        events = list(engine.stream_response("Design QDB probes"))

        # Provider should have been called TWICE (turn 1 + turn 2)
        assert provider.stream_chat.call_count == 2

        # The history should contain the code output as a user message (the followup)
        # History: user → assistant(code) → user(results) → assistant(interpretation)
        assert len(engine.history) == 4
        followup_msg = engine.history[2]
        assert followup_msg.role == Role.USER
        assert "Code execution output" in followup_msg.content
        assert "Interpret these results" in followup_msg.content

        # The interpretation tokens should appear in events
        token_texts = [e.data for e in events if e.event_type == "token"]
        full_text = "".join(token_texts)
        assert "Interpretation" in full_text

    def test_chatengine_no_second_turn_without_code(self):
        """If the LLM response has no code, there should be no second turn."""
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["Just a text response, no code."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        list(engine.stream_response("What is QDB?"))

        # Only one LLM call (no code → no interpretation turn)
        assert provider.stream_chat.call_count == 1

    def test_chatengine_no_second_turn_when_auto_execute_off(self):
        """With auto_execute=False, code is detected but not run, no second turn."""
        from epic.chat import ChatEngine

        llm_response = "```python\nprint('hello')\n```"

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter([llm_response])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)
        events = list(engine.stream_response("Write code"))

        assert provider.stream_chat.call_count == 1
        event_types = [e.event_type for e in events]
        assert "code_detected" in event_types
        assert "code_result" not in event_types

    def test_chatengine_history_with_two_turns(self):
        """After a two-turn design, history should have 4 entries:
        user → assistant(code) → user(results) → assistant(interpretation)."""
        from epic.chat import ChatEngine
        from epic.models import Role

        llm_response = f"Design:\n```python\nprint('done')\n```"

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.side_effect = [
            iter([llm_response]),
            iter(["The output looks correct."]),
        ]

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        list(engine.stream_response("Design probes"))

        assert len(engine.history) == 4
        assert engine.history[0].role == Role.USER
        assert engine.history[1].role == Role.ASSISTANT
        assert engine.history[2].role == Role.USER
        assert "Code execution output" in engine.history[2].content
        assert engine.history[3].role == Role.ASSISTANT
        assert engine.history[3].content == "The output looks correct."

    def test_chatengine_records_history_no_code(self):
        """Without code, history has just 2 entries."""
        from epic.chat import ChatEngine
        from epic.models import Role

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["Design complete."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=False)
        list(engine.stream_response("Design QDB probes"))

        assert len(engine.history) == 2
        assert engine.history[0].role == Role.USER
        assert engine.history[1].role == Role.ASSISTANT
        assert engine.history[1].content == "Design complete."


class TestSystemPromptQDB:
    """Verify the system prompt encodes correct QDB design instructions."""

    def test_prompt_says_len_not_manual_count(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "NEVER count nucleotides by hand" in SYSTEM_PROMPT
        assert "len()" in SYSTEM_PROMPT

    def test_prompt_has_workflow_steps(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "PLAN" in SYSTEM_PROMPT
        assert "CODE" in SYSTEM_PROMPT
        # INTERPRET is handled by the two-turn ChatEngine, not the prompt
        assert "STOP after the code block" in SYSTEM_PROMPT

    def test_prompt_has_correct_probe_assignment(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "Capture probe = reverse complement of second half" in SYSTEM_PROMPT
        assert "Reporter probe = reverse complement of first half" in SYSTEM_PROMPT

    def test_prompt_has_correct_modifications(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "/5AmMC6/" in SYSTEM_PROMPT
        assert "/3Cy5Sp/" in SYSTEM_PROMPT
        assert "5' end of capture probe" in SYSTEM_PROMPT
        assert "3' end of reporter probe" in SYSTEM_PROMPT

    def test_prompt_gc_range(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "35-60%" in SYSTEM_PROMPT

    def test_prompt_tm_range(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "55-72°C" in SYSTEM_PROMPT

    def test_prompt_no_docker_chat(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "docker_chat" not in SYSTEM_PROMPT

    def test_prompt_preloaded_symbols(self):
        from epic.prompts import SYSTEM_PROMPT
        for sym in ("Seq", "gc_fraction", "Tm_NN", "primer3", "run_mafft"):
            assert sym in SYSTEM_PROMPT
