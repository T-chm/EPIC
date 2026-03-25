"""Unit and integration tests for PCR primer design.

Test case:
  "Design PCR primers for the following template sequence:
   5'-CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG
   CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG-3'"
"""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

TARGET_SEQUENCE = (
    "CTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAG"
    "CAGTGCTGACTCAACTCAGGCCTAAACTCATGCAGACCACACAAGG"
)


class TestPrimer3API:
    """Unit tests: verify primer3.design_primers works with this template."""

    def test_primer3_returns_results(self):
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 30,
                "PRIMER_OPT_TM": 60.0,
                "PRIMER_MIN_TM": 55.0,
                "PRIMER_MAX_TM": 72.0,
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        assert result["PRIMER_PAIR_NUM_RETURNED"] >= 1

    def test_primer_sequences_exist(self):
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 30,
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        assert "PRIMER_LEFT_0_SEQUENCE" in result
        assert "PRIMER_RIGHT_0_SEQUENCE" in result

    def test_primer_tm_in_range(self):
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 30,
                "PRIMER_MIN_TM": 55.0,
                "PRIMER_MAX_TM": 72.0,
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        fwd_tm = result["PRIMER_LEFT_0_TM"]
        rev_tm = result["PRIMER_RIGHT_0_TM"]
        assert 50 <= fwd_tm <= 72, f"Fwd Tm {fwd_tm} out of range"
        assert 50 <= rev_tm <= 72, f"Rev Tm {rev_tm} out of range"

    def test_primer_length_in_range(self):
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 30,
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        fwd = result["PRIMER_LEFT_0_SEQUENCE"]
        rev = result["PRIMER_RIGHT_0_SEQUENCE"]
        assert 18 <= len(fwd) <= 30, f"Fwd length {len(fwd)} out of range"
        assert 18 <= len(rev) <= 30, f"Rev length {len(rev)} out of range"

    def test_product_size_in_range(self):
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_OPT_SIZE": 20,
                "PRIMER_MIN_SIZE": 18,
                "PRIMER_MAX_SIZE": 30,
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        size = result["PRIMER_PAIR_0_PRODUCT_SIZE"]
        assert 75 <= size <= 100, f"Product size {size} out of range"

    def test_primers_found_on_template(self):
        """Forward primer should appear in template, reverse primer RC should appear."""
        result = primer3.design_primers(
            seq_args={"SEQUENCE_ID": "test", "SEQUENCE_TEMPLATE": TARGET_SEQUENCE},
            global_args={
                "PRIMER_PRODUCT_SIZE_RANGE": [[75, 100]],
            },
        )
        fwd = result["PRIMER_LEFT_0_SEQUENCE"]
        rev = result["PRIMER_RIGHT_0_SEQUENCE"]
        rev_rc = str(Seq(rev).reverse_complement())

        assert fwd in TARGET_SEQUENCE, f"Fwd primer {fwd} not found in template"
        assert rev_rc in TARGET_SEQUENCE, f"Rev primer RC {rev_rc} not found in template"


class TestDirectDispatchPCR:
    """Integration: verify direct dispatch generates correct PCR code."""

    def test_dispatch_matches_pcr_request(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        code = try_direct_dispatch(
            "Design PCR primers for the following template sequence: "
            f"5'-{TARGET_SEQUENCE}-3'"
        )
        assert code is not None
        assert "primer3.design_primers" in code
        assert TARGET_SEQUENCE in code

    def test_dispatch_does_not_match_vague_request(self):
        from epic.tools.direct_dispatch import try_direct_dispatch

        assert try_direct_dispatch("What are PCR primers?") is None
        assert try_direct_dispatch("Design primers") is None

    def test_dispatch_code_executes_successfully(self):
        from epic.tools.direct_dispatch import try_direct_dispatch
        from epic.tools.interpreter import CodeInterpreter

        code = try_direct_dispatch(
            f"Design PCR primers for the following template sequence: {TARGET_SEQUENCE}"
        )
        assert code is not None

        interp = CodeInterpreter()
        interp.inject_globals({
            "Seq": Seq,
            "gc_fraction": gc_fraction,
            "primer3": primer3,
            "design_primers": primer3.design_primers,
        })

        stdout, stderr, success = interp.execute(code)
        assert success, f"PCR code execution failed: {stderr}"
        assert "Primer pairs found:" in stdout
        assert "PCR PRIMER DESIGN RESULTS" in stdout
        assert "Fwd" in stdout
        assert "Rev" in stdout


class TestChatEnginePCR:
    """Integration: verify ChatEngine handles PCR design via direct dispatch."""

    def test_chatengine_direct_dispatch_pcr(self):
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        # Turn 2: interpretation
        provider.stream_chat.return_value = iter(["The primers look good."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        events = list(engine.stream_response(
            f"Design PCR primers for the following template sequence: {TARGET_SEQUENCE}"
        ))

        event_types = [e.event_type for e in events]
        assert "code_detected" in event_types
        assert "code_result" in event_types

        # The code result should contain primer data
        code_result = next(e for e in events if e.event_type == "code_result")
        assert "Primer pairs found:" in code_result.data
        assert "Fwd" in code_result.data

    def test_chatengine_pcr_falls_back_to_llm(self):
        """A non-obvious PCR request should go through the LLM, not direct dispatch."""
        from epic.chat import ChatEngine

        provider = MagicMock()
        provider.name = "mock"
        provider.model_name = "mock-model"
        provider.stream_chat.return_value = iter(["I can help with primer design."])

        engine = ChatEngine(provider=provider, system_prompt="test", auto_execute=True)
        list(engine.stream_response("How do I design primers with specific constraints?"))

        # Should have gone through LLM (not direct dispatch)
        assert provider.stream_chat.call_count == 1


class TestSystemPromptPCR:
    """Verify the system prompt has correct PCR instructions."""

    def test_prompt_has_correct_primer3_api(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "primer3.design_primers()" in SYSTEM_PROMPT
        # Should NOT reference the old bindings API
        assert "primer3.bindings.design_primers" not in SYSTEM_PROMPT

    def test_prompt_has_primer_parameters(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "PRIMER_OPT_SIZE" in SYSTEM_PROMPT
        assert "PRIMER_PRODUCT_SIZE_RANGE" in SYSTEM_PROMPT
        assert "PRIMER_LEFT_" in SYSTEM_PROMPT
        assert "PRIMER_PAIR_NUM_RETURNED" in SYSTEM_PROMPT

    def test_prompt_has_primer_criteria(self):
        from epic.prompts import SYSTEM_PROMPT
        assert "50-65" in SYSTEM_PROMPT  # Tm range
        assert "40-60" in SYSTEM_PROMPT  # GC range
