"""Direct dispatch: bypass LLM for unambiguous QDB/PCR design requests.

When the user provides a clear sequence and asks for a standard design,
we can generate the Python code from a template instead of waiting for
the LLM to produce it. This saves 3-10 seconds of inference time.
The LLM still interprets the results in Turn 2.
"""

from __future__ import annotations

import re

# Match "design QDB/probes for [sequence]" with a DNA sequence of 20+ chars
_QDB_PATTERN = re.compile(
    r"(?:design|create|make)\s+(?:a\s+)?(?:nucleic\s+acid\s+)?"
    r"(?:quantum\s+dot\s+barcode|qdb)\s+(?:assay\s+)?(?:probes?\s+)?(?:for\s+)?"
    r"(?:(?:this|the|following)\s+)?(?:(?:target\s+)?sequence[:\s]+)?"
    r"(?:5['\u2019]-?)?\s*([ATGCUNRYWSMKHBVD\s\-]{10,})\s*(?:-?3['\u2019])?",
    re.IGNORECASE,
)

# Match "design PCR primers for [sequence]"
_PCR_PATTERN = re.compile(
    r"(?:design|create|make)\s+(?:a\s+)?(?:pcr\s+)?primers?\s+for\s+"
    r"(?:(?:this|the|following|template|target|sequence)\s+)*"  # skip filler words
    r"(?:sequence)?[:\s]*"
    r"(?:5['\u2019]-?)?\s*([ATGCUNRYWSMKHBVD\s\-]{20,})\s*(?:-?3['\u2019])?",
    re.IGNORECASE,
)

_QDB_CODE_TEMPLATE = '''\
target = "{sequence}"
seq = Seq(target)

# Step 1: Count nucleotides
length = len(seq)
print(f"Target length: {{length}} nucleotides")

# Step 2: Split into equal halves
half = length // 2
if length % 2 == 1:
    first_half = str(seq[:half])
    spacer = str(seq[half])
    second_half = str(seq[half+1:])
    print(f"Odd length: spacer nucleotide = {{spacer}}")
else:
    first_half = str(seq[:half])
    second_half = str(seq[half:])

print(f"First half  ({{len(first_half)}} nt): {{first_half}}")
print(f"Second half ({{len(second_half)}} nt): {{second_half}}")

# Step 3: Reverse complement each half
rc_first = str(Seq(first_half).reverse_complement())
rc_second = str(Seq(second_half).reverse_complement())
print(f"RC first half:  {{rc_first}}")
print(f"RC second half: {{rc_second}}")

# Step 4: Assign probes
capture_probe = rc_second
reporter_probe = rc_first

# Step 5: GC content
cap_gc = gc_fraction(Seq(capture_probe)) * 100
rep_gc = gc_fraction(Seq(reporter_probe)) * 100

# Step 6: Melting temperature
cap_tm = Tm_NN(Seq(capture_probe))
rep_tm = Tm_NN(Seq(reporter_probe))

# Step 7: Add modifications
final_capture = "/5AmMC6/" + capture_probe
final_reporter = reporter_probe + "/3Cy5Sp/"

# Step 8: Summary
print(f"\\n{{\'=\' * 60}}")
print(f"QDB ASSAY DESIGN RESULTS")
print(f"{{\'=\' * 60}}")
print(f"{{\'Probe\':<10}} {{\'Sequence\':<50}} {{\'GC%\':<8}} {{\'Tm (C)\':<8}}")
print(f"{{'-' * 76}}")
print(f"{{\'Capture\':<10}} {{final_capture:<50}} {{cap_gc:<8.2f}} {{cap_tm:<8.2f}}")
print(f"{{\'Reporter\':<10}} {{final_reporter:<50}} {{rep_gc:<8.2f}} {{rep_tm:<8.2f}}")
print(f"{{\'=\' * 60}}")

# Validation
issues = []
for name, gc, tm in [("Capture", cap_gc, cap_tm), ("Reporter", rep_gc, rep_tm)]:
    if not (35 <= gc <= 60):
        issues.append(f"{{name}} GC ({{gc:.1f}}%) outside 35-60% range")
    if not (55 <= tm <= 72):
        issues.append(f"{{name}} Tm ({{tm:.1f}}C) outside 55-72C range")

if issues:
    print(f"\\nWARNINGS:")
    for issue in issues:
        print(f"  - {{issue}}")
else:
    print(f"\\nAll probes pass validation criteria.")
'''

_PCR_CODE_TEMPLATE = '''\
template = "{sequence}"
seq = Seq(template)
length = len(seq)
print(f"Template length: {{length}} nucleotides")

# Design primers using primer3
result = primer3.design_primers(
    seq_args={{
        'SEQUENCE_ID': 'target',
        'SEQUENCE_TEMPLATE': template,
    }},
    global_args={{
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 30,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 72.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
    }}
)

num_pairs = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
print(f"\\nPrimer pairs found: {{num_pairs}}")

if num_pairs == 0:
    print("No suitable primers found. Try adjusting parameters or using a longer template.")
else:
    print(f"\\n{{'=' * 90}}")
    print(f"PCR PRIMER DESIGN RESULTS")
    print(f"{{'=' * 90}}")
    print(f"{{'Pair':<6}} {{'Dir':<8}} {{'Sequence':<35}} {{'Len':<5}} {{'Tm(C)':<8}} {{'GC%':<8}} {{'Product':<8}}")
    print(f"{{'-' * 90}}")

    for i in range(min(num_pairs, 5)):
        fwd_seq = result[f'PRIMER_LEFT_{{i}}_SEQUENCE']
        rev_seq = result[f'PRIMER_RIGHT_{{i}}_SEQUENCE']
        fwd_tm = result[f'PRIMER_LEFT_{{i}}_TM']
        rev_tm = result[f'PRIMER_RIGHT_{{i}}_TM']
        fwd_gc = result[f'PRIMER_LEFT_{{i}}_GC_PERCENT']
        rev_gc = result[f'PRIMER_RIGHT_{{i}}_GC_PERCENT']
        product_size = result[f'PRIMER_PAIR_{{i}}_PRODUCT_SIZE']

        print(f"{{i+1:<6}} {{'Fwd':<8}} {{fwd_seq:<35}} {{len(fwd_seq):<5}} {{fwd_tm:<8.2f}} {{fwd_gc:<8.2f}} {{product_size:<8}}")
        print(f"{{'':<6}} {{'Rev':<8}} {{rev_seq:<35}} {{len(rev_seq):<5}} {{rev_tm:<8.2f}} {{rev_gc:<8.2f}}")

        # Validation
        tm_diff = abs(fwd_tm - rev_tm)
        issues = []
        if tm_diff > 5:
            issues.append(f"Tm difference {{tm_diff:.1f}}C > 5C")
        if not fwd_seq[-1] in "GC":
            issues.append(f"Fwd 3\\'end is {{fwd_seq[-1]}}, not G/C")
        if not rev_seq[-1] in "GC":
            issues.append(f"Rev 3\\'end is {{rev_seq[-1]}}, not G/C")
        if issues:
            print(f"{{'':<6}} ** {{\", \".join(issues)}}")
        print()

    print(f"{{'=' * 90}}")
'''


def try_direct_dispatch(user_input: str) -> str | None:
    """Attempt to directly generate design code without LLM.

    Returns the Python code string if the request is unambiguous, else None.
    """
    # Check QDB pattern first
    m = _QDB_PATTERN.search(user_input)
    if m:
        raw = m.group(1).strip().replace(" ", "").replace("-", "").upper()
        if all(c in "ATGC" for c in raw) and len(raw) >= 20:
            return _QDB_CODE_TEMPLATE.format(sequence=raw)

    # Check PCR pattern
    m = _PCR_PATTERN.search(user_input)
    if m:
        raw = m.group(1).strip().replace(" ", "").replace("-", "").upper()
        if all(c in "ATGC" for c in raw) and len(raw) >= 20:
            return _PCR_CODE_TEMPLATE.format(sequence=raw)

    return None
