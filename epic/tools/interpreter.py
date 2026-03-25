"""Python code interpreter for executing LLM-generated code."""

from __future__ import annotations

import code
import io
import logging
import re
import sys
from typing import Any

logger = logging.getLogger(__name__)

# Import lines for modules that are already pre-loaded in the interpreter.
# LLMs frequently generate wrong import paths for these, so we strip them
# and rely on the pre-injected namespace instead.
_PRELOADED_IMPORT_RE = re.compile(
    r"^\s*(?:"
    r"from\s+(?:Bio|Bio\.Seq|Bio\.SeqUtils|Bio\.SeqUtils\.MeltingTemp)\s+import\s+[^\n]+"
    r"|import\s+(?:Bio|primer3)(?:\s*$|[^\S\n]+[^\n]*)"
    r"|from\s+epic\.tools\.mafft\s+import\s+[^\n]+"
    r")$",
    re.MULTILINE,
)


class CodeInterpreter:
    """Sandboxed Python interpreter that captures stdout/stderr."""

    def __init__(self) -> None:
        self._interp = code.InteractiveInterpreter()
        self._interp.locals["__name__"] = "__main__"

    def inject_globals(self, namespace: dict[str, Any]) -> None:
        """Add functions/modules to the interpreter's namespace."""
        self._interp.locals.update(namespace)

    def _sanitize_imports(self, code_str: str) -> str:
        """Strip Bio/primer3 import lines — these are already pre-loaded.

        LLMs often generate incorrect import paths (e.g. `from Bio.SeqUtils import Tm_NN`
        instead of `from Bio.SeqUtils.MeltingTemp import Tm_NN`). Since all relevant
        symbols are already injected into the namespace, we strip these imports entirely
        to avoid ImportErrors.
        """
        return _PRELOADED_IMPORT_RE.sub("# (import handled by runtime)", code_str)

    def execute(self, code_str: str) -> tuple[str, str, bool]:
        """Execute code and return (stdout, stderr, success).

        Returns:
            Tuple of (stdout output, stderr output, success boolean).
        """
        if not code_str.strip():
            return ("No code provided to execute.", "", False)

        code_str = self._sanitize_imports(code_str)
        logger.info("Executing code: %s...", code_str[:200])

        old_stdout, old_stderr = sys.stdout, sys.stderr
        buf_out = io.StringIO()
        buf_err = io.StringIO()
        sys.stdout = buf_out
        sys.stderr = buf_err

        try:
            self._interp.runcode(code_str)
            stdout = buf_out.getvalue()
            stderr = buf_err.getvalue()
            success = not bool(stderr.strip())
            return (
                stdout.strip() if stdout else "",
                stderr.strip() if stderr else "",
                success,
            )
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
