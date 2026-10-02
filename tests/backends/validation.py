"""Shared acceptance checks for opt-in backend validation drivers."""
from pathlib import Path
import re


def check_solver_log(path):
    text = Path(path).read_text(errors='replace')
    failures = re.findall(r'^.*(?:no convergence|HDF5-DIAG|\b(?:nan|[-+]?inf(?:inity)?)\b).*$'
                          , text, flags=re.IGNORECASE | re.MULTILINE)
    if failures:
        raise ValueError(f'Solver diagnostics invalidate {path}: ' + '; '.join(failures[:3]))
