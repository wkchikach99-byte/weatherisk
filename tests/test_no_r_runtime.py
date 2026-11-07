"""Guard: the weatherisk package must not call R at runtime.

weatherisk is a pure-Python package.  R (r_code/ and the generator scripts
under tests/) is used only to produce the frozen reference data in
tests/reference_data/.  This test fails if any file shipped under
weatherisk/ mentions ``Rscript`` or ``rpy2``, the two ways the package could
reach an R installation.
"""

from __future__ import annotations

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent.parent / "weatherisk"
FORBIDDEN = (b"Rscript", b"rpy2")


def find_r_references(root: Path) -> list[str]:
    """Return 'path: token' for every file under root that contains a forbidden token."""
    hits = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        data = path.read_bytes()
        hits.extend(
            f"{path.relative_to(root)}: {token.decode()}"
            for token in FORBIDDEN
            if token in data
        )
    return hits


def test_package_contains_no_r_calls():
    n_files = sum(
        1 for p in PACKAGE_ROOT.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    )
    assert n_files > 0, f"no files found under {PACKAGE_ROOT}"
    hits = find_r_references(PACKAGE_ROOT)
    assert not hits, "weatherisk/ must not call R:\n" + "\n".join(hits)


def test_guard_detects_an_rscript_call(tmp_path):
    (tmp_path / "clean.py").write_text("import numpy\n")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "bridge.py").write_text(
        "import subprocess\nsubprocess.run(['Rscript', 'x.R'])\n"
    )
    assert find_r_references(tmp_path) == [str(Path("sub") / "bridge.py") + ": Rscript"]
