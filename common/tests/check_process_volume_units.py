from __future__ import annotations
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "bin" / "slimmc"

CASES = [
    ("0.100", "0.001"),
    ("0.100 L", "0.001 L"),
    ("0.100 l", "0.001 l"),
    ("100 mL", "1 mL"),
    ("100 ml", "1 ml"),
    ("100 ML", "1 ML"),
]

BASE = '''desc "volume unit contract"
param output_dir "results/main"
param kmc_volume 1e-19
param t_end 0.3
species A 0.2
species B 0
feed F A 0.2
rate k const 1
rxn A -> B k
every 0.1 save
'''

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    for i, (initial, dose) in enumerate(CASES):
        p = td / f"u{i}.model"
        p.write_text(BASE.replace("param t_end", f"param init_volume {initial}\nparam t_end") + f"at 0.2 feed F {dose}\n")
        r = subprocess.run([str(ENGINE), "--check", str(p)], text=True, capture_output=True)
        assert r.returncode == 0, (p, r.stdout, r.stderr)
    bad = td / "bad.model"
    bad.write_text(BASE.replace("param t_end", "param init_volume 100 uL\nparam t_end"))
    r = subprocess.run([str(ENGINE), "--check", str(bad)], text=True, capture_output=True)
    assert r.returncode != 0
    assert "unit must be L or mL" in (r.stdout + r.stderr)
print("[OK] process volume units")
