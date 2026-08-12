#!/usr/bin/env python3
from pathlib import Path
import subprocess, tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "bin" / "slimmc"
BASE = (ROOT / "homo/tests/validation/phase_d/models/H18_rxn_uni.model").read_text()

with tempfile.TemporaryDirectory(prefix="slimmc-stage24-") as td:
    td = Path(td)
    model = td / "stage24.model"
    model.write_text(BASE + "\nparam init_volume 0.100\nfeed F A 1.0\nat 0.20 feed F 0.001\nat 0.30 set_c A 0.005\n")
    check = subprocess.run([str(ENGINE), "--check", str(model)], text=True, capture_output=True)
    assert check.returncode == 0, check.stderr
    for section in ["GENERAL", "DETAILS", "WARNINGS", "ERRORS", "LIMITATIONS", "CHECK RESULT"]:
        assert section in check.stdout
    assert "set_c forces" in check.stdout
    run = subprocess.run([str(ENGINE), str(model)], text=True, capture_output=True)
    assert run.returncode == 0, run.stderr
    assert "WARNING:" in run.stdout and "set_c" in run.stdout
    out = td / "results" / "stage24"
    dose = np.load(out / "feed_events" / "dose_mL.npy")
    assert dose.tolist() == [1.0]
print("stage24 check/units: PASS")
