from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN = (
    "legacy" + "-tsv", "legacy" + "_tsv", "runinfo" + ".json", "start" + ".json",
    "state" + ".table", "chains" + ".table", "moments" + ".table", "memory" + ".table",
    "firings" + ".table", "parameter_states" + ".table", "channel_trace" + ".table",
    "validation" + ".table", "validation_details" + ".table",
    "--trace" + "-actions", "--" + "log",
)
for path in ROOT.rglob("*"):
    if not path.is_file() or ".git" in path.parts:
        continue
    if path.suffix in {".zip", ".png", ".pdf", ".npy", ".pyc"}:
        continue
    text = path.read_text(errors="ignore")
    for token in FORBIDDEN:
        assert token not in text, f"{token!r} in {path}"
print("storage-only repository check: PASS")
