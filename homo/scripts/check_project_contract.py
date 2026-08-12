from __future__ import annotations
from pathlib import Path
import re

# v2.8.5: removed syntax is not allowed in executable model lines.
LEGACY_MODEL_PATTERNS = [
    (re.compile(r'^\s*species\s+[A-Za-z][A-Za-z0-9_]*\s+[^\s#]+\s+[^\s#]+(?:\s|$)'), 'species NAME c0 MW'),
    (re.compile(r'^\s*outputPrefix\b'), 'outputPrefix'),
    (re.compile(r'^\s*param\s+tEnd\b'), 'param tEnd'),
    (re.compile(r'^\s*param\s+maxSteps\b'), 'param maxSteps'),
    (re.compile(r'^\s*param\s+max_events\b'), 'param max_events'),
    (re.compile(r'^\s*param\s+massModel\b'), 'param massModel'),
    (re.compile(r'^\s*param\s+outputDir\b'), 'param outputDir'),
    (re.compile(r'^\s*macro\s+termc\b'), 'macro termc'),
    (re.compile(r'^\s*macro\s+termd\b'), 'macro termd'),
    (re.compile(r'^\s*macro\s+termx\b'), 'macro termx'),
    (re.compile(r'^\s*macro\s+termC\b'), 'macro termC'),
    (re.compile(r'^\s*macro\s+termD\b'), 'macro termD'),
    (re.compile(r'^\s*macro\s+termX\b'), 'macro termX'),
    (re.compile(r'^\s*var\s+[A-Za-z][A-Za-z0-9_]*\s+[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\b'), 'removed var NAME VALUE [UNIT]'),
]

for path in Path('.').rglob('*.model'):
    text = path.read_text(encoding='utf-8', errors='ignore')
    for lineno, line in enumerate(text.splitlines(), 1):
        stripped = line.lstrip()
        if not stripped or stripped.startswith('#'):
            continue
        for pat, label in LEGACY_MODEL_PATTERNS:
            if pat.search(line):
                raise AssertionError(f'{path}:{lineno}: removed syntax is not allowed in v2.8.5: {label}')

# Canonical family-level docs should not advertise old names as supported.
ROOT = Path(__file__).resolve().parents[1]
FAMILY = ROOT.parent
CANONICAL_DOCS = [
    FAMILY / 'README.md',
    FAMILY / 'docs' / 'reference' / 'HOMO.md',
    FAMILY / 'docs' / 'OUTPUT_FORMAT.md',
    FAMILY / 'docs' / 'PYSLIMMC.md',
]
BAD_DOC_TOKENS = ['outputPrefix', 'tEnd', 'maxSteps', 'max_events', 'massModel', 'outputDir', 'termC', 'termD', 'termX', 'axis ']
for path in CANONICAL_DOCS:
    if not path.exists():
        continue
    text = path.read_text(encoding='utf-8', errors='ignore')
    for tok in BAD_DOC_TOKENS:
        if tok in text:
            raise AssertionError(f'{path}: removed token in canonical documentation: {tok}')

print('[project-contract] ok')
