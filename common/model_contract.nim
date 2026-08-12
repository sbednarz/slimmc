## Shared lexical and default-value contract for Slimmc model parsers.
## Chemical validation remains engine-specific.

import strutils

const
  ModelIdentifierPattern* = "[A-Za-z_][A-Za-z0-9_]*"
  DefaultTemperatureK* = 298.15
  DefaultMaxSteps* = 10_000_000_000'i64
  DefaultWhenCheckEvents* = 1'i64
  DefaultSeed* = 12345'i64
  ## Public DP storage and parser contract is int64.  The default remains
  ## deliberately bounded to the portable int32 range until both engines'
  ## chain internals are fully int64-clean.
  DefaultDpMax* = int64(high(int32))
  DefaultSequenceMode* = "composition"
  DefaultMassModel* = "repeat_units"

const ReservedModelIdentifiers* = [
  "param", "desc", "var", "monomer", "species", "endgroup", "polymer",
  "rate", "rxn", "macro", "feed", "every", "from", "repeat", "step", "at", "when",
  "output_dir", "memory_limit", "at_memory", "init", "prop", "deprop",
  "term_c", "term_d", "term_x", "transfer", "transfer_h", "transfer_m",
  "set_k", "add_k", "set_temp", "add_temp", "set_c", "add_c", "print",
  "print_info", "save", "save_chains", "stop", "print_memory", "memory",
  "active", "dead", "const", "arr", "repeat_units", "with_end_groups",
  "full", "composition"
]

proc isModelIdentStart*(c: char): bool {.inline.} =
  (c >= 'A' and c <= 'Z') or (c >= 'a' and c <= 'z') or c == '_'

proc isModelIdentChar*(c: char): bool {.inline.} =
  isModelIdentStart(c) or (c >= '0' and c <= '9')

proc isValidModelIdentifier*(name: string): bool =
  if name.len == 0 or not isModelIdentStart(name[0]):
    return false
  for c in name:
    if not isModelIdentChar(c):
      return false
  true


proc isValidModelPath*(path: string): bool =
  ## Portable public path contract for paths written inside a .model file.
  ## The path itself may be absolute or relative and may use '/' or '\\'
  ## separators. Every non-root segment must follow the model identifier
  ## grammar. Dot segments, extensions, spaces, punctuation and glob
  ## metacharacters are deliberately rejected.
  if path.len == 0:
    return false

  var normalized = path.replace('\\', '/')
  var start = 0

  # Optional Windows drive prefix, e.g. C:/results/run_01.
  if normalized.len >= 2 and normalized[0].isAlphaAscii and normalized[1] == ':':
    start = 2
    if normalized.len == 2:
      return false
    if normalized[2] == '/':
      start = 3

  # POSIX root or UNC prefix. Empty leading segments are structural only.
  while start < normalized.len and normalized[start] == '/':
    inc start
  if start >= normalized.len:
    return false

  let rest = normalized[start .. ^1]
  for segment in rest.split('/'):
    if segment.len == 0 or not isValidModelIdentifier(segment):
      return false
  true

proc isReservedModelIdentifier*(name: string): bool =
  ## Language keywords remain case-sensitive.
  for word in ReservedModelIdentifiers:
    if name == word:
      return true
  false

proc modelNameKey*(name: string): string {.inline.} =
  ## Used only for collision detection. The language itself remains
  ## case-sensitive, but names differing solely by ASCII case are rejected.
  name.toLowerAscii()
