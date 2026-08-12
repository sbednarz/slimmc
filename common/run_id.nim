## Canonical run_id validation shared by the dispatcher and both engines.
##
## A model file stem is the run_id and must be a valid ASCII Python
## identifier. This keeps directory names, interactive ``runs.<run_id>``
## access, and glob-based selection unambiguous and portable.

import std/os

const RunIdPattern* = "[A-Za-z_][A-Za-z0-9_]*"

proc isAsciiLetter(c: char): bool {.inline.} =
  (c >= 'A' and c <= 'Z') or (c >= 'a' and c <= 'z')

proc isRunIdStart(c: char): bool {.inline.} =
  isAsciiLetter(c) or c == '_'

proc isRunIdChar(c: char): bool {.inline.} =
  isRunIdStart(c) or (c >= '0' and c <= '9')

proc isValidRunId*(runId: string): bool =
  if runId.len == 0 or not isRunIdStart(runId[0]):
    return false
  for c in runId:
    if not isRunIdChar(c):
      return false
  true

proc runIdFromModelPath*(modelPath: string): string =
  splitFile(modelPath).name

proc validateRunId*(runId: string) =
  if not isValidRunId(runId):
    raise newException(ValueError,
      "invalid run_id '" & runId & "'; run_id must match " & RunIdPattern)

proc validateModelRunId*(modelPath: string): string =
  result = runIdFromModelPath(modelPath)
  if not isValidRunId(result):
    raise newException(ValueError,
      "invalid model filename '" & extractFilename(modelPath) & "': " &
      "the model filename determines run_id, which must match " & RunIdPattern)
