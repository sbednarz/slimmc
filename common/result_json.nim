## Shared JSON/JSONL writers for Slimmc Slimmc Storage v1.
##
## `run_metadata.json` is the sole human-oriented, pretty-printed JSON file.
## Other structured result streams use JSONL (one object per non-empty line).

import std/[json, os]
import ./atomic_file

proc ensureParentDir(path: string) =
  let parent = parentDir(path)
  if parent.len > 0 and not dirExists(parent):
    createDir(parent)

proc writeTextAtomic*(path, text: string) =
  ## Write complete UTF-8 text through `<path>.tmp`, then rename.
  ensureParentDir(path)
  let tmpPath = path & ".tmp"
  if fileExists(tmpPath):
    removeFile(tmpPath)
  try:
    var f = open(tmpPath, fmWrite)
    try:
      f.write(text)
      f.flushFile()
    finally:
      f.close()
    atomicReplaceFile(tmpPath, path)
  except:
    if fileExists(tmpPath):
      removeFile(tmpPath)
    raise

proc writePrettyJsonAtomic*(path: string; value: JsonNode) =
  ## Write one multiline JSON object with two-space indentation.
  ## A trailing newline keeps the file pleasant in text tools.
  writeTextAtomic(path, pretty(value, indent = 2) & "\n")

proc jsonlText*(records: openArray[JsonNode]): string =
  ## Serialize records as compact JSONL. Blank lines are never emitted.
  for record in records:
    if record.isNil:
      raise newException(ValueError, "JSONL record must not be nil")
    result.add($record)
    result.add('\n')

proc writeJsonlAtomic*(path: string; records: openArray[JsonNode]) =
  writeTextAtomic(path, jsonlText(records))
