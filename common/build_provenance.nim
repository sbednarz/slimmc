import std/[os, json]
import sha256_file

const
  BuildGitCommit* {.strdefine.} = ""
  BuildGitTag* {.strdefine.} = ""
  BuildGitDirty* {.strdefine.} = "unknown"
  BuildTimestampUtc* {.strdefine.} = ""

proc executableProvenance*(): JsonNode =
  result = newJObject()
  let exe = getAppFilename()
  result["binary_name"] = %extractFilename(exe)
  result["binary_hash_algorithm"] = %"sha256"
  if fileExists(exe):
    result["binary_hash"] = %sha256File(exe)
  else:
    result["binary_hash"] = newJNull()
  if BuildGitCommit.len > 0: result["git_commit"] = %BuildGitCommit
  if BuildGitTag.len > 0: result["git_tag"] = %BuildGitTag
  if BuildGitDirty == "true": result["git_dirty"] = %true
  elif BuildGitDirty == "false": result["git_dirty"] = %false
  else: result["git_dirty"] = newJNull()
  if BuildTimestampUtc.len > 0: result["build_timestamp_utc"] = %BuildTimestampUtc
