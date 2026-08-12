## Canonical checksum manifest for one Slimmc Storage run.
##
## The manifest deliberately excludes files that would create a hash cycle or
## are completion/working markers: run_metadata.json, checksums.sha256,
## RESULTS_COMPLETE, and everything below .work/.

import std/[os, algorithm, strutils]
import sha256_file
import result_json

const
  StorageHashAlgorithm* = "sha256"
  StorageHashSchema* = "slimmc-storage-hash-v1"
  StorageChecksumFile* = "checksums.sha256"

type StorageManifest* = object
  hash*: string
  fileCount*: int
  totalBytes*: int64
  text*: string

proc normalizedRelative(root, path: string): string =
  result = relativePath(path, root).replace('\\', '/')

proc isManifestPayload(relative: string): bool =
  if relative in ["run_metadata.json", StorageChecksumFile, "RESULTS_COMPLETE"]:
    return false
  if relative == ".work" or relative.startsWith(".work/"):
    return false
  true

proc buildStorageManifest*(root: string): StorageManifest =
  var paths: seq[string] = @[]
  for path in walkDirRec(root, yieldFilter = {pcFile}, followFilter = {pcDir}):
    let relative = normalizedRelative(root, path)
    if isManifestPayload(relative):
      paths.add relative
  paths.sort(system.cmp[string])

  var lines = newStringOfCap(paths.len * 96)
  var total = 0'i64
  for relative in paths:
    let absolute = root / relative.replace('/', DirSep)
    let digest = sha256File(absolute)
    lines.add digest & "  " & relative & "\n"
    total += getFileSize(absolute)

  let payload = StorageHashSchema & "\n" & lines
  result.hash = sha256Text(payload)
  result.fileCount = paths.len
  result.totalBytes = total
  result.text = lines

proc writeStorageManifest*(root: string): StorageManifest =
  result = buildStorageManifest(root)
  writeTextAtomic(root / StorageChecksumFile, result.text)
