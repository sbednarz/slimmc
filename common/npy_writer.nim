## Minimal one-dimensional NumPy NPY v1.0 writer for Slimmc Storage v1.
## Supported canonical dtypes: uint64, uint32, float64 and bool.

import std/[os, strutils]
import ./atomic_file

const
  NpyMagic = "\x93NUMPY"
  NpyMajor = 1'u8
  NpyMinor = 0'u8

proc ensureParentDir(path: string) =
  let parent = parentDir(path)
  if parent.len > 0 and not dirExists(parent):
    createDir(parent)

proc addLe16(buf: var string; value: uint16) =
  buf.add(char(value and 0xff'u16))
  buf.add(char((value shr 8) and 0xff'u16))

proc addLe32(buf: var string; value: uint32) =
  for shift in countup(0, 24, 8):
    buf.add(char((value shr shift) and 0xff'u32))

proc addLe64(buf: var string; value: uint64) =
  for shift in countup(0, 56, 8):
    buf.add(char((value shr shift) and 0xff'u64))

proc floatBits(value: float64): uint64 =
  copyMem(addr result, unsafeAddr value, sizeof(value))

proc npyHeader(descr: string; length: int): string =
  if length < 0:
    raise newException(ValueError, "NPY array length must be nonnegative")
  var dictionary = "{'descr': '" & descr & "', 'fortran_order': False, 'shape': (" & $length & ",), }"
  # NPY 1.0 requires magic + version + uint16 header length + header
  # to end on a 16-byte boundary. Header itself ends with newline.
  let prefixLength = NpyMagic.len + 2 + 2
  let pad = (16 - ((prefixLength + dictionary.len + 1) mod 16)) mod 16
  dictionary.add(repeat(' ', pad))
  dictionary.add('\n')
  if dictionary.len > int(high(uint16)):
    raise newException(ValueError, "NPY 1.0 header is too large")
  result = NpyMagic
  result.add(char(NpyMajor))
  result.add(char(NpyMinor))
  addLe16(result, uint16(dictionary.len))
  result.add(dictionary)

proc publishAtomic(path, payload: string) =
  ensureParentDir(path)
  let tmpPath = path & ".tmp"
  if fileExists(tmpPath):
    removeFile(tmpPath)
  try:
    var f = open(tmpPath, fmWrite)
    try:
      f.write(payload)
      f.flushFile()
    finally:
      f.close()
    atomicReplaceFile(tmpPath, path)
  except:
    if fileExists(tmpPath):
      removeFile(tmpPath)
    raise

proc writeNpyUint64*(path: string; values: openArray[uint64]) =
  var payload = npyHeader("<u8", values.len)
  for value in values:
    addLe64(payload, value)
  publishAtomic(path, payload)

proc writeNpyUint32*(path: string; values: openArray[uint32]) =
  var payload = npyHeader("<u4", values.len)
  for value in values:
    addLe32(payload, value)
  publishAtomic(path, payload)

proc writeNpyFloat64*(path: string; values: openArray[float64]) =
  var payload = npyHeader("<f8", values.len)
  for value in values:
    addLe64(payload, floatBits(value))
  publishAtomic(path, payload)

proc writeNpyBool*(path: string; values: openArray[bool]) =
  var payload = npyHeader("|b1", values.len)
  for value in values:
    payload.add(if value: '\x01' else: '\x00')
  publishAtomic(path, payload)
