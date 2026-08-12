# copo_sequence.nim
# Linear sequence operations for slimmc-copo v0.2.

import strutils
import copo_types

proc initSequence*(): LinearSequence =
  result.data = @[]

proc initSequence*(m: int): LinearSequence =
  result.data = @[uint8(m)]

proc len*(s: LinearSequence): int {.inline.} = s.data.len

proc push*(s: var LinearSequence; m: int) {.inline.} =
  s.data.add uint8(m)

proc pop*(s: var LinearSequence): int {.inline.} =
  result = int(s.data[^1])
  s.data.setLen(s.data.len - 1)

proc last*(s: LinearSequence): int {.inline.} =
  if s.data.len == 0: return -1
  int(s.data[^1])

proc prev*(s: LinearSequence): int {.inline.} =
  if s.data.len < 2: return -1
  int(s.data[^2])

proc first*(s: LinearSequence): int {.inline.} =
  if s.data.len == 0: return -1
  int(s.data[0])

proc reverseCopy*(s: LinearSequence): LinearSequence =
  result.data = newSeq[uint8](s.data.len)
  for i in 0 ..< s.data.len:
    result.data[i] = s.data[s.data.len - 1 - i]

proc append*(a: var LinearSequence; b: LinearSequence) =
  for x in b.data:
    a.data.add x

proc appendReverse*(a: var LinearSequence; b: LinearSequence) =
  if b.data.len == 0: return
  for i in countdown(b.data.len - 1, 0):
    a.data.add b.data[i]

proc toText*(s: LinearSequence; monomerNames: seq[string]): string =
  ## Renders as monomer names joined by "|" (e.g. "MAA|AA|MAA"), not bare
  ## concatenation ("MAAAAMAA"). Direct concatenation was provably
  ## ambiguous for any monomer-name set with prefix overlaps (e.g. names
  ## "A"/"AB"/"B": the string "AB" could mean either token "AB" alone or
  ## tokens "A" then "B" -- no downstream reader, however clever, can
  ## recover the real answer from bare-concatenated text alone). The "|"
  ## character is disallowed in monomer names (see parser validation) so
  ## it can never collide with a real name.
  var parts: seq[string] = @[]
  for x in s.data:
    let i = int(x)
    if i >= 0 and i < monomerNames.len:
      parts.add monomerNames[i]
    else:
      parts.add "?"
  result = parts.join("|")

proc countMers*(s: LinearSequence; nMonomers: int): seq[int32] =
  result = newSeq[int32](nMonomers)
  for x in s.data:
    let i = int(x)
    if i >= 0 and i < nMonomers:
      inc result[i]

proc dyadCounts*(s: LinearSequence; nMonomers: int): seq[int64] =
  result = newSeq[int64](nMonomers * nMonomers)
  if s.data.len < 2: return
  for i in 0 ..< s.data.len - 1:
    let a = int(s.data[i])
    let b = int(s.data[i + 1])
    inc result[a * nMonomers + b]

proc triadCounts*(s: LinearSequence; nMonomers: int): seq[int64] =
  result = newSeq[int64](nMonomers * nMonomers * nMonomers)
  if s.data.len < 3: return
  for i in 0 ..< s.data.len - 2:
    let a = int(s.data[i])
    let b = int(s.data[i + 1])
    let c = int(s.data[i + 2])
    inc result[(a * nMonomers + b) * nMonomers + c]


proc addBlockCount(acc: var seq[BlockCount]; monomerId: int; length: int32) =
  if length <= 0: return
  for i in 0 ..< acc.len:
    if acc[i].monomerId == monomerId and acc[i].length == length:
      acc[i].count += 1
      return
  acc.add BlockCount(monomerId: monomerId, length: length, count: 1)

proc blockCounts*(s: LinearSequence; nMonomers: int): seq[BlockCount] =
  ## Exact contiguous-run histogram by monomer.  Unlike sequenceText, this is
  ## intended to be preserved in dead summaries even when full sequence text is
  ## omitted when sequence_mode=composition.
  if s.data.len == 0: return
  var cur = int(s.data[0])
  var runLen = 1'i32
  for i in 1 ..< s.data.len:
    let m = int(s.data[i])
    if m == cur:
      inc runLen
    else:
      if cur >= 0 and cur < nMonomers:
        addBlockCount(result, cur, runLen)
      cur = m
      runLen = 1
  if cur >= 0 and cur < nMonomers:
    addBlockCount(result, cur, runLen)
