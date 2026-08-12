## Pure dispatch/routing logic for the slimmc CLI: which engine (homo or
## copo) should handle a given .model file, based only on counting
## declared `monomer NAME ...` lines. No engine code is imported here --
## this stays a small, standalone module so the logic-only unit tests in
## tests/test_dispatch.nim don't need to link homo or copo at all.
##
## Dispatch rule: count non-comment `monomer NAME ...` declarations.
##   0-1 monomers -> homo (0 is a legitimate case: pure-kinetics models
##                  with no polymer chain growth at all -- confirmed by
##                  homo's own validation suite, where 20 of 44 engine
##                  validation cases declare zero monomers; an earlier
##                  version of this rule treated 0 as a dispatcher-level
##                  error, which silently worked only because nothing
##                  routed those models through the unified dispatcher)
##   2-3 monomers -> copo (copo hard-requires >=2, confirmed in copo_parser.nim,
##                  so 0 or 1 can never legitimately mean "this was a copo model")
##   >3 monomers  -> clear dispatcher-level error, before either engine runs
##
## This is a sniff, not a parser: it does not validate the model beyond
## counting monomer lines. Full validation happens in whichever engine
## gets dispatched to -- this module's only job is picking the right one.

import std/strutils

proc stripComment*(line: string): string =
  ## `monomer` lines never contain quoted strings, so a plain "cut at the
  ## first #" is sufficient here (unlike full model-line parsing
  ## elsewhere in the codebase, which has to respect quotes for `desc`).
  let idx = line.find('#')
  result = (if idx >= 0: line[0 ..< idx] else: line).strip()

proc countMonomers*(path: string): int =
  ## Count non-comment `monomer NAME ...` declaration lines. Returns -1
  ## if the file cannot be read at all (distinct from 0 monomers found).
  var handle: File
  if not open(handle, path, fmRead):
    return -1
  defer: handle.close()
  var count = 0
  for rawLine in handle.lines:
    let line = stripComment(rawLine)
    if line.len == 0:
      continue
    let firstToken = line.split(WhiteSpace, maxsplit = 1)[0]
    if firstToken == "monomer":
      inc count
  result = count

proc engineForCount*(n: int): string =
  case n
  of 0, 1: "homo"
  of 2, 3: "copo"
  else: ""

