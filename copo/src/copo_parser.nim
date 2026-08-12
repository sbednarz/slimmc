# copo_parser.nim
# Model parser for slimmc-copo v0.6.13.

import os, strutils, tables, math
import copo_types
import ../../common/safe_numeric
import ../../common/run_id
import ../../common/model_contract

proc parseFail(lineNo: int; msg: string) =
  if lineNo > 0:
    raise newException(ValueError, "line " & $lineNo & ": " & msg)
  else:
    raise newException(ValueError, msg)

proc parseF(s0: string; lineNo: int = 0; what: string = "number"): float =
  ## Parse a locale-independent floating-point number with classic slimmc-style diagnostics.
  let s = s0.strip()
  if s.len == 0:
    parseFail(lineNo, "empty " & what)

  var i = 0
  var sign = 1.0

  if s[i] == '+':
    inc i
  elif s[i] == '-':
    sign = -1.0
    inc i

  var intPart = 0.0
  var hasDigit = false

  while i < s.len and s[i].isDigit:
    hasDigit = true
    intPart = intPart * 10.0 + float(ord(s[i]) - ord('0'))
    inc i

  var fracPart = 0.0
  var fracScale = 1.0

  if i < s.len and s[i] == '.':
    inc i
    while i < s.len and s[i].isDigit:
      hasDigit = true
      fracPart = fracPart * 10.0 + float(ord(s[i]) - ord('0'))
      fracScale *= 10.0
      inc i

  if not hasDigit:
    parseFail(lineNo, "invalid " & what & ": " & s0)

  var expVal = 0
  if i < s.len and (s[i] == 'e' or s[i] == 'E'):
    inc i
    var expSign = 1
    if i < s.len and s[i] == '+':
      inc i
    elif i < s.len and s[i] == '-':
      expSign = -1
      inc i

    var hasExpDigit = false
    while i < s.len and s[i].isDigit:
      hasExpDigit = true
      expVal = expVal * 10 + (ord(s[i]) - ord('0'))
      inc i

    if not hasExpDigit:
      parseFail(lineNo, "invalid exponent in " & what & ": " & s0)

    expVal *= expSign

  if i != s.len:
    parseFail(lineNo, "invalid trailing characters in " & what & ": " & s0)

  result = sign * (intPart + fracPart / fracScale) * pow(10.0, float(expVal))

proc parseI64(s0: string; lineNo: int = 0; what: string = "integer"): int64 =
  try:
    result = checkedParseInt64(s0, what)
  except ValueError as exc:
    parseFail(lineNo, exc.msg)

proc parseSizeBytes*(s0: string): int64 =
  var s = s0.strip().toLowerAscii()
  var mult = 1.0
  if s.endsWith("kb"):
    mult = 1024.0; s = s[0 .. ^3]
  elif s.endsWith("mb"):
    mult = 1024.0 * 1024.0; s = s[0 .. ^3]
  elif s.endsWith("gb"):
    mult = 1024.0 * 1024.0 * 1024.0; s = s[0 .. ^3]
  elif s.endsWith("b"):
    mult = 1.0; s = s[0 .. ^2]
  result = int64(parseF(s.strip(), 0, "size") * mult)

proc stripComment(line: string): string =
  var inQuote = false
  var quoteChar = '\0'
  var escaped = false
  for i, ch in line:
    if escaped:
      escaped = false
    elif inQuote and ch == '\\':
      escaped = true
    elif ch == '"' or ch == '\'':
      if inQuote and ch == quoteChar:
        inQuote = false
      elif not inQuote:
        inQuote = true
        quoteChar = ch
    elif ch == '#' and not inQuote:
      return line[0 ..< i].strip()
  result = line.strip()

proc require(cond: bool; lineNo: int; msg: string) =
  if not cond:
    raise newException(ValueError, "line " & $lineNo & ": " & msg)

proc tokenize*(line: string; lineNo: int): seq[string] =
  ## Splits a model line into tokens. Comments start at "#" only outside
  ## quoted strings (item 35). "->" and "+" are always split into their
  ## own tokens regardless of surrounding whitespace (item 33/34: "A+B->C",
  ## "A + B -> C", and whitespace-separated variants must all tokenize
  ## identically) -- this parser's statement handlers expect "+"/"->" as
  ## standalone tokens at fixed positions, so the split has to happen
  ## here, not via a later join-and-resplit trick.
  var tokens: seq[string] = @[]
  var cur = ""
  var inQuote = false
  var quoteChar = '\0'
  var escaped = false
  var i = 0

  proc flush() =
    if cur.len > 0:
      tokens.add cur
      cur = ""

  while i < line.len:
    let ch = line[i]
    if escaped:
      case ch
      of 'n': cur.add('\n')
      of 'r': cur.add('\r')
      of 't': cur.add('\t')
      of '"': cur.add('"')
      of '\'': cur.add('\'')
      of '\\': cur.add('\\')
      else: cur.add(ch)
      escaped = false
      inc i
      continue
    if inQuote:
      if ch == '\\':
        escaped = true
      elif ch == quoteChar:
        inQuote = false
        flush()
      else:
        cur.add ch
      inc i
      continue
    if ch == '#':
      break
    elif ch == '"' or ch == '\'':
      inQuote = true
      quoteChar = ch
    elif ch == ' ' or ch == '\t' or ch == '\r' or ch == '\n':
      flush()
    elif ch == '-' and i + 1 < line.len and line[i + 1] == '>':
      flush()
      tokens.add "->"
      inc i
    elif ch == '+' and cur.len >= 2 and (cur[^1] == 'e' or cur[^1] == 'E') and
         (cur[^2].isDigit or cur[^2] == '.'):
      # Part of a float's exponent (e.g. "1.0e+5"), not a term separator
      # -- distinguished from "Fe+Cl"-style identifier-then-plus by
      # checking what's actually accumulated in `cur` so far: a numeric
      # exponent has a digit/'.' right before the 'e'/'E', an identifier
      # like "Fe" does not.
      cur.add ch
    elif ch == '+':
      flush()
      tokens.add "+"
    else:
      cur.add ch
    inc i
  require(not inQuote, lineNo, "unterminated quoted string")
  require(not escaped, lineNo, "unterminated escape in quoted string")
  flush()
  result = tokens

proc validateName(name: string; lineNo: int; what: string) =
  require(isValidModelIdentifier(name), lineNo,
    what & " name must match " & ModelIdentifierPattern & ": " & name)
  require(not isReservedModelIdentifier(name), lineNo,
    what & " name is reserved: " & name)

proc findDeclaredNameCaseInsensitive(m: Model; name: string): string =
  let key = modelNameKey(name)
  for existing in m.monomerByName.keys:
    if modelNameKey(existing) == key: return existing
  for existing in m.speciesByName.keys:
    if modelNameKey(existing) == key: return existing
  for existing in m.poolByName.keys:
    if modelNameKey(existing) == key: return existing
  for existing in m.rateByName.keys:
    if modelNameKey(existing) == key: return existing
  ""

proc requireCoreNameFree(m: Model; name: string; lineNo: int; what: string) =
  let existing = m.findDeclaredNameCaseInsensitive(name)
  require(existing.len == 0, lineNo,
    what & " name conflicts with already declared model name: " & existing)

proc parseMassModelValue*(value: string; lineNo: int): MassModel =
  ## Only the exact, case-sensitive model-language spellings are accepted.
  case value
  of "repeat_units":
    return mmRepeatUnits
  of "with_end_groups":
    return mmWithEndgroups
  else:
    raise newException(ValueError, "line " & $lineNo & ": mass_model must be repeat_units or with_end_groups")

proc defaultOutputDir(filename: string): string =
  let sf = splitFile(filename)
  let base = if sf.name.len > 0: sf.name else: "run"
  if sf.dir.len > 0:
    result = sf.dir / "results" / base
  else:
    result = "results" / base

proc addMonomer(m: var Model; name: string; c0, mw: float; lineNo: int) =
  validateName(name, lineNo, "monomer")
  m.requireCoreNameFree(name, lineNo, "monomer")
  m.monomerByName[name] = m.monomers.len
  m.monomers.add MonomerDef(name: name, c0: c0, mw: mw)

proc addSpecies(m: var Model; name: string; c0: float; lineNo: int) =
  validateName(name, lineNo, "species")
  m.requireCoreNameFree(name, lineNo, "species")
  m.speciesByName[name] = m.species.len
  m.species.add SpeciesDef(name: name, c0: c0)

proc addEndGroup(m: var Model; name: string; mw: float; lineNo: int) =
  validateName(name, lineNo, "endgroup")
  let egKey = modelNameKey(name)
  for existing in m.endgroupByName.keys:
    require(modelNameKey(existing) != egKey, lineNo,
      "endgroup name conflicts with already declared endgroup: " & existing)
  require(not m.endgroupByName.hasKey(name), lineNo, "duplicate endgroup: " & name)
  m.endgroupByName[name] = m.endgroups.len
  m.endgroups.add EndGroupDef(name: name, mw: mw)

proc addPool(m: var Model; name: string; kind: PoolKind; lineNo: int) =
  validateName(name, lineNo, "polymer pool")
  m.requireCoreNameFree(name, lineNo, "polymer pool")
  m.poolByName[name] = m.pools.len
  if kind == pkDead:
    m.deadPoolId = m.pools.len
  m.pools.add PoolDef(name: name, kind: kind)
  m.poolTerminalMer.add -1
  m.poolPenultimateMer.add -1

proc addRate(m: var Model; r: RateDef; lineNo: int) =
  validateName(r.name, lineNo, "rate")
  m.requireCoreNameFree(r.name, lineNo, "rate")
  m.rateByName[r.name] = m.rates.len
  m.rates.add r

proc addFixedRate(m: var Model; name: string; k: float; lineNo: int) =
  require(k >= 0.0, lineNo, "rate value must be >= 0")
  m.addRate(RateDef(name: name, kind: rkFixed, kConst: k), lineNo)

proc addArrheniusRate(m: var Model; name: string; apre, ea: float; lineNo: int) =
  require(apre >= 0.0, lineNo, "Arrhenius pre-exponential factor must be >= 0")
  m.addRate(RateDef(name: name, kind: rkArr, Apre: apre, Ea: ea, declaredArrhenius: true), lineNo)

proc monomerId(m: Model; name: string; lineNo: int): int =
  require(m.monomerByName.hasKey(name), lineNo, "unknown monomer: " & name)
  m.monomerByName[name]

proc speciesId(m: Model; name: string; lineNo: int): int =
  require(m.speciesByName.hasKey(name), lineNo, "unknown species: " & name)
  m.speciesByName[name]

proc poolId(m: Model; name: string; lineNo: int): int =
  require(m.poolByName.hasKey(name), lineNo, "unknown polymer pool: " & name)
  m.poolByName[name]

proc rateId(m: Model; name: string; lineNo: int): int =
  require(m.rateByName.hasKey(name), lineNo, "unknown rate: " & name)
  m.rateByName[name]

proc setPoolTerminal(m: var Model; poolId, merId, lineNo: int) =
  if poolId < 0 or poolId >= m.poolTerminalMer.len: return
  if m.pools[poolId].kind != pkActive: return
  let old = m.poolTerminalMer[poolId]
  if old == -1:
    m.poolTerminalMer[poolId] = merId
  else:
    require(old == merId, lineNo, "active pool " & m.pools[poolId].name &
      " is assigned inconsistent terminal monomers")

proc setPoolPenultimate(m: var Model; poolId, merId, lineNo: int) =
  if poolId < 0 or poolId >= m.poolPenultimateMer.len: return
  if m.pools[poolId].kind != pkActive: return
  let old = m.poolPenultimateMer[poolId]
  if old == -1:
    m.poolPenultimateMer[poolId] = merId
  elif old != merId:
    # -2 means mixed penultimate pool. This is expected for terminal-model
    # pools such as PA receiving both AA* and BA* chains. Explicit
    # penultimate pools remain single-valued when the model uses them.
    m.poolPenultimateMer[poolId] = -2

proc parseActionKind*(s: string; lineNo: int): ActionKind =
  ## Canonical action names only (item 55) -- "snapshot"/"chains"/"memory"
  ## used to be accepted as aliases of "save"/"save_chains"/"print_memory";
  ## removed per the agreed contract cleanup (matches homo's identical
  ## removal). Also made case-sensitive (item 39): this used to match on
  ## s.toLowerAscii(), so "SAVE"/"Save" were silently accepted as actions
  ## too -- an unjustified parser asymmetry with homo's plain,
  ## case-sensitive `case s` match, not a chemical difference.
  case s
  of "print": result = eaPrint
  of "print_info": result = eaPrintInfo
  of "save": result = eaSave
  of "save_chains": result = eaSaveChains
  of "stop": result = eaStop
  of "print_memory": result = eaPrintMemory
  of "set_k": result = eaSetK
  of "add_k": result = eaAddK
  of "set_temp": result = eaSetTemp
  of "add_temp": result = eaAddTemp
  of "set_c": result = eaSetC
  of "add_c": result = eaAddC
  of "feed": result = eaFeed
  else:
    raise newException(ValueError, "line " & $lineNo & ": unknown action: " & s)

proc validateActionArgs(action: ActionKind; args: seq[string]; lineNo: int) =
  ## Items 56/57/58: "print" bare (no message) is now a parser error --
  ## the progress-printout behavior it used to double as moved to the
  ## new, distinct "print_info" action (matches homo's identical fix).
  ## "print" now always requires exactly one message argument.
  case action
  of eaPrint:
    if args.len != 1:
      raise newException(ValueError, "line " & $lineNo & ": print syntax: print \"message\" (bare print with no message is no longer valid -- use print_info for a progress printout)")
  of eaPrintInfo:
    if args.len != 0:
      raise newException(ValueError, "line " & $lineNo & ": print_info takes no arguments")
  of eaFeed:
    if args.len notin [2, 3]:
      raise newException(ValueError, "line " & $lineNo & ": feed action syntax: feed NAME VOLUME [L|l|mL|ml|ML]; without a unit VOLUME is in L")
  of eaStop:
    if args.len != 0:
      raise newException(ValueError, "line " & $lineNo & ": stop takes no arguments")
  else:
    discard

proc parseScheduledAction(m: var Model; toks: seq[string]; lineNo: int) =
  require(toks.len >= 3, lineNo, "scheduled action syntax: every PERIOD ACTION ... OR from START step PERIOD ACTION ... OR at TIME ACTION ...")
  case toks[0]
  of "every":
    let period = parseF(toks[1], lineNo, "scheduled action numeric argument")
    require(period > 0.0, lineNo, "every PERIOD must be > 0")
    let action = parseActionKind(toks[2], lineNo)
    if action == eaStop:
      raise newException(ValueError, "line " & $lineNo & ": stop is only valid as a when action")
    let args = if toks.len > 3: toks[3 .. ^1] else: @[]
    validateActionArgs(action, args, lineNo)
    # Periodic schedules start at t=0. `every DT ACTION` is the compact
    # spelling of `from 0 step DT ACTION`.
    m.scheduledActions.add ScheduledAction(
      startTime: 0.0, nextTime: 0.0, period: period, repeat: true, remaining: -1, active: true,
      action: action, args: args, lineNo: lineNo)
  of "from":
    let start = parseF(toks[1], lineNo, "scheduled action start")
    require(start >= 0.0, lineNo, "from START must be >= 0")
    if toks.len >= 7 and toks[2] == "repeat" and toks[4] == "every":
      let count = parseI64(toks[3], lineNo, "repeat count")
      require(count > 0, lineNo, "repeat COUNT must be > 0")
      let period = parseF(toks[5], lineNo, "scheduled action period")
      require(period > 0.0, lineNo, "every PERIOD must be > 0")
      let action = parseActionKind(toks[6], lineNo)
      if action == eaStop: raise newException(ValueError, "line " & $lineNo & ": stop is only valid as a when action")
      let args = if toks.len > 7: toks[7 .. ^1] else: @[]
      validateActionArgs(action, args, lineNo)
      m.scheduledActions.add ScheduledAction(startTime: start, nextTime: start, period: period, repeat: true, remaining: count, active: true, action: action, args: args, lineNo: lineNo)
    else:
      require(toks.len >= 5 and toks[2] == "step", lineNo,
        "from syntax: from START step PERIOD ACTION ... OR from START repeat COUNT every PERIOD ACTION ...")
      let period = parseF(toks[3], lineNo, "scheduled action period")
      require(period > 0.0, lineNo, "from step PERIOD must be > 0")
      let action = parseActionKind(toks[4], lineNo)
      if action == eaStop: raise newException(ValueError, "line " & $lineNo & ": stop is only valid as a when action")
      let args = if toks.len > 5: toks[5 .. ^1] else: @[]
      validateActionArgs(action, args, lineNo)
      m.scheduledActions.add ScheduledAction(startTime: start, nextTime: start, period: period, repeat: true, remaining: -1, active: true, action: action, args: args, lineNo: lineNo)
  of "at":
    let t = parseF(toks[1], lineNo, "scheduled action numeric argument")
    require(t >= 0.0, lineNo, "at TIME must be >= 0")
    let action = parseActionKind(toks[2], lineNo)
    if action == eaStop:
      raise newException(ValueError, "line " & $lineNo & ": stop is only valid as a when action")
    let args = if toks.len > 3: toks[3 .. ^1] else: @[]
    validateActionArgs(action, args, lineNo)
    m.scheduledActions.add ScheduledAction(
      startTime: t, nextTime: t, period: 0.0, repeat: false, remaining: 1, active: true,
      action: action, args: args, lineNo: lineNo)
  else:
    raise newException(ValueError, "line " & $lineNo & ": internal parser error in scheduled action")

proc parseConditionalActionLine(m: var Model; raw: RawLine) =
  let toks = tokenize(raw.text, raw.lineNo)
  require(toks.len >= 5, raw.lineNo,
    "when syntax: when CONDITION [and CONDITION ...] ACTION ...")
  require(toks[0] == "when", raw.lineNo, "internal parser error: expected when")

  var pos = 1
  var conditions: seq[AtomicCondition] = @[]
  while true:
    var observable: ConditionalObservableKind
    var targetId = -1
    require(pos < toks.len, raw.lineNo, "missing condition in when")
    if toks[pos] == "X":
      if pos + 1 < toks.len and toks[pos + 1] notin [">", "<"]:
        targetId = m.monomerId(toks[pos + 1], raw.lineNo)
        observable = woMonomerConversion
        pos += 2
      else:
        observable = woTotalConversion
        pos += 1
    elif toks[pos] == "c":
      require(pos + 1 < toks.len, raw.lineNo, "when c requires a species or monomer name")
      let target = toks[pos + 1]
      if m.monomerByName.hasKey(target):
        observable = woMonomerConc
        targetId = m.monomerByName[target]
      else:
        observable = woSpeciesConc
        targetId = m.speciesId(target, raw.lineNo)
      pos += 2
    else:
      raise newException(ValueError, "line " & $raw.lineNo & ": when condition must start with X or c NAME")

    require(pos + 1 < toks.len, raw.lineNo, "incomplete when condition")
    var comparison: ComparisonKind
    case toks[pos]
    of ">": comparison = coGreater
    of "<": comparison = coLess
    else: raise newException(ValueError, "line " & $raw.lineNo & ": when supports only > or <")
    inc pos
    let threshold = parseF(toks[pos], raw.lineNo, "when threshold")
    inc pos
    conditions.add AtomicCondition(observable: observable, targetId: targetId,
                                   comparison: comparison, threshold: threshold)
    if pos < toks.len and toks[pos] == "and":
      inc pos
      continue
    break

  require(pos < toks.len, raw.lineNo, "missing action in when")
  let action = parseActionKind(toks[pos], raw.lineNo)
  inc pos
  let args = if pos < toks.len: toks[pos .. ^1] else: @[]
  validateActionArgs(action, args, raw.lineNo)
  let first = conditions[0]
  m.conditionalActions.add ConditionalAction(
    active: true, conditions: conditions,
    observable: first.observable, targetId: first.targetId,
    comparison: first.comparison, threshold: first.threshold,
    action: action, args: args, lineNo: raw.lineNo)

proc normalizeSideTokens(toks: seq[string]): seq[string] =
  ## Split compact A+B forms while keeping classic spaced A + B syntax.
  for t in toks:
    if t == "+":
      result.add t
    elif t.contains("+"):
      let parts = t.split("+")
      for i, part in parts:
        if part.len > 0:
          result.add part
        if i < parts.len - 1:
          result.add "+"
    else:
      result.add t

proc parseStoichToken(tok: string; lineNo: int): tuple[name: string, stoich: int] =
  var i = 0
  while i < tok.len and tok[i].isDigit:
    inc i
  if i > 0:
    result.stoich = int(parseI64(tok[0 ..< i], lineNo, "stoichiometric coefficient"))
    result.name = tok[i .. ^1]
  else:
    result.stoich = 1
    result.name = tok
  require(result.stoich > 0, lineNo, "stoichiometric coefficient must be > 0")
  require(result.name.len > 0, lineNo, "missing species name in reaction side")

proc smallKindAndId(m: Model; name: string; lineNo: int): tuple[kind: SmallKind, id: int] =
  if m.speciesByName.hasKey(name):
    return (kind: skSpecies, id: m.speciesByName[name])
  if m.monomerByName.hasKey(name):
    return (kind: skMonomer, id: m.monomerByName[name])
  if m.poolByName.hasKey(name):
    parseFail(lineNo, name & " is a polymer pool; use macro, not rxn")
  parseFail(lineNo, "unknown species or monomer in rxn: " & name)

proc allCharsAreDigits(s: string): bool =
  if s.len == 0: return false
  for ch in s:
    if not ch.isDigit: return false
  result = true

proc parseSmallSide(m: Model; sideToks0: seq[string]; lineNo: int; allowEmpty: bool): seq[SmallRef] =
  var sideToks = normalizeSideTokens(sideToks0)
  # Item 33: a stoichiometric coefficient may be written as its own token
  # before the species name ("2 R0_") as well as concatenated ("2R0_") --
  # arbitrary whitespace between tokens must not change meaning. Merge a
  # standalone all-digit token with the following non-"+" token before
  # parseStoichToken runs, so both spellings parse identically.
  block:
    var merged: seq[string] = @[]
    var i = 0
    while i < sideToks.len:
      if sideToks[i].len > 0 and allCharsAreDigits(sideToks[i]) and
         i + 1 < sideToks.len and sideToks[i + 1] != "+":
        merged.add sideToks[i] & sideToks[i + 1]
        i += 2
      else:
        merged.add sideToks[i]
        i += 1
    sideToks = merged
  if sideToks.len == 0:
    require(allowEmpty, lineNo, "reaction side cannot be empty")
    return
  var expectItem = true
  for tok in sideToks:
    if tok == "+":
      require(not expectItem, lineNo, "unexpected + in reaction side")
      expectItem = true
    else:
      require(expectItem, lineNo, "missing + between reaction species")
      let st = parseStoichToken(tok, lineNo)
      let sid = m.smallKindAndId(st.name, lineNo)
      result.add SmallRef(kind: sid.kind, id: sid.id, stoich: st.stoich)
      expectItem = false
  require(not expectItem, lineNo, "reaction side ends with +")

proc addDelta(deltas: var seq[SmallDelta]; kind: SmallKind; id: int; delta: int) =
  if delta == 0: return
  for i in 0 ..< deltas.len:
    if deltas[i].kind == kind and deltas[i].id == id:
      deltas[i].delta += delta
      return
  deltas.add SmallDelta(kind: kind, id: id, delta: delta)

proc tokenIndex(toks: seq[string]; target: string): int =
  result = -1
  for i, t in toks:
    if t == target:
      return i

proc sliceTokens(toks: seq[string]; a, b: int): seq[string] =
  result = @[]
  if b <= a:
    return
  for i in a ..< b:
    result.add toks[i]

proc parseRxn(m: var Model; toks: seq[string]; lineNo: int) =
  ## Elementary small-molecule reaction. Species and monomers are allowed;
  ## polymer pools are deliberately rejected and must use macro reactions.
  ## Syntax: rxn lhs -> rhs rate [eff]
  let arrow = tokenIndex(toks, "->")
  require(arrow >= 2, lineNo, "rxn syntax: rxn lhs -> rhs rate [eff]")
  require(arrow < toks.len - 2, lineNo, "missing rate in rxn")
  var ratePos = toks.len - 1
  var eff = 1.0
  if toks.len - arrow >= 4:
    # Last token can be an optional efficiency only if the preceding token is a known rate.
    if m.rateByName.hasKey(toks[^2]):
      ratePos = toks.len - 2
      eff = parseF(toks[^1], lineNo, "rxn efficiency")
      require(eff >= 0.0 and eff <= 1.0, lineNo, "rxn efficiency must be in [0,1]")
  let rateName = toks[ratePos]
  let lhs = m.parseSmallSide(sliceTokens(toks, 1, arrow), lineNo, false)
  let rhs = m.parseSmallSide(sliceTokens(toks, arrow + 1, ratePos), lineNo, true)
  require(lhs.len > 0, lineNo, "rxn lhs cannot be empty")

  var deltas: seq[SmallDelta] = @[]
  for r in lhs:
    deltas.addDelta(r.kind, r.id, -r.stoich)
  for r in rhs:
    deltas.addDelta(r.kind, r.id, r.stoich)

  var kind: ChannelKind
  if lhs.len == 1 and lhs[0].stoich == 1:
    kind = chRxnUni
  elif lhs.len == 1 and lhs[0].stoich == 2:
    kind = chRxnBiSame
  elif lhs.len == 2 and lhs[0].stoich == 1 and lhs[1].stoich == 1:
    if lhs[0].kind == lhs[1].kind and lhs[0].id == lhs[1].id:
      kind = chRxnBiSame
    else:
      kind = chRxnBiDiff
  else:
    parseFail(lineNo, "rxn supports only A -> ..., A+B -> ..., or 2A -> ...")

  m.channels.add KmcChannel(
    name: "rxn_" & $m.channels.len, lineNo: lineNo, kind: kind, kId: m.rateId(rateName, lineNo),
    smallReactants: lhs, smallProducts: rhs, smallDeltas: deltas, efficiency: eff)

proc parseMacro(m: var Model; toks: seq[string]; lineNo: int) =
  require(toks.len >= 2, lineNo, "macro kind missing")
  let kind = toks[1].toLowerAscii()

  case kind
  of "init":
    # macro init R + A -> PA ki_a
    require(toks.len == 8 and toks[3] == "+" and toks[5] == "->", lineNo,
      "usage: macro init R + A -> PA k")
    let sp = m.speciesId(toks[2], lineNo)
    let mon = m.monomerId(toks[4], lineNo)
    let outp = m.poolId(toks[6], lineNo)
    require(m.pools[outp].kind == pkActive, lineNo, "init output pool must be active")
    m.setPoolTerminal(outp, mon, lineNo)
    m.setPoolPenultimate(outp, -1, lineNo)
    m.channels.add KmcChannel(
      name: "init_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: chMacroInit, kId: m.rateId(toks[7], lineNo),
      speciesId: sp, monomerId: mon, poolOut: outp
    )

  of "prop":
    # macro prop PA + B -> PB kp_ab
    require(toks.len == 8 and toks[3] == "+" and toks[5] == "->", lineNo,
      "usage: macro prop PA + B -> PB k")
    let inp = m.poolId(toks[2], lineNo)
    let mon = m.monomerId(toks[4], lineNo)
    let outp = m.poolId(toks[6], lineNo)
    require(m.pools[inp].kind == pkActive, lineNo, "prop input pool must be active")
    require(m.pools[outp].kind == pkActive, lineNo, "prop output pool must be active")
    m.setPoolTerminal(outp, mon, lineNo)
    if inp >= 0 and inp < m.poolTerminalMer.len:
      m.setPoolPenultimate(outp, m.poolTerminalMer[inp], lineNo)
    m.channels.add KmcChannel(
      name: "prop_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: chMacroProp, kId: m.rateId(toks[7], lineNo),
      pool1: inp, monomerId: mon,
      poolOut: outp
    )

  of "term_c", "term_d":
    # macro term_c PA + PB -> D ktc
    # macro term_d PA + PB -> D ktd
    # macro term_d PA + PB -> D + D ktd   # classic slimmc-compatible alias
    var outName: string
    var rateName: string
    if kind == "term_d" and toks.len == 10:
      require(toks[3] == "+" and toks[5] == "->" and toks[7] == "+", lineNo,
        "usage: macro term_d PA + PB -> D k OR macro term_d PA + PB -> D + D k")
      require(toks[6] == toks[8], lineNo, "term_d D + D alias requires the same dead pool on both sides")
      outName = toks[6]
      rateName = toks[9]
    else:
      require(toks.len == 8 and toks[3] == "+" and toks[5] == "->", lineNo,
        "usage: macro term_c PA + PB -> D k OR macro term_d PA + PB -> D [ + D ] k")
      outName = toks[6]
      rateName = toks[7]
    let ck = if kind == "term_c": chMacroTermC else: chMacroTermD
    m.channels.add KmcChannel(
      name: kind & "_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: ck, kId: m.rateId(rateName, lineNo),
      pool1: m.poolId(toks[2], lineNo), pool2: m.poolId(toks[4], lineNo),
      poolOut: m.poolId(outName, lineNo)
    )

  of "term_x":
    # macro term_x PA + Cap -> D ktx
    require(toks.len == 8 and toks[3] == "+" and toks[5] == "->", lineNo,
      "usage: macro term_x PA + Cap -> D k")
    let inp = m.poolId(toks[2], lineNo)
    let cap = m.speciesId(toks[4], lineNo)
    let outp = m.poolId(toks[6], lineNo)
    require(m.pools[inp].kind == pkActive, lineNo, "term_x input pool must be active")
    require(m.pools[outp].kind == pkDead, lineNo, "term_x output pool must be dead")
    m.channels.add KmcChannel(
      name: "term_x_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: chMacroTermX, kId: m.rateId(toks[7], lineNo),
      pool1: inp, speciesId: cap, poolOut: outp
    )

  of "transfer":
    # macro transfer PA + CTA -> D + Rcta ktr_a
    # "transfer_h" used to be accepted as an alias here (same
    # chMacroTransfer channel, different spelling only) -- removed per
    # item 55's action/macro-alias cleanup, matching homo's identical
    # removal.
    let macroName = kind
    require(toks.len == 10 and toks[3] == "+" and toks[5] == "->" and toks[7] == "+", lineNo,
      "usage: macro " & macroName & " PA + CTA -> D + Rcta k")
    let inp = m.poolId(toks[2], lineNo)
    let acceptor = m.speciesId(toks[4], lineNo)
    let outp = m.poolId(toks[6], lineNo)
    let radical = m.speciesId(toks[8], lineNo)
    require(m.pools[inp].kind == pkActive, lineNo, macroName & " input pool must be active")
    require(m.pools[outp].kind == pkDead, lineNo, macroName & " output pool must be dead")
    m.channels.add KmcChannel(
      name: "transfer_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: chMacroTransfer, kId: m.rateId(toks[9], lineNo),
      pool1: inp, speciesId: acceptor, poolOut: outp, speciesOutId: radical
    )

  of "transfer_m":
    # macro transfer_m PA + A -> D + PA ktrm_a
    require(toks.len == 10 and toks[3] == "+" and toks[5] == "->" and toks[7] == "+", lineNo,
      "usage: macro transfer_m PA + A -> D + PA k")
    let inp = m.poolId(toks[2], lineNo)
    let mon = m.monomerId(toks[4], lineNo)
    let outp = m.poolId(toks[6], lineNo)
    let newp = m.poolId(toks[8], lineNo)
    require(m.pools[inp].kind == pkActive, lineNo, "transfer_m input pool must be active")
    require(m.pools[outp].kind == pkDead, lineNo, "transfer_m dead output pool must be dead")
    require(m.pools[newp].kind == pkActive, lineNo, "transfer_m active product pool must be active")
    m.setPoolTerminal(newp, mon, lineNo)
    m.setPoolPenultimate(newp, -1, lineNo)
    m.channels.add KmcChannel(
      name: "transfer_m_" & toks[2] & "_" & toks[4], lineNo: lineNo,
      kind: chMacroTransferM, kId: m.rateId(toks[9], lineNo),
      pool1: inp, monomerId: mon, poolOut: outp, pool2: newp
    )

  of "deprop":
    # macro deprop PA -> PB + A kdeprop_BA
    require(toks.len == 8 and toks[3] == "->" and toks[5] == "+", lineNo,
      "usage: macro deprop PA -> PB + A k")
    let inp = m.poolId(toks[2], lineNo)
    let outp = m.poolId(toks[4], lineNo)
    let mon = m.monomerId(toks[6], lineNo)
    require(m.pools[inp].kind == pkActive, lineNo, "deprop input pool must be active")
    require(m.pools[outp].kind == pkActive, lineNo, "deprop output pool must be active")
    m.setPoolTerminal(inp, mon, lineNo)
    m.channels.add KmcChannel(
      name: "deprop_" & toks[2] & "_" & toks[6], lineNo: lineNo,
      kind: chMacroDeprop, kId: m.rateId(toks[7], lineNo),
      pool1: inp, monomerId: mon, poolOut: outp
    )

  else:
    raise newException(ValueError, "line " & $lineNo & ": unsupported macro reaction in v0.6: " & kind)


proc inferPoolMetadata(m: var Model) =
  ## Re-run terminal/penultimate inference after all channels have been parsed.
  ## This removes a hidden ordering assumption for penultimate propagation models:
  ## an output pool may be referenced before the parser has seen the channel that
  ## defines the input pool terminal metadata.
  var changed = true
  while changed:
    changed = false
    for ch in m.channels:
      case ch.kind
      of chMacroInit:
        let oldT = m.poolTerminalMer[ch.poolOut]
        m.setPoolTerminal(ch.poolOut, ch.monomerId, ch.lineNo)
        if m.poolTerminalMer[ch.poolOut] != oldT: changed = true
      of chMacroProp:
        let oldT = m.poolTerminalMer[ch.poolOut]
        let oldP = m.poolPenultimateMer[ch.poolOut]
        m.setPoolTerminal(ch.poolOut, ch.monomerId, ch.lineNo)
        if ch.pool1 >= 0 and ch.pool1 < m.poolTerminalMer.len:
          let prevTerminal = m.poolTerminalMer[ch.pool1]
          if prevTerminal >= 0:
            m.setPoolPenultimate(ch.poolOut, prevTerminal, ch.lineNo)
        if m.poolTerminalMer[ch.poolOut] != oldT or m.poolPenultimateMer[ch.poolOut] != oldP:
          changed = true
      of chMacroTransferM:
        let oldT = m.poolTerminalMer[ch.pool2]
        m.setPoolTerminal(ch.pool2, ch.monomerId, ch.lineNo)
        if m.poolTerminalMer[ch.pool2] != oldT: changed = true
      of chMacroDeprop:
        let oldT = m.poolTerminalMer[ch.pool1]
        m.setPoolTerminal(ch.pool1, ch.monomerId, ch.lineNo)
        if m.poolTerminalMer[ch.pool1] != oldT: changed = true
      else:
        discard

proc initDefaults(m: var Model) =
  # V and t_end deliberately default to invalid sentinels (item 44:
  # volume and t_end are mandatory parameters) -- matches homo's
  # approach exactly. Previously defaulted to 1.0e-18/10.0, silently
  # letting a model with neither "param kmc_volume" nor "param t_end" run
  # to completion with arbitrary values instead of erroring.
  m.V = 0.0
  m.initVolumeMl = 0.0
  m.currentVolumeMl = 0.0
  m.hasInitVolume = false
  m.T = DefaultTemperatureK
  m.t_end = -1.0
  m.max_steps = DefaultMaxSteps
  m.whenCheckEvents = DefaultWhenCheckEvents
  m.seed = DefaultSeed
  m.output_dir = ""
  m.output_dirWasSet = false
  m.description = ""
  m.hasDescription = false
  m.variables = @[]
  m.dp_max = DefaultDpMax
  m.sequence_mode = DefaultSequenceMode
  m.mass_model = mmRepeatUnits
  m.deadPoolId = -1
  m.monomerByName = initTable[string, int]()
  m.speciesByName = initTable[string, int]()
  m.feeds = @[]
  m.feedByName = initTable[string, int]()
  m.endgroupByName = initTable[string, int]()
  m.poolByName = initTable[string, int]()
  m.poolTerminalMer = @[]
  m.poolPenultimateMer = @[]
  m.rateByName = initTable[string, int]()
  m.scheduledActions = @[]
  m.conditionalActions = @[]
  m.rawConditionalActions = @[]

proc resolveVars(m: var Model) =
  for i in 0 ..< m.variables.len:
    case m.variables[i].kind
    of "rate":
      let rid = m.rateId(m.variables[i].name, 0)
      m.variables[i].value = m.rateValue(rid)
    of "param":
      case m.variables[i].name
      of "kmc_volume": m.variables[i].value = m.V
      of "init_volume": m.variables[i].value = m.initVolumeMl
      of "temperature": m.variables[i].value = m.T
      of "t_end": m.variables[i].value = m.t_end
      of "max_steps": m.variables[i].value = float(m.max_steps)
      of "when_check_events": m.variables[i].value = float(m.whenCheckEvents)
      of "seed": m.variables[i].value = float(m.seed)
      of "dp_max": m.variables[i].value = float(m.dp_max)
      else:
        raise newException(ValueError, "var param target is not a numeric built-in parameter: " & m.variables[i].name)
    of "species":
      let sid = m.speciesId(m.variables[i].name, 0)
      m.variables[i].value = m.species[sid].c0
    of "monomer":
      let mid = m.monomerId(m.variables[i].name, 0)
      m.variables[i].value = m.monomers[mid].c0

    of "endgroup":
      require(m.endgroupByName.hasKey(m.variables[i].name), 0, "var endgroup target not found: " & m.variables[i].name)
      let eid = m.endgroupByName[m.variables[i].name]
      m.variables[i].value = m.endgroups[eid].mw
    else:
      raise newException(ValueError, "unknown var kind: " & m.variables[i].kind)


proc parseProcessVolumeMl(toks: seq[string]; valueIndex: int; lineNo: int; context: string): float =
  ## Process volumes default to litres. Explicit L/l and mL/ml/ML are accepted.
  require(toks.len in [valueIndex + 1, valueIndex + 2], lineNo, context & " syntax: VALUE [L|l|mL|ml|ML]")
  let value = parseF(toks[valueIndex], lineNo, context)
  require(value > 0.0, lineNo, context & " must be > 0")
  if toks.len == valueIndex + 1:
    return value * 1000.0
  let unit = toks[valueIndex + 1].toLowerAscii()
  case unit
  of "l": value * 1000.0
  of "ml": value
  else:
    require(false, lineNo, context & " unit must be L or mL (case-insensitive); without a unit the value is in L")
    0.0

proc parseModel*(filename: string): Model =
  result.initDefaults()
  result.modelFile = filename
  result.modelStem = validateModelRunId(filename)
  let lines = readFile(filename).splitLines()
  for lineNo0, raw in lines:
    let lineNo = lineNo0 + 1
    let line = stripComment(raw)
    if line.len == 0: continue
    result.rawLines.add RawLine(text: raw, lineNo: lineNo)
    let toks = tokenize(line, lineNo)
    if toks.len == 0: continue
    let key = toks[0]

    case key
    of "param":
      require(toks.len in [3, 4], lineNo, "usage: param NAME VALUE [unit]")
      let p = toks[1]
      case p
      of "volume":
        raise newException(ValueError, "line " & $lineNo & ": param volume was renamed to param kmc_volume; replace: param volume VALUE -> param kmc_volume VALUE")
      of "kmc_volume":
        require(toks.len == 3, lineNo, "param kmc_volume does not accept a unit; its unit is L")
        result.V = parseF(toks[2], lineNo, "numeric value")
      of "init_volume":
        result.initVolumeMl = parseProcessVolumeMl(toks, 2, lineNo, "parameter init_volume")
        result.hasInitVolume = true
      of "temperature":
        require(toks.len == 3, lineNo, "only param init_volume accepts an explicit unit (L or mL, case-insensitive)")
        result.T = parseF(toks[2], lineNo, "numeric value")
      of "t_end": result.t_end = parseF(toks[2], lineNo, "numeric value")
      of "max_steps": result.max_steps = parseI64(toks[2], lineNo, "integer parameter")
      of "when_check_events": result.whenCheckEvents = parseI64(toks[2], lineNo, "integer parameter")
      of "seed": result.seed = parseI64(toks[2], lineNo, "integer parameter")
      of "output_dir":
        block:
          let idx = raw.find("output_dir")
          require(idx >= 0, lineNo, "param output_dir syntax: param output_dir \"PATH\"")
          let afterKeyword = raw[idx + "output_dir".len .. ^1].strip(leading = true, trailing = false)
          let wasQuoted = afterKeyword.len >= 2 and (afterKeyword[0] == '"' or afterKeyword[0] == '\'')
          require(wasQuoted, lineNo, "output_dir must be quoted: param output_dir \"results/run_01\"")
        require(isValidModelPath(toks[2]), lineNo,
          "invalid output_dir path: each path segment must match " & ModelIdentifierPattern)
        result.output_dir = toks[2]
        result.output_dirWasSet = true
      of "dp_max":
        result.dp_max = parseI64(toks[2], lineNo, "integer parameter dp_max")
        require(result.dp_max > 0 and result.dp_max <= int64(high(int32)), lineNo, "dp_max must be in 1..high(int32)")
      of "sequence_mode":
        require(toks[2] in ["composition", "full"], lineNo, "sequence_mode must be composition or full")
        result.sequence_mode = toks[2]
      of "mass_model": result.mass_model = parseMassModelValue(toks[2], lineNo)
      else:
        raise newException(ValueError, "line " & $lineNo & ": unknown param: " & p & " (only canonical snake_case parameters are supported)")

    of "desc":
      require(not result.hasDescription, lineNo, "desc may appear only once")
      let descStripped = line.strip()
      let quotedOk = descStripped.startsWith("desc \"") or descStripped.startsWith("desc '")
      require(toks.len == 2 and quotedOk, lineNo, "desc syntax: desc \"text\"")
      result.description = toks[1]
      result.hasDescription = true

    of "var":
      require(toks.len == 4 and toks[1] in ["rate", "param", "species", "monomer", "endgroup"],
        lineNo, "usage: var rate|param|species|monomer|endgroup NAME UNIT")
      for v in result.variables:
        require(v.name != toks[2], lineNo, "var target name may be declared only once: " & toks[2])
      result.variables.add VarDef(kind: toks[1], name: toks[2], value: 0.0, unit: toks[3])

    of "monomer":
      require(toks.len == 4, lineNo, "usage: monomer NAME c0 MW")
      result.addMonomer(toks[1], parseF(toks[2], lineNo, "numeric value"), parseF(toks[3], lineNo, "numeric value"), lineNo)

    of "species":
      require(toks.len == 3, lineNo, "usage: species NAME c0")
      result.addSpecies(toks[1], parseF(toks[2], lineNo, "numeric value"), lineNo)

    of "feed":
      require(toks.len == 4, lineNo, "usage: feed NAME SPECIES_OR_MONOMER CONCENTRATION")
      let fname = toks[1]
      require(isValidModelIdentifier(fname), lineNo, "invalid feed name: " & fname)
      let c = parseF(toks[3], lineNo, "feed concentration")
      require(c >= 0.0, lineNo, "feed concentration must be >= 0")
      var fid: int
      if result.feedByName.hasKey(fname):
        fid = result.feedByName[fname]
      else:
        fid = result.feeds.len
        result.feedByName[fname] = fid
        result.feeds.add FeedDef(name: fname, monomerConcentrations: newSeq[float](result.monomers.len), speciesConcentrations: newSeq[float](result.species.len))
      if result.feeds[fid].monomerConcentrations.len < result.monomers.len: result.feeds[fid].monomerConcentrations.setLen(result.monomers.len)
      if result.feeds[fid].speciesConcentrations.len < result.species.len: result.feeds[fid].speciesConcentrations.setLen(result.species.len)
      if result.monomerByName.hasKey(toks[2]):
        let mid = result.monomerByName[toks[2]]
        require(result.feeds[fid].monomerConcentrations[mid] == 0.0, lineNo, "feed component already declared")
        result.feeds[fid].monomerConcentrations[mid] = c
      elif result.speciesByName.hasKey(toks[2]):
        let sid = result.speciesByName[toks[2]]
        require(result.feeds[fid].speciesConcentrations[sid] == 0.0, lineNo, "feed component already declared")
        result.feeds[fid].speciesConcentrations[sid] = c
      else:
        raise newException(ValueError, "line " & $lineNo & ": unknown species/monomer in feed: " & toks[2])



    of "endgroup":
      require(toks.len == 3, lineNo, "usage: endgroup NAME MW")
      result.addEndGroup(toks[1], parseF(toks[2], lineNo, "numeric value"), lineNo)

    of "polymer":
      require(toks.len == 3, lineNo, "usage: polymer NAME active|dead")
      let kind = case toks[2]
        of "active": pkActive
        of "dead": pkDead
        else: raise newException(ValueError, "line " & $lineNo & ": polymer kind must be active/dead")
      result.addPool(toks[1], kind, lineNo)

    of "rate":
      require(toks.len == 3 or toks.len == 4 or toks.len == 5, lineNo, "usage: rate NAME VALUE, rate NAME const VALUE, or rate NAME arr Apre Ea")
      if toks.len == 3:
        result.addFixedRate(toks[1], parseF(toks[2], lineNo, "numeric value"), lineNo)
      elif toks.len == 4 and toks[2] == "const":
        result.addFixedRate(toks[1], parseF(toks[3], lineNo, "numeric value"), lineNo)
      elif toks.len == 5 and toks[2] == "arr":
        result.addArrheniusRate(toks[1], parseF(toks[3], lineNo, "numeric value"), parseF(toks[4], lineNo, "numeric value"), lineNo)
      else:
        raise newException(ValueError, "line " & $lineNo & ": usage: rate NAME VALUE, rate NAME const VALUE, or rate NAME arr Apre Ea")

    of "rxn":
      result.parseRxn(toks, lineNo)

    of "macro":
      result.parseMacro(toks, lineNo)

    of "every", "from", "at":
      result.parseScheduledAction(toks, lineNo)

    of "when":
      result.rawConditionalActions.add RawLine(text: line, lineNo: lineNo)

    of "memory_limit":
      require(toks.len == 2, lineNo, "usage: memory_limit 50GB")
      result.memoryPolicy.hasLimit = true
      result.memoryPolicy.limitBytes = parseSizeBytes(toks[1])

    of "at_memory":
      # Items 52/53/39: unknown tokens must raise (previously silently
      # discarded); "compact_dead"/"drop_dead_seq" removed entirely
      # (permanently inert no-ops, never implemented -- copo already
      # always stores DEAD summaries); matching is case-sensitive
      # (previously .toLowerAscii()).
      require(toks.len == 3 or toks.len == 4, lineNo, "usage: at_memory 50GB save|stop")
      result.memoryPolicy.hasLimit = true
      result.memoryPolicy.limitBytes = parseSizeBytes(toks[1])
      for a in toks[2 .. ^1]:
        case a
        of "save": result.memoryPolicy.snapshotOnLimit = true
        of "stop": result.memoryPolicy.stopOnLimit = true
        else: raise newException(ValueError, "line " & $lineNo & ": at_memory: unknown token: " & a & " (expected save and/or stop)")

    else:
      raise newException(ValueError, "line " & $lineNo & ": unknown directive: " & toks[0])

  require(result.V > 0, 0, "kmc_volume must be > 0")
  require(result.T > 0, 0, "temperature must be > 0")
  require(result.t_end >= 0, 0, "t_end must be >= 0")
  require(result.max_steps > 0, 0, "max_steps must be > 0")
  require(result.monomers.len >= 2 and result.monomers.len <= 3, 0,
    "v0.3 expects 2 or 3 monomers")
  require(result.deadPoolId >= 0, 0, "one dead polymer pool is required, e.g. polymer D dead")
  require(result.channels.len > 0, 0, "no reaction channels defined")
  require(result.whenCheckEvents > 0, 0, "when_check_events must be > 0")
  result.inferPoolMetadata()
  for ch in result.channels:
    if ch.kind == chMacroDeprop:
      require(result.sequence_mode == "full", ch.lineNo,
        "sequence_mode=composition cannot be used with copolymer depropagation; set sequence_mode=full or remove depropagation channels")
      require(result.poolTerminalMer[ch.poolOut] >= 0, ch.lineNo,
        "deprop output pool terminal is unknown; define init/prop channels for it before use")
  for raw in result.rawConditionalActions:
    result.parseConditionalActionLine(raw)
  if result.feeds.len > 0:
    require(result.hasInitVolume and result.initVolumeMl > 0.0, 0, "param init_volume must be > 0 when feed is used")
  elif result.hasInitVolume:
    require(result.initVolumeMl > 0.0, 0, "init_volume must be > 0")
  result.currentVolumeMl = result.initVolumeMl
  result.resolveVars()
  if not result.output_dirWasSet:
    result.output_dir = defaultOutputDir(filename)
  else:
    # Parity fix (see audit report, "output_dir CWD-vs-file-relative"
    # finding): an explicit relative `param output_dir X` used to be used
    # completely as-is, resolved by the OS against the process's current
    # working directory at runtime. Classic slimmc instead resolves a
    # relative output_dir against the .model file's own directory
    # (slimmc_parser.nim: parseOutputDir), which is predictable regardless
    # of the caller's CWD. Unified here to match slimmc's behavior.
    if not result.output_dir.isAbsolute:
      let parent = parentDir(filename)
      if parent.len > 0 and parent != ".":
        result.output_dir = parent / result.output_dir
