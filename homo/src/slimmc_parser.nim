# slimmc_parser.nim
# Parser and model validation for slimmc.

import strutils, tables, os
import slimmc_types
import ../../common/run_id
import ../../common/model_contract

proc isFloatToken(s: string): bool =
  try:
    discard parseFloat(s)
    result = true
  except ValueError:
    result = false

proc validateIdentifier(name: string; lineNo: int; what: string) =
  if not isValidModelIdentifier(name):
    fail(lineNo, what & " must match " & ModelIdentifierPattern & ": " & name)
  if isReservedModelIdentifier(name):
    fail(lineNo, what & " is reserved: " & name)

proc findDeclaredNameCaseInsensitive(m: Model; name: string): string =
  let key = modelNameKey(name)
  for existing in m.speciesByName.keys:
    if modelNameKey(existing) == key: return existing
  for existing in m.poolByName.keys:
    if modelNameKey(existing) == key: return existing
  for existing in m.rateByName.keys:
    if modelNameKey(existing) == key: return existing
  ""

proc requireGlobalNameFree(m: Model; name: string; lineNo: int; what: string) =
  let existing = m.findDeclaredNameCaseInsensitive(name)
  if existing.len > 0:
    fail(lineNo, what & " name conflicts with already declared model name: " & existing)

proc indexOfToken(tokens: seq[string]; val: string): int =
  result = -1
  for i, t in tokens:
    if t == val:
      return i

proc sliceTokens(tokens: seq[string]; a, b: int): seq[string] =
  result = @[]
  if b <= a:
    return
  for i in a ..< b:
    result.add tokens[i]

proc tokenize*(line: string; lineNo: int = 0): seq[string] =
  ## Splits a model line, keeping quoted strings as one token. Both
  ## single and double quotes are accepted as string delimiters (item
  ## 36) -- whichever one opens a string, the same one must close it
  ## (quoteChar tracks which). Comments start at # only outside quoted
  ## strings. The -> token is separated.
  var tokens: seq[string] = @[]
  var current = ""
  var inQuote = false
  var quoteChar = '\0'
  var escaped = false
  var i = 0

  proc flush() =
    if current.len > 0:
      tokens.add current
      current = ""

  while i < line.len:
    let ch = line[i]

    if inQuote:
      if escaped:
        # Common escape sequences in both string kinds (item 38): \n,
        # \r, \t, \\, \", \'.
        case ch
        of 'n': current.add '\n'
        of 'r': current.add '\r'
        of 't': current.add '\t'
        of '"': current.add '"'
        of '\'': current.add '\''
        of '\\': current.add '\\'
        else: current.add ch
        escaped = false
      elif ch == '\\':
        escaped = true
      elif ch == quoteChar:
        inQuote = false
        flush()
      else:
        current.add ch
      inc i
      continue

    if ch == '#':
      break
    elif ch == '"' or ch == '\'':
      if current.len > 0:
        fail(lineNo, "quoted text must start after whitespace")
      inQuote = true
      quoteChar = ch
    elif ch.isSpaceAscii:
      flush()
    elif ch == '-' and i + 1 < line.len and line[i + 1] == '>':
      flush()
      tokens.add "->"
      inc i
    else:
      current.add ch
    inc i

  if inQuote:
    fail(lineNo, "unterminated quoted string")
  if escaped:
    fail(lineNo, "unterminated escape in quoted string")
  flush()
  result = tokens

proc parseSide(tokens: seq[string]; lineNo: int): seq[RawTerm] =
  result = @[]
  let s = tokens.join("")
  if s.len == 0 or s == "0":
    return

  for part0 in s.split("+"):
    let part = part0.strip()
    if part.len == 0:
      continue

    var i = 0
    while i < part.len and part[i].isDigit:
      inc i

    var st = 1
    var name = part
    if i > 0:
      st = parseInt(part[0 ..< i])
      name = part[i .. ^1]

    if st <= 0:
      fail(lineNo, "invalid stoichiometry in term: " & part)
    if name.len == 0:
      fail(lineNo, "missing name in term: " & part)

    result.add RawTerm(name: name, stoich: st)

proc expandTerms(raw: seq[RawTerm]): seq[string] =
  result = @[]
  for rt in raw:
    for i in 0 ..< rt.stoich:
      result.add rt.name

proc addSpecies(
  m: var Model;
  name: string;
  kind: SpeciesKind;
  c0: float;
  mw: float;
  hasMw: bool;
  lineNo: int
): int =
  validateIdentifier(name, lineNo, "species name")
  m.requireGlobalNameFree(name, lineNo, "species/monomer")

  if m.poolByName.hasKey(name):
    fail(lineNo, name & " is already declared as polymer pool")
  if m.speciesByName.hasKey(name):
    fail(lineNo, name & " is already declared as species/monomer")
  if kind == skMonomer and m.monomerId >= 0:
    fail(lineNo, "only one monomer is supported")

  result = m.species.len
  m.species.add SpeciesDef(name: name, kind: kind, c0: c0, mw: mw, hasMw: hasMw)
  m.speciesByName[name] = result

  if kind == skMonomer:
    m.monomerId = result

proc addPool(m: var Model; name: string; kind: PoolKind; lineNo: int): int =
  validateIdentifier(name, lineNo, "polymer pool name")
  m.requireGlobalNameFree(name, lineNo, "polymer pool")

  if m.speciesByName.hasKey(name):
    fail(lineNo, name & " is already declared as species/monomer")
  if m.poolByName.hasKey(name):
    fail(lineNo, name & " is already declared as polymer pool")

  result = m.pools.len
  m.pools.add PoolDef(name: name, kind: kind)
  m.poolByName[name] = result

proc addRate(m: var Model; r: RateDef; lineNo: int): int =
  validateIdentifier(r.name, lineNo, "rate name")
  m.requireGlobalNameFree(r.name, lineNo, "rate")
  if m.rateByName.hasKey(r.name):
    fail(lineNo, "rate already declared: " & r.name)

  result = m.rates.len
  m.rates.add r
  m.rateByName[r.name] = result

proc getRateId*(m: Model; name: string; lineNo: int): int =
  if not m.rateByName.hasKey(name):
    fail(lineNo, "unknown rate: " & name)
  result = m.rateByName[name]

proc getSpeciesId*(m: Model; name: string; lineNo: int): int =
  if not m.speciesByName.hasKey(name):
    fail(lineNo, "unknown species/monomer: " & name)
  result = m.speciesByName[name]

proc getPoolId*(m: Model; name: string; lineNo: int): int =
  if not m.poolByName.hasKey(name):
    fail(lineNo, "unknown polymer pool: " & name)
  result = m.poolByName[name]

proc requirePoolKind*(m: Model; name: string; kind: PoolKind; lineNo: int): int =
  result = m.getPoolId(name, lineNo)
  if m.pools[result].kind != kind:
    let expected = if kind == pkActive: "active" else: "dead"
    fail(lineNo, "polymer pool " & name & " must be " & expected)

proc parseActionKind(s: string; lineNo: int): ActionKind =
  ## Canonical action names only (item 55) -- "snapshot"/"chains"/"memory"
  ## were accepted as aliases of "save"/"save_chains"/"print_memory" for
  ## parity with an older slimmc-copo spelling; removed per the agreed
  ## contract cleanup (both engines now require the canonical spelling
  ## only -- see copo's parseActionKind for the matching removal).
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
    fail(lineNo, "unknown action: " & s)

proc sourceValueIsQuoted(line: string; keyword: string): bool =
  ## Tokenization deliberately removes quote delimiters. This small source check
  ## enforces the v2.8 rule that free text and output paths are quoted --
  ## either quote kind is accepted (item 36: single and double quotes).
  let stripped = line.strip()
  if not stripped.startsWith(keyword):
    return false
  if stripped.len <= keyword.len:
    return false
  if not stripped[keyword.len].isSpaceAscii:
    return false
  let rest = stripped[keyword.len .. ^1].strip()
  result = rest.len >= 2 and (rest[0] == '"' or rest[0] == '\'')

proc parseDesc(m: var Model; tokens: seq[string]; lineNo: int; source: string) =
  if m.hasDescription:
    fail(lineNo, "desc may be declared only once")
  if tokens.len != 2 or not sourceValueIsQuoted(source, "desc"):
    fail(lineNo, "desc syntax: desc \"text\"")
  m.description = tokens[1]
  m.hasDescription = true

proc parseVar(m: var Model; tokens: seq[string]; lineNo: int) =
  # Canonical syntax: var rate|param|species|monomer|endgroup NAME UNIT.
  if tokens.len == 4 and tokens[1] in ["rate", "param", "species", "monomer", "endgroup"]:
    validateIdentifier(tokens[2], lineNo, "var target name")
    for v in m.variables:
      if v.name == tokens[2]:
        fail(lineNo, "var target name may be declared only once: " & tokens[2])
    m.variables.add VarDef(kind: tokens[1], name: tokens[2], value: 0.0, unit: tokens[3])
    return

  fail(lineNo, "var syntax: var rate|param|species|monomer|endgroup NAME UNIT")

proc parseOutputDir(m: var Model; tokens: seq[string]; lineNo: int; sourcePath: string; source: string) =
  if tokens.len != 3:
    fail(lineNo, "param output_dir syntax: param output_dir \"PATH\"")

  let idx = source.find("output_dir")
  if idx < 0:
    fail(lineNo, "param output_dir syntax: param output_dir \"PATH\"")
  let afterKeyword = source[idx + "output_dir".len .. ^1].strip(leading = true, trailing = false)
  let wasQuoted = afterKeyword.len >= 2 and (afterKeyword[0] == '"' or afterKeyword[0] == '\'')
  if not wasQuoted:
    fail(lineNo, "output_dir must be quoted: param output_dir \"results/run_01\"")

  let raw = tokens[2]
  if not isValidModelPath(raw):
    fail(lineNo, "invalid output_dir path: each path segment must match " & ModelIdentifierPattern)

  if raw.isAbsolute:
    m.outputDir = raw
  else:
    let parent = parentDir(sourcePath)
    m.outputDir = if parent.len == 0 or parent == ".": raw else: parent / raw

proc resolveSeriesVars(m: var Model) =
  for i in 0 ..< m.variables.len:
    case m.variables[i].kind
    of "rate":
      let k = m.getRateId(m.variables[i].name, 0)
      m.variables[i].value = rateValue(m, k)
    of "species", "monomer":
      let sp = m.getSpeciesId(m.variables[i].name, 0)
      if m.variables[i].kind == "monomer" and m.species[sp].kind != skMonomer:
        fail(0, "var monomer target is not a monomer: " & m.variables[i].name)
      if m.variables[i].kind == "species" and m.species[sp].kind != skSpecies:
        fail(0, "var species target is not a species: " & m.variables[i].name)
      m.variables[i].value = m.species[sp].c0
    of "endgroup":
      if not m.egByName.hasKey(m.variables[i].name):
        fail(0, "var endgroup target not found: " & m.variables[i].name)
      let eg = m.egByName[m.variables[i].name]
      m.variables[i].value = m.egMass[eg]
    of "param":
      case m.variables[i].name
      of "kmc_volume": m.variables[i].value = m.V
      of "init_volume": m.variables[i].value = m.initVolumeMl
      of "temperature": m.variables[i].value = m.T
      of "t_end": m.variables[i].value = m.tEnd
      of "max_steps": m.variables[i].value = float(m.maxEvents)
      of "when_check_events": m.variables[i].value = float(m.whenCheckEvents)
      of "seed": m.variables[i].value = float(m.seed)
      of "dp_max": m.variables[i].value = float(m.dpMax)
      else: fail(0, "var param target is not a numeric built-in parameter: " & m.variables[i].name)
    else:
      fail(0, "unknown var kind: " & m.variables[i].kind)


proc parseProcessVolumeMl(tokens: seq[string]; valueIndex: int; lineNo: int; context: string): float =
  ## Process volumes default to litres. Explicit L/l and mL/ml/ML are accepted.
  if tokens.len notin [valueIndex + 1, valueIndex + 2]:
    fail(lineNo, context & " syntax: VALUE [L|l|mL|ml|ML]")
  let value = parseF(tokens[valueIndex], lineNo, context)
  if value <= 0.0:
    fail(lineNo, context & " must be > 0")
  if tokens.len == valueIndex + 1:
    return value * 1000.0
  let unit = tokens[valueIndex + 1].toLowerAscii()
  case unit
  of "l": value * 1000.0
  of "ml": value
  else:
    fail(lineNo, context & " unit must be L or mL (case-insensitive); without a unit the value is in L")
    0.0

proc parseParam(m: var Model; tokens: seq[string]; lineNo: int; source: string) =
  if tokens.len < 3 or tokens.len > 4:
    fail(lineNo, "param requires: param name value [unit]")

  let name = tokens[1]
  let value = tokens[2]
  if name != "init_volume" and tokens.len != 3:
    fail(lineNo, "only param init_volume accepts an explicit unit (L or mL, case-insensitive) (L or mL, case-insensitive)")

  case name
  of "volume":
    fail(lineNo, "param volume was renamed to param kmc_volume; replace: param volume VALUE -> param kmc_volume VALUE")
  of "kmc_volume": m.V = parseF(value, lineNo, "parameter kmc_volume")
  of "init_volume":
    m.initVolumeMl = parseProcessVolumeMl(tokens, 2, lineNo, "parameter init_volume")
    m.hasInitVolume = true
  of "temperature": m.T = parseF(value, lineNo, "parameter temperature")
  of "t_end": m.tEnd = parseF(value, lineNo, "parameter t_end")
  of "max_steps": m.maxEvents = parseI64(value, lineNo, "integer parameter max_steps")
  of "when_check_events": m.whenCheckEvents = parseI64(value, lineNo, "integer parameter when_check_events")
  of "seed": m.seed = parseI64(value, lineNo, "integer parameter seed")
  of "mass_model":
    if value notin ["repeat_units", "with_end_groups"]:
      fail(lineNo, "mass_model must be repeat_units or with_end_groups")
    m.massModel = value
  of "dp_max":
    m.dpMax = parseI64(value, lineNo, "integer parameter dp_max")
    if m.dpMax <= 0:
      fail(lineNo, "dp_max must be > 0")
  of "sequence_mode":
    if value notin ["composition", "full"]:
      fail(lineNo, "sequence_mode must be composition or full")
    m.sequenceMode = value
  of "output_dir": m.parseOutputDir(tokens, lineNo, m.modelSourceFile, source)
  else:
    fail(lineNo, "unknown param: " & name & " (only canonical snake_case parameters are supported)")

proc parseSpeciesDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 3:
    fail(lineNo, "species syntax: species NAME c0; species MW is no longer supported, use endgroup NAME MW for chain-mass contributions")

  let name = tokens[1]
  let c0 = parseF(tokens[2], lineNo, "species concentration")
  if c0 < 0.0:
    fail(lineNo, "species concentration must be >= 0")

  discard m.addSpecies(name, skSpecies, c0, 0.0, false, lineNo)

proc parseMonomerDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 4:
    fail(lineNo, "monomer syntax: monomer NAME c0 Mw")

  let name = tokens[1]
  let c0 = parseF(tokens[2], lineNo, "monomer concentration")
  let mw = parseF(tokens[3], lineNo, "monomer molecular weight")

  if c0 < 0.0:
    fail(lineNo, "monomer concentration must be >= 0")
  if mw <= 0.0:
    fail(lineNo, "monomer molecular weight must be > 0")

  discard m.addSpecies(name, skMonomer, c0, mw, true, lineNo)


proc parseFeedDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 4:
    fail(lineNo, "feed definition syntax: feed NAME SPECIES CONCENTRATION")
  let feedName = tokens[1]
  validateIdentifier(feedName, lineNo, "feed name")
  let sp = m.getSpeciesId(tokens[2], lineNo)
  let c = parseF(tokens[3], lineNo, "feed concentration")
  if c < 0.0:
    fail(lineNo, "feed concentration must be >= 0")
  var fid: int
  if m.feedByName.hasKey(feedName):
    fid = m.feedByName[feedName]
  else:
    fid = m.feeds.len
    m.feedByName[feedName] = fid
    m.feeds.add FeedDef(name: feedName, concentrations: newSeq[float](m.species.len))
  if m.feeds[fid].concentrations.len < m.species.len:
    m.feeds[fid].concentrations.setLen(m.species.len)
  if m.feeds[fid].concentrations[sp] != 0.0:
    fail(lineNo, "feed component already declared: " & feedName & " " & tokens[2])
  m.feeds[fid].concentrations[sp] = c

proc parseEndGroupDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 3:
    fail(lineNo, "endgroup syntax: endgroup NAME CONTRIBUTION")
  validateIdentifier(tokens[1], lineNo, "endgroup name")
  let egKey = modelNameKey(tokens[1])
  for existing in m.egByName.keys:
    if modelNameKey(existing) == egKey:
      let id = m.egByName[existing]
      if m.egMassSource[id] != "builtin":
        fail(lineNo, "endgroup name conflicts with already declared endgroup: " & existing)
      let contribution = parseF(tokens[2], lineNo, "endgroup contribution")
      m.setEgMass(existing, contribution, "explicit")
      return
  let contribution = parseF(tokens[2], lineNo, "endgroup contribution")
  m.setEgMass(tokens[1], contribution, "explicit")

proc parseSizeBytes(s0: string; lineNo: int): int64 =
  ## Parity port from slimmc-copo: parse a size literal like "50GB"/"512MB"
  ## into a raw byte count.
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
  result = int64(parseF(s.strip(), lineNo, "size") * mult)

proc parseMemoryLimitDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 2:
    fail(lineNo, "usage: memory_limit 50GB")
  m.memoryPolicy.hasLimit = true
  m.memoryPolicy.limitBytes = parseSizeBytes(tokens[1], lineNo)

proc parseAtMemoryDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  ## Items 52/53/39: unknown tokens must raise (previously silently
  ## discarded); "compact_dead"/"drop_dead_seq" removed entirely (they
  ## were permanently inert no-ops, never implemented -- slimmc's Chain
  ## is already a fixed-size record with nothing to compact); matching
  ## is case-sensitive (previously .toLowerAscii(), an unjustified
  ## departure from every other keyword in this parser).
  if tokens.len < 3 or tokens.len > 4:
    fail(lineNo, "usage: at_memory 50GB save|stop")
  m.memoryPolicy.hasLimit = true
  m.memoryPolicy.limitBytes = parseSizeBytes(tokens[1], lineNo)
  for a in tokens[2 .. ^1]:
    case a
    of "save": m.memoryPolicy.snapshotOnLimit = true
    of "stop": m.memoryPolicy.stopOnLimit = true
    else: fail(lineNo, "at_memory: unknown token: " & a & " (expected save and/or stop)")

proc parsePolymerDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len != 3:
    fail(lineNo, "polymer syntax: polymer NAME active|dead")

  var kind: PoolKind
  case tokens[2]
  of "active": kind = pkActive
  of "dead": kind = pkDead
  else: fail(lineNo, "polymer kind must be active or dead")

  discard m.addPool(tokens[1], kind, lineNo)

proc parseRateDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len notin [3, 4, 5]:
    fail(lineNo, "rate syntax: rate NAME value OR rate NAME const value OR rate NAME arr Apre Ea")

  let name = tokens[1]
  if tokens.len == 3:
    let kval = parseF(tokens[2], lineNo, "rate value")
    if kval < 0.0:
      fail(lineNo, "rate value must be >= 0")
    discard m.addRate(RateDef(name: name, kind: rkFixed, kConst: kval), lineNo)
    return

  if tokens.len == 4 and tokens[2] == "const":
    let kval = parseF(tokens[3], lineNo, "rate value")
    if kval < 0.0:
      fail(lineNo, "rate value must be >= 0")
    discard m.addRate(RateDef(name: name, kind: rkFixed, kConst: kval), lineNo)
    return

  if tokens.len == 5 and tokens[2] == "arr":
    let apre = parseF(tokens[3], lineNo, "Arrhenius pre-exponential factor")
    let ea = parseF(tokens[4], lineNo, "Arrhenius activation energy")
    if apre < 0.0:
      fail(lineNo, "Arrhenius pre-exponential factor must be >= 0")
    discard m.addRate(RateDef(name: name, kind: rkArr, Apre: apre, Ea: ea, declaredArrhenius: true), lineNo)
    return

  fail(lineNo, "rate syntax: rate NAME value OR rate NAME const value OR rate NAME arr Apre Ea")

proc validateActionArgs(action: ActionKind; args: seq[string]; lineNo: int) =
  ## Items 56/57/58: "print" bare (no message) is now a parser error --
  ## the progress-printout behavior it used to double as moved to the
  ## new, distinct "print_info" action. "print" now always requires
  ## exactly one message argument.
  case action
  of eaPrint:
    if args.len != 1:
      fail(lineNo, "print syntax: print \"message\" (bare print with no message is no longer valid -- use print_info for a progress printout)")
  of eaPrintInfo:
    if args.len != 0:
      fail(lineNo, "print_info takes no arguments")
  of eaFeed:
    if args.len notin [2, 3]:
      fail(lineNo, "feed action syntax: feed NAME VOLUME [L|l|mL|ml|ML]; without a unit VOLUME is in L")
  of eaStop:
    if args.len != 0:
      fail(lineNo, "stop takes no arguments")
  else:
    discard

proc parseScheduledActionDecl(m: var Model; tokens: seq[string]; lineNo: int) =
  if tokens.len < 3:
    fail(lineNo, "scheduled action syntax: every PERIOD ACTION ... OR from START step PERIOD ACTION ... OR at TIME ACTION ...")

  case tokens[0]
  of "every":
    let dt = parseF(tokens[1], lineNo, "scheduled action period")
    if dt <= 0.0:
      fail(lineNo, "every period must be > 0")
    let action = parseActionKind(tokens[2], lineNo)
    if action == eaStop: fail(lineNo, "stop is only valid as a when action")
    let args = sliceTokens(tokens, 3, tokens.len)
    validateActionArgs(action, args, lineNo)
    # Periodic schedules start at t=0. `every DT ACTION` is the compact
    # spelling of `from 0 step DT ACTION`.
    m.scheduledActions.add ScheduledAction(
      startTime: 0.0, nextTime: 0.0, period: dt, repeat: true, remaining: -1, active: true,
      action: action, args: args, lineNo: lineNo
    )
  of "from":
    let start = parseF(tokens[1], lineNo, "scheduled action start")
    if start < 0.0:
      fail(lineNo, "from start must be >= 0")
    if tokens.len >= 7 and tokens[2] == "repeat" and tokens[4] == "every":
      let count = parseI64(tokens[3], lineNo, "repeat count")
      if count <= 0: fail(lineNo, "repeat COUNT must be > 0")
      let dt = parseF(tokens[5], lineNo, "scheduled action period")
      if dt <= 0.0: fail(lineNo, "every PERIOD must be > 0")
      let action = parseActionKind(tokens[6], lineNo)
      if action == eaStop: fail(lineNo, "stop is only valid as a when action")
      let args = sliceTokens(tokens, 7, tokens.len)
      validateActionArgs(action, args, lineNo)
      m.scheduledActions.add ScheduledAction(startTime: start, nextTime: start, period: dt, repeat: true, remaining: count, active: true, action: action, args: args, lineNo: lineNo)
    else:
      if tokens.len < 5 or tokens[2] != "step":
        fail(lineNo, "from syntax: from START step PERIOD ACTION ... OR from START repeat COUNT every PERIOD ACTION ...")
      let dt = parseF(tokens[3], lineNo, "scheduled action period")
      if dt <= 0.0: fail(lineNo, "from step period must be > 0")
      let action = parseActionKind(tokens[4], lineNo)
      if action == eaStop: fail(lineNo, "stop is only valid as a when action")
      let args = sliceTokens(tokens, 5, tokens.len)
      validateActionArgs(action, args, lineNo)
      m.scheduledActions.add ScheduledAction(startTime: start, nextTime: start, period: dt, repeat: true, remaining: -1, active: true, action: action, args: args, lineNo: lineNo)
  of "at":
    let t = parseF(tokens[1], lineNo, "scheduled action time")
    if t < 0.0:
      fail(lineNo, "at time must be >= 0")
    let action = parseActionKind(tokens[2], lineNo)
    if action == eaStop: fail(lineNo, "stop is only valid as a when action")
    let args = sliceTokens(tokens, 3, tokens.len)
    validateActionArgs(action, args, lineNo)
    m.scheduledActions.add ScheduledAction(
      startTime: t, nextTime: t, period: 0.0, repeat: false, remaining: 1, active: true,
      action: action, args: args, lineNo: lineNo
    )
  else:
    fail(lineNo, "internal parser error in scheduled action declaration")

proc toSpeciesTerm(m: Model; rt: RawTerm; lineNo: int): Term =
  if m.poolByName.hasKey(rt.name):
    fail(lineNo, rt.name & " is a polymer pool; use macro, not rxn")
  result = Term(sp: m.getSpeciesId(rt.name, lineNo), stoich: rt.stoich)

proc parseRxnLine(m: var Model; raw: RawLine) =
  let tokens = tokenize(raw.text, raw.lineNo)
  if tokens.len < 5:
    fail(raw.lineNo, "rxn syntax: rxn lhs -> rhs rate [eff]")

  let arrow = indexOfToken(tokens, "->")
  if arrow < 0:
    fail(raw.lineNo, "missing -> in rxn")

  var eff = 1.0
  var kIdx = tokens.len - 1
  if isFloatToken(tokens[^1]):
    eff = parseF(tokens[^1], raw.lineNo, "reaction efficiency")
    kIdx = tokens.len - 2

  if eff < 0.0 or eff > 1.0:
    fail(raw.lineNo, "rxn efficiency must be in [0,1]")
  if kIdx <= arrow:
    fail(raw.lineNo, "missing rate in rxn")

  let lhsRaw = parseSide(sliceTokens(tokens, 1, arrow), raw.lineNo)
  let rhsRaw = parseSide(sliceTokens(tokens, arrow + 1, kIdx), raw.lineNo)
  if lhsRaw.len == 0:
    fail(raw.lineNo, "rxn lhs cannot be empty")

  var lhs: seq[Term] = @[]
  var rhs: seq[Term] = @[]
  for rt in lhsRaw:
    let tr = m.toSpeciesTerm(rt, raw.lineNo)
    if tr.sp == m.monomerId:
      warn(raw.lineNo, "rxn uses monomer; monomer balance may not hold")
    lhs.add tr
  for rt in rhsRaw:
    let tr = m.toSpeciesTerm(rt, raw.lineNo)
    if tr.sp == m.monomerId:
      warn(raw.lineNo, "rxn produces monomer; monomer balance may not hold")
    rhs.add tr

  var kind: ChannelKind
  if lhs.len == 1 and lhs[0].stoich == 1:
    kind = chRxnUni
  elif lhs.len == 1 and lhs[0].stoich == 2:
    kind = chRxnBiSame
  elif lhs.len == 2 and lhs[0].stoich == 1 and lhs[1].stoich == 1:
    if lhs[0].sp == lhs[1].sp:
      kind = chRxnBiSame
      lhs = @[Term(sp: lhs[0].sp, stoich: 2)]
    else:
      kind = chRxnBiDiff
  else:
    fail(raw.lineNo, "rxn supports only A -> ..., A+B -> ..., or 2A -> ...")

  let kId = m.getRateId(tokens[kIdx], raw.lineNo)
  m.channels.add KmcChannel(
    name: "channel_" & $(m.channels.len + 1), expr: raw.text.strip(),
    lineNo: raw.lineNo, kind: kind, kId: kId, eff: eff,
    lhs: lhs, rhs: rhs, sp1: -1, sp2: -1,
    pool1: -1, pool2: -1, poolOut: -1
  )

proc exactlyOneMonomerAndOneSpecies(
  m: Model; names: seq[string]; lineNo: int
): tuple[sp: int, mon: int] =
  if names.len != 2:
    fail(lineNo, "expected two lhs terms")

  var mon = -1
  var sp = -1
  for nm in names:
    if m.poolByName.hasKey(nm):
      fail(lineNo, "unexpected polymer pool in this macro pattern: " & nm)
    let id = m.getSpeciesId(nm, lineNo)
    if id == m.monomerId:
      if mon >= 0: fail(lineNo, "two monomers found")
      mon = id
    else:
      if sp >= 0: fail(lineNo, "two non-monomer species found")
      sp = id

  if mon < 0: fail(lineNo, "macro init requires monomer")
  if sp < 0: fail(lineNo, "macro init requires initiating species")
  result = (sp: sp, mon: mon)

proc parseMacroLine(m: var Model; raw: RawLine) =
  let tokens = tokenize(raw.text, raw.lineNo)
  if tokens.len < 6:
    fail(raw.lineNo, "macro syntax: macro kind lhs -> rhs rate")

  let macroKind = tokens[1]
  if macroKind in ["termc", "termd", "termx", "termC", "termD", "termX"]:
    fail(raw.lineNo, "removed macro kind " & macroKind & " is not supported in v2.8.5; use term_c, term_d, term_x, or transfer")

  let arrow = indexOfToken(tokens, "->")
  if arrow < 0:
    fail(raw.lineNo, "missing -> in macro")

  let kIdx = tokens.len - 1
  let kId = m.getRateId(tokens[kIdx], raw.lineNo)
  let lhs = expandTerms(parseSide(sliceTokens(tokens, 2, arrow), raw.lineNo))
  let rhs = expandTerms(parseSide(sliceTokens(tokens, arrow + 1, kIdx), raw.lineNo))

  var ch = KmcChannel(
    name: "channel_" & $(m.channels.len + 1), expr: raw.text.strip(),
    lineNo: raw.lineNo, kId: kId, eff: 1.0,
    lhs: @[], rhs: @[], sp1: -1, sp2: -1,
    pool1: -1, pool2: -1, poolOut: -1
  )

  case macroKind
  of "init":
    if rhs.len != 1:
      fail(raw.lineNo, "macro init rhs must be one active polymer pool")
    let x = exactlyOneMonomerAndOneSpecies(m, lhs, raw.lineNo)
    ch.kind = chMacroInit
    ch.sp1 = x.sp
    ch.sp2 = x.mon
    ch.poolOut = m.requirePoolKind(rhs[0], pkActive, raw.lineNo)
    discard m.ensureEg(endGroupName(m.species[x.sp].name))

  of "prop":
    if lhs.len != 2 or rhs.len != 1:
      fail(raw.lineNo, "macro prop syntax: macro prop P+M -> P k")
    var p = -1
    var mon = -1
    for nm in lhs:
      if m.poolByName.hasKey(nm):
        p = m.requirePoolKind(nm, pkActive, raw.lineNo)
      else:
        let sp = m.getSpeciesId(nm, raw.lineNo)
        if sp != m.monomerId:
          fail(raw.lineNo, "macro prop small reactant must be monomer")
        mon = sp
    if p < 0 or mon < 0:
      fail(raw.lineNo, "macro prop requires active polymer and monomer")
    let rhsP = m.requirePoolKind(rhs[0], pkActive, raw.lineNo)
    if rhsP != p:
      fail(raw.lineNo, "macro prop must keep the same active pool")
    ch.kind = chMacroProp
    ch.pool1 = p
    ch.sp1 = mon

  of "deprop":
    if lhs.len != 1 or rhs.len != 2:
      fail(raw.lineNo, "macro deprop syntax: macro deprop P -> P+M k")
    let p = m.requirePoolKind(lhs[0], pkActive, raw.lineNo)
    var rhsP = -1
    var mon = -1
    for nm in rhs:
      if m.poolByName.hasKey(nm):
        rhsP = m.requirePoolKind(nm, pkActive, raw.lineNo)
      else:
        let sp = m.getSpeciesId(nm, raw.lineNo)
        if sp != m.monomerId:
          fail(raw.lineNo, "macro deprop product species must be monomer")
        mon = sp
    if rhsP != p:
      fail(raw.lineNo, "macro deprop must return to the same active pool")
    if mon < 0:
      fail(raw.lineNo, "macro deprop requires monomer product")
    ch.kind = chMacroDeprop
    ch.pool1 = p
    ch.sp1 = mon

  of "term_c":
    if lhs.len != 2 or rhs.len != 1:
      fail(raw.lineNo, "macro term_c syntax: macro term_c P+Q -> D k")
    ch.kind = chMacroTermC
    ch.pool1 = m.requirePoolKind(lhs[0], pkActive, raw.lineNo)
    ch.pool2 = m.requirePoolKind(lhs[1], pkActive, raw.lineNo)
    ch.poolOut = m.requirePoolKind(rhs[0], pkDead, raw.lineNo)

  of "term_d":
    if lhs.len != 2 or rhs.len != 2:
      fail(raw.lineNo, "macro term_d syntax: macro term_d P+Q -> D+D k")
    let d1 = m.requirePoolKind(rhs[0], pkDead, raw.lineNo)
    let d2 = m.requirePoolKind(rhs[1], pkDead, raw.lineNo)
    if d1 != d2:
      fail(raw.lineNo, "macro term_d must write both chains to the same dead pool")
    ch.kind = chMacroTermD
    ch.pool1 = m.requirePoolKind(lhs[0], pkActive, raw.lineNo)
    ch.pool2 = m.requirePoolKind(lhs[1], pkActive, raw.lineNo)
    ch.poolOut = d1

  of "term_x":
    if lhs.len != 2 or rhs.len != 1:
      fail(raw.lineNo, "macro term_x syntax: macro term_x P+S -> D k")
    var p = -1
    var sp = -1
    for nm in lhs:
      if m.poolByName.hasKey(nm):
        p = m.requirePoolKind(nm, pkActive, raw.lineNo)
      else:
        let id = m.getSpeciesId(nm, raw.lineNo)
        if id == m.monomerId:
          fail(raw.lineNo, "macro term_x small reactant cannot be monomer")
        sp = id
    if p < 0 or sp < 0:
      fail(raw.lineNo, "macro term_x requires active polymer and small species")
    ch.kind = chMacroTermX
    ch.pool1 = p
    ch.sp1 = sp
    ch.poolOut = m.requirePoolKind(rhs[0], pkDead, raw.lineNo)
    discard m.ensureEg(endGroupName(m.species[sp].name))

  of "transfer":
    # "transfer_h" used to be accepted as an alias here (same
    # chMacroTransferH channel, different spelling only) -- removed per
    # item 55's action/macro-alias cleanup. The formed_by *output* label
    # "transfer_h" (fbTransferH in slimmc_types.nim) is unrelated and
    # unchanged -- it records which chemical channel actually formed a
    # given chain, not which spelling the user typed in the model file.
    let transferName = "transfer"
    if lhs.len != 2 or rhs.len != 2:
      fail(raw.lineNo, "macro " & transferName & " syntax: macro " & transferName & " P+X -> D+R k")
    var p = -1
    var donor = -1
    for nm in lhs:
      if m.poolByName.hasKey(nm):
        p = m.requirePoolKind(nm, pkActive, raw.lineNo)
      else:
        let id = m.getSpeciesId(nm, raw.lineNo)
        if id == m.monomerId:
          fail(raw.lineNo, "macro " & transferName & " donor cannot be monomer")
        donor = id
    var d = -1
    var radical = -1
    for nm in rhs:
      if m.poolByName.hasKey(nm):
        d = m.requirePoolKind(nm, pkDead, raw.lineNo)
      else:
        let id = m.getSpeciesId(nm, raw.lineNo)
        if id == m.monomerId:
          fail(raw.lineNo, "macro " & transferName & " radical product cannot be monomer")
        radical = id
    if p < 0 or donor < 0 or d < 0 or radical < 0:
      fail(raw.lineNo, "macro " & transferName & " requires active pool + donor -> dead pool + radical species")
    ch.kind = chMacroTransferH
    ch.pool1 = p
    ch.poolOut = d
    ch.sp1 = donor
    ch.sp2 = radical
    discard m.ensureEg(endGroupName(m.species[radical].name))

  of "transfer_m":
    if lhs.len != 2 or rhs.len != 2:
      fail(raw.lineNo, "macro transfer_m syntax: macro transfer_m P+M -> D+P k")

    var p = -1
    var mon = -1
    for nm in lhs:
      if m.poolByName.hasKey(nm):
        p = m.requirePoolKind(nm, pkActive, raw.lineNo)
      else:
        let id = m.getSpeciesId(nm, raw.lineNo)
        if id != m.monomerId:
          fail(raw.lineNo, "macro transfer_m small reactant must be the declared monomer")
        mon = id

    var d = -1
    var pNew = -1
    for nm in rhs:
      if not m.poolByName.hasKey(nm):
        fail(raw.lineNo, "macro transfer_m products must be one dead and one active polymer pool")
      let pid = m.getPoolId(nm, raw.lineNo)
      case m.pools[pid].kind
      of pkDead:
        if d >= 0:
          fail(raw.lineNo, "macro transfer_m rhs must contain exactly one dead polymer pool")
        d = pid
      of pkActive:
        if pNew >= 0:
          fail(raw.lineNo, "macro transfer_m rhs must contain exactly one active polymer pool")
        pNew = pid

    if p < 0 or mon < 0 or d < 0 or pNew < 0:
      fail(raw.lineNo, "macro transfer_m requires active pool + monomer -> dead pool + active pool")
    ch.kind = chMacroTransferM
    ch.pool1 = p
    ch.pool2 = pNew
    ch.poolOut = d
    ch.sp1 = mon
    discard m.ensureBuiltinEg(transferMonomerEndGroupName(m.species[mon].name), -1.008, "builtin_transfer_m")

  else:
    fail(raw.lineNo, "unknown macro kind: " & macroKind)

  m.channels.add ch

proc parseReactionLine(m: var Model; raw: RawLine) =
  let tokens = tokenize(raw.text, raw.lineNo)
  if tokens.len == 0:
    return
  case tokens[0]
  of "rxn": m.parseRxnLine(raw)
  of "macro": m.parseMacroLine(raw)
  else: fail(raw.lineNo, "internal parser error: expected rxn or macro")

proc parseConditionalActionLine(m: var Model; raw: RawLine) =
  let tokens = tokenize(raw.text, raw.lineNo)
  if tokens.len < 5:
    fail(raw.lineNo, "when syntax: when CONDITION [and CONDITION ...] ACTION ...")
  if tokens[0] != "when":
    fail(raw.lineNo, "internal parser error: expected when")

  var pos = 1
  var conditions: seq[AtomicCondition] = @[]
  while true:
    var observable: ConditionalObservableKind
    var speciesId = -1
    if pos >= tokens.len:
      fail(raw.lineNo, "missing condition in when")
    if tokens[pos] == "X":
      if m.monomerId < 0:
        fail(raw.lineNo, "when X requires a declared monomer")
      if pos + 1 >= tokens.len:
        fail(raw.lineNo, "when X requires the monomer name: when X MONOMER ...")
      let monomerName = tokens[pos + 1]
      if monomerName != m.species[m.monomerId].name:
        fail(raw.lineNo, "unknown monomer in when X: " & monomerName & "; expected " & m.species[m.monomerId].name)
      observable = woConversion
      pos += 2
    elif tokens[pos] == "c":
      if pos + 1 >= tokens.len:
        fail(raw.lineNo, "when c requires a species name")
      observable = woSpeciesConc
      speciesId = m.getSpeciesId(tokens[pos + 1], raw.lineNo)
      pos += 2
    else:
      fail(raw.lineNo, "when condition must start with X or c SPECIES")

    if pos + 1 >= tokens.len:
      fail(raw.lineNo, "incomplete when condition")
    var comparison: ComparisonKind
    case tokens[pos]
    of ">": comparison = coGreater
    of "<": comparison = coLess
    else: fail(raw.lineNo, "when supports only > or <")
    inc pos
    let threshold = parseF(tokens[pos], raw.lineNo, "when threshold")
    inc pos
    conditions.add AtomicCondition(observable: observable, speciesId: speciesId,
                                   comparison: comparison, threshold: threshold)
    if pos < tokens.len and tokens[pos] == "and":
      inc pos
      continue
    break

  if pos >= tokens.len:
    fail(raw.lineNo, "missing action in when")
  let action = parseActionKind(tokens[pos], raw.lineNo)
  inc pos
  let args = sliceTokens(tokens, pos, tokens.len)
  validateActionArgs(action, args, raw.lineNo)
  let first = conditions[0]
  m.conditionalActions.add ConditionalAction(
    active: true, conditions: conditions,
    observable: first.observable, speciesId: first.speciesId,
    comparison: first.comparison, threshold: first.threshold,
    action: action, args: args, lineNo: raw.lineNo
  )

proc parseModel*(path: string): Model =
  var m = initModel()
  m.modelSourceFile = path
  let lines = readFile(path).splitLines()

  for i, line in lines:
    let lineNo = i + 1
    let tokens = tokenize(line, lineNo)
    if tokens.len == 0:
      continue

    case tokens[0]
    of "desc": m.parseDesc(tokens, lineNo, line)
    of "var": m.parseVar(tokens, lineNo)
    of "param": m.parseParam(tokens, lineNo, line)
    of "species": m.parseSpeciesDecl(tokens, lineNo)
    of "feed": m.parseFeedDecl(tokens, lineNo)
    of "monomer": m.parseMonomerDecl(tokens, lineNo)
    of "endgroup": m.parseEndGroupDecl(tokens, lineNo)
    of "polymer": m.parsePolymerDecl(tokens, lineNo)
    of "rate": m.parseRateDecl(tokens, lineNo)
    of "rxn", "macro": m.rawReactions.add RawLine(text: line, lineNo: lineNo)
    of "every", "from", "at": m.parseScheduledActionDecl(tokens, lineNo)
    of "when": m.rawConditionalActions.add RawLine(text: line, lineNo: lineNo)
    of "memory_limit": m.parseMemoryLimitDecl(tokens, lineNo)
    of "at_memory": m.parseAtMemoryDecl(tokens, lineNo)
    else: fail(lineNo, "unknown model statement: " & tokens[0])

  if m.V <= 0.0:
    fail(0, "kmc_volume must be > 0")
  if m.T <= 0.0:
    fail(0, "T must be > 0")
  if m.tEnd < 0.0:
    fail(0, "t_end must be >= 0")
  if m.maxEvents <= 0:
    fail(0, "max_steps must be > 0")
  if m.whenCheckEvents <= 0:
    fail(0, "when_check_events must be > 0")
  for raw in m.rawReactions:
    m.parseReactionLine(raw)

  for raw in m.rawConditionalActions:
    m.parseConditionalActionLine(raw)

  if m.feeds.len > 0:
    if not m.hasInitVolume or m.initVolumeMl <= 0.0:
      fail(0, "param init_volume must be > 0 when feed is used")
  elif m.hasInitVolume and m.initVolumeMl <= 0.0:
    fail(0, "init_volume must be > 0")
  m.currentVolumeMl = m.initVolumeMl
  m.resolveSeriesVars()

  let parts = splitFile(path)
  m.runId = validateModelRunId(path)

  if m.outputDir.len == 0:
    let base = if parts.dir.len == 0 or parts.dir == ".": "results" else: parts.dir / "results"
    m.outputDir = base / m.runId

  result = m
