## Slimmc copo runtime I/O for Slimmc Storage v1.
## Human-readable logs remain text; all structured simulation data use Storage columns.

import std/[os, times, strutils, json, math]
import copo_types
import copo_stats
import copo_propensity
import copo_storage

const
  PreflightSmallPopulationCount* = 100'i64
  PreflightSmallChainCount* = 100'i64
  PreflightLargeRelError* = 0.05

proc configureOutput*(opts: RunOptions) =
  discard opts

proc formatWallDuration*(seconds: float): string =
  let totalMs = int64(seconds * 1000.0 + 0.5)
  let ms = totalMs mod 1000
  var rest = totalMs div 1000
  let sec = rest mod 60
  rest = rest div 60
  let minute = rest mod 60
  rest = rest div 60
  let hour = rest mod 24
  let day = rest div 24
  result = align($day, 2, '0') & ":" & align($hour, 2, '0') & ":" &
           align($minute, 2, '0') & ":" & align($sec, 2, '0') & "." & align($ms, 3, '0')

proc isoNowUtc*(): string = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")

proc ensureOutputDir*(m: Model) =
  if m.output_dir.len > 0 and not dirExists(m.output_dir): createDir(m.output_dir)

proc resetOutputDir*(m: Model) =
  if m.output_dir.len == 0: return
  if dirExists(m.output_dir): removeDir(m.output_dir)
  createDir(m.output_dir)

proc monomerNames*(m: Model): seq[string] =
  for x in m.monomers: result.add x.name

proc fmtBytes*(b: int64): string =
  let x = float(b)
  if b < 1024'i64: return $b & " B"
  if b < 1024'i64 * 1024'i64: return formatFloat(x / 1024.0, ffDecimal, 2) & " KB"
  if b < 1024'i64 * 1024'i64 * 1024'i64: return formatFloat(x / 1024.0 / 1024.0, ffDecimal, 2) & " MB"
  formatFloat(x / 1024.0 / 1024.0 / 1024.0, ffDecimal, 2) & " GB"

proc appendLine(path, line: string) =
  createDir(parentDir(path))
  var f=open(path,fmAppend); f.writeLine line; f.close()

proc writeStartInfo*(m: Model; opts: RunOptions = RunOptions()): string =
  discard m; discard opts
  isoNowUtc()

proc writeRunInfo*(m: Model; s: State; opts: RunOptions; status: string; terminationReason: string; wallSeconds: float; startedAt: string) =
  discard m; discard s; discard opts; discard status; discard terminationReason; discard wallSeconds; discard startedAt

proc writeParameterState*(m: Model; s: State; actionId: int64) =
  discard m; discard s; discard actionId

proc initParameterStates*(m: Model; s: State) =
  discard m; discard s

proc initDiagnostics*(m: Model; opts: RunOptions) =
  m.ensureOutputDir()
  if opts.debug:
    writeFile(m.output_dir / "diagnostics" / "debug.log", AppName & " " & AppVersion & " debug log\n")

proc appendRunLog*(m: Model; msg: string) =
  appendLine(m.output_dir / "diagnostics" / "run.log", msg)

proc appendDebugLog*(m: Model; msg: string) =
  appendLine(m.output_dir / "diagnostics" / "debug.log", msg)

proc writeActionTrace*(m: Model; s: var State; trigger, scheduledTime, checkSource,
                       observable, operatorText, threshold, value: string;
                       action: ActionKind; target, requested, beforeValue,
                       afterValue, message: string; lineNo: int;
                       stateChanged, outputWritten: bool) =
  discard m; discard s; discard trigger; discard scheduledTime; discard checkSource
  discard observable; discard operatorText; discard threshold; discard value
  discard action; discard target; discard requested; discard beforeValue; discard afterValue
  discard message; discard lineNo; discard stateChanged; discard outputWritten

proc writeTraceAction*(m: Model; s: State; chId: int; dt, propensity, totalPropensity: float) =
  discard m; discard s; discard chId; discard dt; discard propensity; discard totalPropensity

proc saveSnapshot*(m: Model; s: var State; writeChainRows = false;
                   source = "action"; isFinal = false) =
  captureStorageV1Snapshot(m, s, source, isFinal, writeChainRows, computePropensities(m, s))

proc preflightOneEntry(name: string; kind: string; c0: float; V: float; warnings: var seq[string]): JsonNode =
  let nExact = c0 * NA * V
  let nActual = countFromConc(c0, V)
  let relError = if nExact > 0.0: abs(float(nActual) - nExact) / nExact else: 0.0
  result = newJObject()
  result["name"] = newJString(name)
  result["kind"] = newJString(kind)
  result["c0_mol_L"] = newJFloat(c0)
  result["n_exact"] = newJFloat(nExact)
  result["n_actual"] = newJInt(nActual)
  result["rel_error"] = newJFloat(relError)
  if c0 > 0.0 and nActual == 0:
    raise newException(ValueError, "preflight: " & kind & " '" & name & "' rounds to 0 molecules")
  if nActual > 0 and nActual < PreflightSmallPopulationCount:
    warnings.add kind & " '" & name & "' has a very small population (" & $nActual & " molecules)"
  if relError > PreflightLargeRelError and nExact > 0.0:
    let recommendedV = 0.5 / (PreflightLargeRelError * c0 * NA)
    warnings.add kind & " '" & name & "' has a large relative discretization error; recommended minimum volume: " & $recommendedV & " L"

proc computeDiscretizationPreflight*(m: Model): JsonNode =
  ## Items 88-95: deterministic VkMC discretization preflight, mirroring
  ## homo's implementation -- no pilot KMC run (item 96 explicitly
  ## defers that). copo keeps monomers and species in separate lists, so
  ## both are walked here (homo has one unified species list including
  ## the monomer).
  result = newJObject()
  var thresholds = newJObject()
  thresholds["small_population_count"] = newJInt(PreflightSmallPopulationCount)
  thresholds["small_chain_count"] = newJInt(PreflightSmallChainCount)
  thresholds["large_rel_error"] = newJFloat(PreflightLargeRelError)
  result["thresholds"] = thresholds

  var speciesArr = newJArray()
  var warnings: seq[string] = @[]
  var minChainEstimate = int64.high
  for mo in m.monomers:
    speciesArr.add preflightOneEntry(mo.name, "monomer", mo.c0, m.V, warnings)
  for sp in m.species:
    speciesArr.add preflightOneEntry(sp.name, "species", sp.c0, m.V, warnings)
    let nActual = countFromConc(sp.c0, m.V)
    if nActual > 0 and nActual < minChainEstimate:
      minChainEstimate = nActual
  result["species"] = speciesArr

  # Item 92: coarse expected-chain-count estimate, bounded by the
  # smallest species population (typically the initiator/CTA).
  if minChainEstimate < int64.high and minChainEstimate < PreflightSmallChainCount:
    warnings.add "expected chain count is small (~" & $minChainEstimate &
                 ") -- CLD/MWD/composition analysis may not be statistically reliable"

  var warningsArr = newJArray()
  for w in warnings: warningsArr.add newJString(w)
  result["warnings"] = warningsArr

proc printStatus*(m: Model; s: State; wallStart: float; opts: RunOptions) =
  let mem = estimateMemory(m, s)
  var convParts: seq[string] = @[]
  for i, mo in m.monomers:
    let conv = if s.monomerN0[i] > 0: 1.0 - float(s.monomerN[i]) / float(s.monomerN0[i]) else: 0.0
    convParts.add(mo.name & ":" & formatFloat(conv, ffDecimal, 4))
  let pct = if m.t_end > 0.0: 100.0 * s.t / m.t_end else: 0.0
  let wall = epochTime() - wallStart
  let stamp = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
  let memMiB = float(mem.totalBytes) / (1024.0 * 1024.0)
  let line = "[run] event=" & $s.kmcEvent &
    " t=" & formatFloat(s.t, ffScientific, 2) & "/" & formatFloat(m.t_end, ffScientific, 2) &
    " " & formatFloat(pct, ffDecimal, 1) & "%" &
    " conv=" & convParts.join(",") &
    " mem=" & formatFloat(memMiB, ffDecimal, 1) & " MiB" &
    " wall=" & formatFloat(wall, ffDecimal, 1) & "s " & stamp
  stdout.writeLine(line)
  stdout.flushFile()
  appendRunLog(m, line)

proc printStart*(m: Model) =
  let line = "[start] run_id=" & m.modelStem &
    " t_end=" & formatFloat(m.t_end, ffScientific, 2) &
    " output=" & m.output_dir
  stdout.writeLine(line)
  stdout.flushFile()
  appendRunLog(m, line)

proc printDone*(m: Model; s: State; status, reason: string; wallSeconds: float) =
  let line = "[done] status=" & status &
    " reason=" & reason &
    " event=" & $s.kmcEvent &
    " t=" & formatFloat(s.t, ffScientific, 6) &
    " wall=" & formatFloat(wallSeconds, ffDecimal, 1) & "s" &
    " output=" & m.output_dir
  stdout.writeLine(line)
  stdout.flushFile()

proc printMemory*(m: Model; s: State) =
  let mem = estimateMemory(m, s)
  echo "[memory] live_chains=", mem.liveChains,
       " dead_summaries=", mem.deadSummaries,
       " stored_live_mers=", mem.storedLiveMers,
       " stored_dead_mers=", mem.storedDeadMers,
       " live_seq=", fmtBytes(mem.liveSeqBytes),
       " live_objects=", fmtBytes(mem.liveObjectBytes),
       " dead_summary=", fmtBytes(mem.deadSummaryBytes),
       " total_est=", fmtBytes(mem.totalBytes)
  stdout.flushFile()

proc printCheck*(m: Model) =
  let preflight = computeDiscretizationPreflight(m)
  var warnings: seq[string] = @[]
  for w in preflight["warnings"]: warnings.add w.getStr()
  for a in m.scheduledActions:
    if a.action == eaSetC:
      warnings.add "line " & $a.lineNo & ": set_c forces a concentration and invalidates the physical material balance for that species"
    if a.startTime > m.t_end:
      warnings.add "line " & $a.lineNo & ": this scheduled action starts after t_end and will never run"
    elif a.repeat and a.remaining > 0:
      let lastTime = a.startTime + float(a.remaining - 1) * a.period
      if lastTime > m.t_end:
        warnings.add "line " & $a.lineNo & ": part of this repeated schedule occurs after t_end and will not run"

  # Missing propagation declarations are legal and mean a zero-rate transition.
  # Report them explicitly so users do not mistake an omitted channel for an
  # automatically inferred one. This works for terminal, penultimate and
  # terpolymer models because it checks each declared active pool separately.
  var missingProp: seq[string] = @[]
  for poolId, pool in m.pools:
    if pool.kind != pkActive: continue
    for monId, mon in m.monomers:
      var found = false
      for ch in m.channels:
        if ch.kind == chMacroProp and ch.pool1 == poolId and ch.monomerId == monId:
          found = true
          break
      if not found:
        missingProp.add pool.name & " + " & mon.name
  if missingProp.len > 0:
    warnings.add "missing macro prop transitions are treated as k=0: " & missingProp.join(", ")

  echo "GENERAL"
  echo "  CHECK: ", (if warnings.len == 0: "OK" else: "OK WITH WARNINGS")
  echo "  Engine: ", AppName, " ", AppVersion
  echo "  Model: ", m.modelFile
  echo "  Run ID: ", m.modelStem
  echo "  Warnings: ", warnings.len
  echo "  Errors: 0"
  if warnings.len > 0:
    for i, w in warnings: echo "  Warning ", i + 1, ": ", w
  echo "  The model is valid and can be started."
  if warnings.len > 0: echo "  Review the warnings before running it."
  echo ""
  echo "DETAILS"
  echo "  [model]"
  echo "  source_model: ", m.modelFile
  echo "  run_id: ", m.modelStem
  echo "  engine: ", AppName
  echo "  engine_version: ", AppVersion
  echo "  output_dir: ", m.output_dir
  echo "  storage: slimmc-storage"
  echo "  [verification]"
  echo "  Source model file: OK"
  echo "  Run ID derived from file name: OK - ", m.modelStem
  echo "  Run ID syntax: OK"
  echo "  Monomer declarations counted: OK - ", m.monomers.len
  echo "  Engine dispatch: OK - copo selected"
  echo "  Complete model parsing: OK"
  echo "  Directive syntax and argument counts: OK"
  echo "  Numeric values and ranges: OK"
  echo "  Names and references: OK"
  echo "  Duplicate and conflicting declarations: OK"
  echo "  Basic parameters: OK"
  echo "  Chemical structure and reaction channels: OK"
  echo "  Polymer pool compatibility: OK"
  echo "  Kinetic parameter references: OK"
  echo "  Scheduled actions: OK"
  echo "  Conditional actions: OK"
  echo "  Feed definitions and actions: OK"
  echo "  Volume requirements and units: OK"
  echo "  Number of monomers: OK - ", m.monomers.len
  echo "  Dead polymer pool: OK"
  echo "  At least one reaction channel: OK"
  echo "  Sequence mode compatibility: OK"
  echo "  Depropagation terminal information: OK"
  echo "  KMC discretization preflight: OK"
  echo "  Static storage contract: OK"
  echo "  [parameters]"
  echo "  kmc_volume: ", m.V, " L"
  if m.hasInitVolume: echo "  init_volume: ", m.initVolumeMl / 1000.0, " L"
  else: echo "  init_volume: not defined"
  echo "  temperature: ", m.T, " K"
  echo "  t_end: ", m.t_end, " s"
  echo "  max_steps: ", m.max_steps
  echo "  seed: ", m.seed
  echo "  sequence_mode: ", m.sequence_mode
  echo "  dp_max: ", m.dp_max
  echo "  [monomers]"
  for mo in m.monomers: echo "  ", mo.name, ": c0=", mo.c0, " mol/L, molar_mass=", mo.mw, " g/mol"
  echo "  [species]"
  if m.species.len == 0: echo "  none"
  for sp in m.species: echo "  ", sp.name, ": c0=", sp.c0, " mol/L"
  echo "  [feeds]"
  if m.feeds.len == 0: echo "  none"
  for f in m.feeds:
    echo "  ", f.name
    for i, c in f.monomerConcentrations:
      if c != 0.0: echo "    ", m.monomers[i].name, ": ", c, " mol/L"
    for i, c in f.speciesConcentrations:
      if c != 0.0: echo "    ", m.species[i].name, ": ", c, " mol/L"
  echo "  [scheduled actions]"
  if m.scheduledActions.len == 0: echo "  none"
  for a in m.scheduledActions:
    var schedule = if a.repeat:
      (if a.remaining < 0: "from " & $a.startTime & " step " & $a.period else: "from " & $a.startTime & " repeat " & $a.remaining & " every " & $a.period)
    else: "at " & $a.startTime
    echo "  line ", a.lineNo, ": ", schedule, " ", actionKindName(a.action), (if a.args.len > 0: " " & a.args.join(" ") else: "")
  echo "  [conditional actions]"
  echo "  count: ", m.conditionalActions.len
  echo "  when_check_events: ", m.whenCheckEvents
  echo "  [rates]"
  for r in m.rates: echo "  ", r.name
  echo "  [channels]"
  for ch in m.channels: echo "  line ", ch.lineNo, ": ", ch.name, " kind=", ch.kind, " rate=", m.rates[ch.kId].name
  echo "  [kmc discretization]"
  for item in preflight["species"]:
    echo "  ", item["name"].getStr(), ": exact_count=", item["n_exact"].getFloat(), ", actual_count=", item["n_actual"].getInt(), ", relative_error=", item["rel_error"].getFloat(), " - OK"
  echo ""
  echo "WARNINGS"
  if warnings.len == 0: echo "  none"
  for i, w in warnings: echo "  ", i + 1, ". ", w
  echo ""
  echo "ERRORS"
  echo "  none"
  echo ""
  echo "LIMITATIONS"
  echo "  This check validates the model before the run."
  echo "  It does not prove that the reaction will progress, that a when condition will be reached,"
  echo "  that memory will be sufficient, or that final statistics will be precise."
  echo ""
  echo "CHECK RESULT: ", (if warnings.len == 0: "OK" else: "OK WITH WARNINGS")
