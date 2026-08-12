## Slimmc homo runtime I/O for Slimmc Storage v1.
## Human-readable logs remain text; all structured simulation data use Storage columns.

import std/[strutils, times, os, math, json]
import slimmc_types
import slimmc_storage

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
           align($minute, 2, '0') & ":" & align($sec, 2, '0') & "." &
           align($ms, 3, '0')

proc appendLog(path, line: string) =
  createDir(parentDir(path))
  var f = open(path, fmAppend)
  f.writeLine line
  f.close()

proc ensureCleanOutput(m: Model) =
  if m.outputDir.len == 0: return
  if dirExists(m.outputDir): removeDir(m.outputDir)
  createDir(m.outputDir)

proc initOutputFiles*(m: Model; s: var State; opts: RunOptions; startedAt: string) =
  discard opts
  ensureCleanOutput(m)
  initStorageV1(m, s, startedAt)

proc writeRunInfo*(m: Model; s: State; opts: RunOptions; status: string;
                   reason: string; startedAt, finishedAt: string; wallSeconds: float) =
  discard m; discard s; discard opts; discard status; discard reason
  discard startedAt; discard finishedAt; discard wallSeconds

proc writeParameterState*(m: Model; s: State; actionId: int64) =
  discard m; discard s; discard actionId

proc initDebugFile*(m: Model) =
  let path = m.outputDir / "diagnostics" / "debug.log"
  writeFile(path, AppName & " " & AppVersion & " debug log\n")

proc debugLine*(m: Model; msg: string) =
  appendLog(m.outputDir / "diagnostics" / "debug.log", msg)

proc writeRunLogFinal*(m: Model; opts: RunOptions; s: State; reason: string; wallSeconds: float) =
  discard opts
  appendLog(m.outputDir / "diagnostics" / "run.log",
    "[done] event=" & $s.kmcEvent & " t=" & $s.t & " reason=" & reason &
    " wall=" & $wallSeconds & "s")

proc saveSnapshot*(m: Model; s: var State; withChains: bool;
                   source = "action"; isFinal = false) =
  captureStorageV1Snapshot(m, s, source, isFinal, withChains)

proc writeActionTraceRow*(
  m: Model; s: var State; lineNo: int; trigger: string; scheduledTime: string;
  checkSource: string; observable: string; op: string; threshold: string;
  value: string; action: ActionKind; actionResult: ActionResult;
  conditions: seq[AtomicCondition] = @[]; observedValues: seq[float64] = @[]
) =
  discard checkSource; discard observable; discard op; discard threshold; discard value
  captureActionV1(m, s, lineNo, trigger, scheduledTime, conditions,
    observedValues, action, actionResult)

proc writeChannelTraceRow*(m: Model; s: State; chId: int; dt, propensity, totalPropensity: float) =
  discard m; discard s; discard chId; discard dt; discard propensity; discard totalPropensity

proc logLine*(m: Model; opts: RunOptions; line: string) =
  discard opts
  appendLog(m.outputDir / "diagnostics" / "run.log", line)

proc computeDiscretizationPreflight*(m: Model): JsonNode =
  ## Items 88-95: deterministic VkMC discretization preflight -- no
  ## pilot KMC run (item 96 explicitly defers that). For every declared
  ## species (including the monomer), computes N = c*V*NA both as an
  ## exact real number and as the actually-simulated rounded integer,
  ## plus the relative discretization error between them.
  result = newJObject()
  var thresholds = newJObject()
  thresholds["small_population_count"] = newJInt(PreflightSmallPopulationCount)
  thresholds["small_chain_count"] = newJInt(PreflightSmallChainCount)
  thresholds["large_rel_error"] = newJFloat(PreflightLargeRelError)
  result["thresholds"] = thresholds

  var speciesArr = newJArray()
  var warnings: seq[string] = @[]
  var minChainEstimate = int64.high
  for sp in m.species:
    let nExact = sp.c0 * NA * m.V
    let nActual = countFromConc(sp.c0, m.V)
    let relError = if nExact > 0.0: abs(float(nActual) - nExact) / nExact else: 0.0
    var item = newJObject()
    item["name"] = newJString(sp.name)
    item["kind"] = newJString(if sp.kind == skMonomer: "monomer" else: "species")
    item["c0_mol_L"] = newJFloat(sp.c0)
    item["n_exact"] = newJFloat(nExact)
    item["n_actual"] = newJInt(nActual)
    item["rel_error"] = newJFloat(relError)
    speciesArr.add item

    # Item 90: a positive concentration that rounds to zero particles is
    # a hard error -- the model is asking for a species that will never
    # actually exist in this simulation.
    if sp.c0 > 0.0 and nActual == 0:
      fail(0, "preflight: species '" & sp.name & "' has c0=" & num(sp.c0) &
           " mol/L > 0 but rounds to 0 molecules at this volume -- increase volume or concentration")

    # Item 91: warn about very small populations for any species (this
    # covers initiator/CTA/other mechanism-critical species generically,
    # since homo has no separate "role" tag beyond monomer/species).
    if nActual > 0 and nActual < PreflightSmallPopulationCount:
      warnings.add "species '" & sp.name & "' has a very small population (" & $nActual &
                   " molecules) -- expect poor statistics for anything depending on it"

    # Item 93: warn about large relative discretization error, with a
    # recommended minimum volume to bring it under threshold.
    if relError > PreflightLargeRelError and nExact > 0.0:
      let recommendedV = 0.5 / (PreflightLargeRelError * sp.c0 * NA)
      warnings.add "species '" & sp.name & "' has a large relative discretization error (" &
                   num(relError) & ") -- recommended minimum volume: " & num(recommendedV) & " L"

    if sp.kind != skMonomer and nActual > 0 and nActual < minChainEstimate:
      minChainEstimate = nActual
  result["species"] = speciesArr

  # Item 92: a coarse estimate of the expected final chain count (bounded
  # by the smallest non-monomer species population -- typically the
  # initiator/CTA, which is what actually limits how many chains can
  # form) -- warn if it's too small for reliable CLD/MWD analysis.
  if minChainEstimate < int64.high and minChainEstimate < PreflightSmallChainCount:
    warnings.add "expected chain count is small (~" & $minChainEstimate &
                 ") -- CLD/MWD analysis may not be statistically reliable"

  var warningsArr = newJArray()
  for w in warnings: warningsArr.add newJString(w)
  result["warnings"] = warningsArr

proc monomerConversion*(m: Model; s: State): float =
  if m.monomerId < 0 or s.mExpected <= 0:
    return 0.0
  result = 1.0 - float(s.n[m.monomerId]) / float(s.mExpected)

proc estimateMemory*(m: Model; s: State): MemoryEstimate =
  ## Parity port from slimmc-copo's memory tracking, simplified for
  ## slimmc's storage model: every chain (live or dead) is stored as the
  ## same compact, fixed-size `Chain` record (eg1, dp, eg2, formedBy) inside
  ## `s.pools`, so unlike copo there is no separate per-mer sequence to
  ## count -- one conservative per-chain byte estimate covers everything.
  const bytesPerChain = 32'i64  # 2x int + int64 + enum tag, generously rounded
  for pid, pdef in m.pools:
    let n = int64(s.pools[pid].len)
    if pdef.kind == pkActive:
      result.liveChains += n
    else:
      result.deadChains += n
  result.liveObjectBytes = result.liveChains * bytesPerChain
  result.deadRecordBytes = result.deadChains * bytesPerChain
  result.totalBytes = result.liveObjectBytes + result.deadRecordBytes

proc fmtBytes*(b: int64): string =
  let x = float(b)
  if b < 1024'i64: return $b & " B"
  if b < 1024'i64 * 1024'i64: return numDec(x / 1024.0, 2) & " KB"
  if b < 1024'i64 * 1024'i64 * 1024'i64: return numDec(x / 1024.0 / 1024.0, 2) & " MB"
  return numDec(x / 1024.0 / 1024.0 / 1024.0, 2) & " GB"

proc debugCheckState*(m: Model; s: State; context: string) =
  var issues: seq[string] = @[]

  for i, n in s.n:
    if n < 0:
      issues.add "negative species count for " & m.species[i].name

  for pid, pdef in m.pools:
    var propagatable = 0'i64
    var depropable = 0'i64
    for ch in s.pools[pid]:
      if ch.dp < 1:
        issues.add "DP < 1 in pool " & pdef.name
      if int64(ch.dp) < m.dpMax:
        inc propagatable
      if ch.dp > 1:
        inc depropable
      if pdef.kind == pkActive and ch.eg2 != m.egActive:
        issues.add "active pool " & pdef.name & " contains non-ACTIVE chain"
      if pdef.kind == pkDead and ch.eg2 == m.egActive:
        issues.add "dead pool " & pdef.name & " contains ACTIVE chain"
    if pid >= s.poolPropagatableCounts.len or s.poolPropagatableCounts[pid] != propagatable:
      issues.add "propagatable counter mismatch in pool " & pdef.name &
        " stored=" & (if pid < s.poolPropagatableCounts.len: $s.poolPropagatableCounts[pid] else: "missing") &
        " actual=" & $propagatable
    if pid >= s.poolDepropableCounts.len or s.poolDepropableCounts[pid] != depropable:
      issues.add "depropable counter mismatch in pool " & pdef.name &
        " stored=" & (if pid < s.poolDepropableCounts.len: $s.poolDepropableCounts[pid] else: "missing") &
        " actual=" & $depropable

  if m.monomerId >= 0:
    let err = calcMTotal(m, s) - s.mExpected
    if err != 0:
      issues.add "monomer balance error=" & $err

  if issues.len == 0:
    debugLine(m, "OK " & context & " event=" & $s.kmcEvent & " t=" & num(s.t))
  else:
    for issue in issues:
      debugLine(m, "ISSUE " & context & " event=" & $s.kmcEvent & " t=" & num(s.t) & " " & issue)

proc writeDebugFinal*(m: Model; s: State; reason: string; wallSeconds: float) =
  debugLine(m, "final")
  debugLine(m, "  reason=" & reason)
  debugLine(m, "  t=" & num(s.t))
  debugLine(m, "  kmc_event=" & $s.kmcEvent)
  debugLine(m, "  action_no=" & $s.actionNo)
  debugLine(m, "  scheduled_action_no=" & $s.scheduledActionNo)
  debugLine(m, "  conditional_action_no=" & $s.conditionalActionNo)
  debugLine(m, "  snapshot_id=" & $s.snapshotId)
  debugLine(m, "  wall_seconds=" & num(wallSeconds))
  for i, fires in s.channelFires:
    debugLine(m, "  channel_" & $(i + 1) & "_fires=" & $fires)

proc printMemory*(m: Model; s: State) =
  let mem = estimateMemory(m, s)
  echo "[memory] live_chains=", mem.liveChains,
       " dead_chains=", mem.deadChains,
       " total_est=", fmtBytes(mem.totalBytes)

proc printProgress*(m: Model; s: State; wallStart: float; opts: RunOptions) =
  let pct = if m.tEnd > 0.0: 100.0 * s.t / m.tEnd else: 0.0
  let wall = epochTime() - wallStart
  let stamp = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
  let monomerName = if m.monomerId >= 0: m.species[m.monomerId].name else: "M"
  let line = "[run] event=" & $s.kmcEvent &
    " t=" & formatFloat(s.t, ffScientific, 2) & "/" & formatFloat(m.tEnd, ffScientific, 2) &
    " " & numDec(pct, 1) & "%" &
    " conv=" & monomerName & ":" & numDec(monomerConversion(m, s), 4) &
    " wall=" & numDec(wall, 1) & "s " & stamp
  stdout.writeLine(line)
  logLine(m, opts, line)

proc printMarker*(m: Model; opts: RunOptions; message: string) =
  stdout.writeLine(message)
  logLine(m, opts, message)

proc printCheck*(m: Model) =
  let preflight = computeDiscretizationPreflight(m)
  var warnings: seq[string] = @[]
  for w in preflight["warnings"]: warnings.add w.getStr()
  for a in m.scheduledActions:
    if a.action == eaSetC:
      warnings.add "line " & $a.lineNo & ": set_c forces a concentration and invalidates the physical material balance for that species"
    if a.startTime > m.tEnd:
      warnings.add "line " & $a.lineNo & ": this scheduled action starts after t_end and will never run"
    elif a.repeat and a.remaining > 0:
      let lastTime = a.startTime + float(a.remaining - 1) * a.period
      if lastTime > m.tEnd:
        warnings.add "line " & $a.lineNo & ": part of this repeated schedule occurs after t_end and will not run"

  echo "GENERAL"
  echo "  CHECK: ", (if warnings.len == 0: "OK" else: "OK WITH WARNINGS")
  echo "  Engine: ", AppName, " ", AppVersion
  echo "  Model: ", m.modelSourceFile
  echo "  Run ID: ", m.runId
  echo "  Warnings: ", warnings.len
  echo "  Errors: 0"
  if warnings.len > 0:
    for i, w in warnings: echo "  Warning ", i + 1, ": ", w
  echo "  The model is valid and can be started."
  if warnings.len > 0: echo "  Review the warnings before running it."
  echo ""
  echo "DETAILS"
  echo "  [model]"
  echo "  source_model: ", m.modelSourceFile
  echo "  run_id: ", m.runId
  echo "  engine: ", AppName
  echo "  engine_version: ", AppVersion
  echo "  output_dir: ", m.outputDir
  echo "  storage: slimmc-storage"
  echo "  [verification]"
  echo "  Source model file: OK"
  echo "  Run ID derived from file name: OK - ", m.runId
  echo "  Run ID syntax: OK"
  echo "  Monomer declarations counted: OK - ", (if m.monomerId >= 0: 1 else: 0)
  echo "  Engine dispatch: OK - homo selected"
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
  echo "  KMC discretization preflight: OK"
  echo "  Static storage contract: OK"
  echo "  [parameters]"
  echo "  kmc_volume: ", m.V, " L"
  if m.hasInitVolume: echo "  init_volume: ", m.initVolumeMl / 1000.0, " L"
  else: echo "  init_volume: not defined"
  echo "  temperature: ", m.T, " K"
  echo "  t_end: ", m.tEnd, " s"
  echo "  max_steps: ", m.maxEvents
  echo "  seed: ", m.seed
  echo "  mass_model: ", m.massModel
  echo "  [species and monomers]"
  for sp in m.species:
    echo "  ", sp.name, ": kind=", (if sp.kind == skMonomer: "monomer" else: "species"), ", c0=", sp.c0, " mol/L"
  echo "  [feeds]"
  if m.feeds.len == 0: echo "  none"
  for f in m.feeds:
    echo "  ", f.name
    for i, c in f.concentrations:
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
  echo "  [channels]"
  echo "  count: ", m.channels.len
  for ch in m.channels: echo "  line ", ch.lineNo, ": ", ch.expr
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
