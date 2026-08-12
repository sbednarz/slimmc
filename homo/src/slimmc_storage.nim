## Slimmc Storage v1 stage C for the homo engine: run container,
## snapshots and dense state columns. Structured output is written exclusively through Slimmc Storage.

import std/[os, json, times, strutils, tables, algorithm]
import slimmc_types
import slimmc_storage_validator
import storage/homo_schema
import ../../common/[npy_writer, result_json, results_types, sha256_file, storage_manifest, build_provenance]

proc utcNow(): string = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")

proc reasonId(reason: string): uint32 =
  case reason
  of "initial": 0'u32
  of "scheduled": 1'u32
  of "action", "save", "save_chains": 2'u32
  of "final": 4'u32
  else: 3'u32

proc jsonArgs(args: seq[string]): JsonNode =
  result = newJArray()
  for arg in args: result.add %arg

proc homoConditionNode(m: Model; c: AtomicCondition): JsonNode =
  result = newJObject()
  result["observable"] = %conditionalObservableName(m, c)
  result["comparison"] = %comparisonName(c.comparison)
  result["threshold"] = %c.threshold

proc homoResolvedModelNode(m: Model): JsonNode =
  result = newJObject()
  result["schema"] = %"slimmc-resolved-model-v1"
  result["kinetic_model"] = %"homo"
  result["desc"] = %(if m.hasDescription: m.description else: "")

  result["parameters"] = %*{
    "kmc_volume_L": (if m.hasInitVolume and m.currentVolumeMl > 0.0: m.V * m.initVolumeMl / m.currentVolumeMl else: m.V),
    "init_volume_mL": (if m.hasInitVolume: m.initVolumeMl else: 0.0),
    "initial_temperature_K": m.T,
    "t_end_s": m.tEnd,
    "max_events": m.maxEvents,
    "when_check_events": m.whenCheckEvents,
    "seed": m.seed,
    "dp_max": m.dpMax,
    "sequence_mode": m.sequenceMode,
    "mass_model": m.massModel
  }

  result["memory_policy"] = %*{
    "has_limit": m.memoryPolicy.hasLimit,
    "limit_bytes": m.memoryPolicy.limitBytes,
    "save_on_limit": m.memoryPolicy.snapshotOnLimit,
    "stop_on_limit": m.memoryPolicy.stopOnLimit
  }

  var monomers = newJArray()
  var species = newJArray()
  for item in m.species:
    var node = newJObject()
    node["name"] = %item.name
    node["initial_concentration_mol_L"] = %item.c0
    if item.hasMw: node["molar_mass_g_mol"] = %item.mw
    else: node["molar_mass_g_mol"] = newJNull()
    if item.kind == skMonomer: monomers.add node else: species.add node
  result["monomers"] = monomers
  result["species"] = species

  var feeds = newJArray()
  for feed in m.feeds:
    var components = newJArray()
    for i, c in feed.concentrations:
      if c > 0.0:
        components.add %*{"name": m.species[i].name, "concentration_mol_L": c}
    feeds.add %*{"name": feed.name, "components": components}
  result["feeds"] = feeds

  var polymers = newJArray()
  for pool in m.pools:
    polymers.add %*{"name": pool.name, "population_activity": (if pool.kind == pkActive: "live" else: "dead")}
  result["polymers"] = polymers

  var endgroups = newJArray()
  for i, name in m.egNames:
    var node = newJObject()
    node["name"] = %name
    node["mass_known"] = %m.egMassKnown[i]
    if m.egMassKnown[i]: node["molar_mass_g_mol"] = %m.egMass[i]
    else: node["molar_mass_g_mol"] = newJNull()
    node["source"] = %m.egMassSource[i]
    node["builtin"] = %(m.egMassSource[i] == "builtin")
    endgroups.add node
  result["endgroups"] = endgroups

  var rates = newJArray()
  for rate in m.rates:
    var node = newJObject()
    node["name"] = %rate.name
    node["kind"] = %(if rate.kind == rkFixed: "constant" else: "arrhenius")
    if rate.kind == rkFixed:
      node["value"] = %rate.kConst
    else:
      node["A"] = %rate.Apre
      node["Ea_J_mol"] = %rate.Ea
    rates.add node
  result["rates"] = rates

  var reactions = newJArray()
  for channel in m.channels:
    var node = newJObject()
    node["name"] = %channel.name
    node["kind"] = %channelKindName(channel.kind)
    node["source_line"] = %channel.lineNo
    node["rate"] = %(if channel.kId >= 0: m.rates[channel.kId].name else: "")
    node["efficiency"] = %channel.eff
    reactions.add node
  result["reactions"] = reactions

  var orderedActions: seq[tuple[lineNo: int, node: JsonNode]] = @[]
  for action in m.scheduledActions:
    var trigger = newJObject()
    if action.repeat:
      trigger["kind"] = %"periodic"
      trigger["start_s"] = %action.startTime
      trigger["step_s"] = %action.period
      if action.remaining > 0: trigger["count"] = %action.remaining
    else:
      trigger["kind"] = %"at"
      trigger["time_s"] = %action.startTime
    let node = %*{
      "source_order": action.lineNo,
      "trigger": trigger,
      "action": {"kind": actionKindName(action.action), "args": jsonArgs(action.args)}
    }
    orderedActions.add (action.lineNo, node)
  for action in m.conditionalActions:
    var conditions = newJArray()
    for condition in action.conditions: conditions.add homoConditionNode(m, condition)
    let trigger = %*{"kind": "when", "operator": "and", "conditions": conditions}
    let node = %*{
      "source_order": action.lineNo,
      "trigger": trigger,
      "action": {"kind": actionKindName(action.action), "args": jsonArgs(action.args)}
    }
    orderedActions.add (action.lineNo, node)
  orderedActions.sort(proc(a, b: tuple[lineNo: int, node: JsonNode]): int = cmp(a.lineNo, b.lineNo))
  var actions = newJArray()
  for item in orderedActions: actions.add item.node
  result["actions"] = actions

  var variables = newJArray()
  for variable in m.variables:
    variables.add %*{"kind": variable.kind, "name": variable.name, "value": variable.value, "unit": variable.unit}
  result["variables"] = variables

  result["engine_specific"] = %*{
    "chain_dp_dtype": "int",
    "chain_dp_max": $high(int),
    "builtin_endgroups": ["H", "U", "ACTIVE"]
  }

proc metadataNode(m: Model; s: State; status: RunStatus; startedAt: string;
                  finishedAt = ""; wallTime = 0.0; exitCode = 0; terminationReason = "";
                  validationStatus = "not_run"; validationWarnings = 0;
                  validationErrors = 0; storageManifest: StorageManifest = StorageManifest()): JsonNode =
  result = newJObject()
  result["run_id"] = %m.runId
  result["input_model_file"] = %"input.model"
  result["source_model_name"] = %extractFilename(m.modelSourceFile)
  let inputPath = m.outputDir / "input.model"
  let schemaPath = m.outputDir / "schema.jsonl"
  if fileExists(inputPath): result["input_model_sha256"] = %sha256File(inputPath)
  result["storage"] = %StorageName
  result["storage_format_version"] = %StorageFormatVersion
  if fileExists(schemaPath): result["schema_sha256"] = %sha256File(schemaPath)
  result["seed"] = %($uint64(m.seed))
  result["engine"] = %AppName
  result["cli_version"] = %AppVersion
  result["engine_version"] = %AppVersion
  result["started_at_utc"] = %startedAt
  if finishedAt.len > 0: result["finished_at_utc"] = %finishedAt
  if finishedAt.len > 0: result["wall_time_s"] = %wallTime
  result["run_status"] = %($status)
  result["termination_reason"] = %terminationReason
  result["exit_code"] = %exitCode
  let initialKmcVolumeL = if m.hasInitVolume and m.currentVolumeMl > 0.0: m.V * m.initVolumeMl / m.currentVolumeMl else: m.V
  result["kmc_volume_L"] = %initialKmcVolumeL
  result["initial_kmc_volume_L"] = %initialKmcVolumeL
  result["current_kmc_volume_L"] = %m.V
  if m.hasInitVolume:
    result["initial_volume_mL"] = %m.initVolumeMl
    result["current_volume_mL"] = %m.currentVolumeMl
  result["volume_mode"] = %(if m.feeds.len > 0: "variable" else: "constant")
  result["avogadro_constant_mol_inv"] = %AvogadroConstantMolInv
  result["initial_temperature_K"] = %m.T
  result["kinetic_model"] = %"homo"
  result["n_monomers"] = %(if m.monomerId >= 0: 1 else: 0)
  result["dp_max"] = %($m.dpMax)
  result["sequence_mode"] = %m.sequenceMode
  result["engine_chain_dp_dtype"] = %"int"
  result["engine_chain_dp_max"] = %($high(int))
  result["platform"] = %(hostOS & "-" & hostCPU)
  result["threads"] = %1
  result["compiler"] = %("Nim " & NimVersion)
  result["build_mode"] = %(when defined(release): "release" else: "debug")
  let provenance = executableProvenance()
  for key, value in provenance: result[key] = value
  var executionInfo = %*{
    "engine": AppName,
    "engine_version": AppVersion,
    "cli_version": AppVersion,
    "started_at_utc": startedAt,
    "status": $status,
    "exit_code": exitCode,
    "platform": hostOS & "-" & hostCPU,
    "threads": 1,
    "compiler": "Nim " & NimVersion,
    "build_mode": (when defined(release): "release" else: "debug")
  }
  if finishedAt.len > 0:
    executionInfo["finished_at_utc"] = %finishedAt
    executionInfo["wall_time_s"] = %wallTime
  for key, value in provenance: executionInfo[key] = value
  result["execution"] = executionInfo
  result["validation_status"] = %validationStatus
  result["validation_warning_count"] = %validationWarnings
  result["validation_error_count"] = %validationErrors
  result["channel_trace_enabled"] = %(s.storageTraceKmcEvents.len > 0 or s.channelTraceRowsWritten > 0)
  result["channel_trace_complete"] = %(not s.channelTraceTruncated)
  result["channel_trace_rows"] = %s.channelTraceRowsWritten
  result["channel_trace_truncated"] = %s.channelTraceTruncated
  var variables = newJArray()
  for v in m.variables:
    variables.add(%*{"kind": v.kind, "name": v.name, "value": v.value, "unit": v.unit})
  result["variables"] = variables
  result["model"] = homoResolvedModelNode(m)
  var storageInfo = %*{
    "name": StorageName,
    "format_version": StorageFormatVersion,
    "complete": status == rsCompleted and validationErrors == 0,
    "hash_algorithm": StorageHashAlgorithm,
    "hash_schema": StorageHashSchema
  }
  if storageManifest.hash.len > 0:
    storageInfo["hash"] = %storageManifest.hash
    storageInfo["manifest_file"] = %StorageChecksumFile
    storageInfo["file_count"] = %storageManifest.fileCount
    storageInfo["total_bytes"] = %storageManifest.totalBytes
  else:
    storageInfo["hash"] = newJNull()
  result["storage_info"] = storageInfo


proc captureKineticParameterSetV1*(m: Model; s: var State; hasSourceAction: bool; sourceActionId: uint64 = 0'u64) =
  let setId = uint64(s.storageV1KineticSetIds.len)
  s.storageV1KineticSetIds.add setId
  s.storageV1KineticStartEvents.add uint64(s.kmcEvent)
  s.storageV1KineticStartTimes.add s.t
  s.storageV1KineticHasSourceAction.add hasSourceAction
  s.storageV1KineticSourceActionIds.add(if hasSourceAction: sourceActionId else: 0'u64)
  var pid = 0'u32
  template addValue(v: float64) =
    doAssert v == v and v != Inf and v != -Inf
    s.storageV1KineticValueSetIds.add setId
    s.storageV1KineticParameterIds.add pid
    s.storageV1KineticValues.add v
    inc pid
  addValue(m.T)
  for rid, rate in m.rates:
    addValue(m.rateValue(rid))
    if rate.declaredArrhenius:
      addValue(rate.Apre)
      addValue(rate.Ea)

proc parseOptionalFloat(text: string): float64 =
  if text.len == 0 or text == "nan": return NaN
  try: result = parseFloat(text)
  except ValueError: result = NaN

proc actionTypeId(action: ActionKind): uint32 = uint32(ord(action))

proc triggerTypeId(trigger: string): uint32 =
  case trigger
  of "at": 1'u32
  of "every": 2'u32
  of "when": 3'u32
  else: 0'u32

proc conditionObservableId(m: Model; a: AtomicCondition): uint32 =
  case a.observable
  of woConversion: 1'u32
  of woSpeciesConc: uint32(2 + a.speciesId)

proc conditionOperatorId(a: AtomicCondition): uint32 =
  case a.comparison
  of coGreater: 1'u32
  of coLess: 2'u32

proc targetId(m: Model; target: string): uint32 =
  if target.len == 0: return 0'u32
  if target == "T": return 1'u32
  for i, rate in m.rates:
    if rate.name == target: return uint32(2 + i)
  let base = 2 + m.rates.len
  for i, sp in m.species:
    if sp.name == target: return uint32(base + i)
  result = 0'u32

proc captureActionV1*(m: Model; s: var State; lineNo: int; trigger,
                      scheduledTime: string; conditions: seq[AtomicCondition];
                      observedValues: seq[float64]; action: ActionKind;
                      actionResult: ActionResult) =
  let actionId = uint64(max(0'i64, s.actionNo - 1))
  doAssert actionId == uint64(s.storageV1ActionIds.len)
  let hasSnapshot = action in {eaSave, eaSaveChains, eaSetK, eaAddK, eaSetTemp, eaAddTemp}
  let hasKinetic = action in {eaSetK, eaAddK, eaSetTemp, eaAddTemp}
  s.storageV1ActionIds.add actionId
  s.storageV1ActionKmcEvents.add uint64(s.kmcEvent)
  s.storageV1ActionTimes.add s.t
  s.storageV1ActionSourceLines.add uint32(lineNo)
  s.storageV1ActionTypeIds.add actionTypeId(action)
  s.storageV1ActionTriggerTypeIds.add triggerTypeId(trigger)
  s.storageV1ActionScheduledTimes.add parseOptionalFloat(scheduledTime)
  s.storageV1ActionTargetIds.add targetId(m, actionResult.target)
  s.storageV1ActionRequestedValues.add parseOptionalFloat(actionResult.requested)
  s.storageV1ActionBeforeValues.add(if actionResult.hasNumeric: actionResult.before else: NaN)
  s.storageV1ActionAfterValues.add(if actionResult.hasNumeric: actionResult.after else: NaN)
  s.storageV1ActionStateChanged.add actionResult.stateChanged
  s.storageV1ActionOutputWritten.add(actionResult.outputWritten.len > 0 or
    action in {eaSave, eaSaveChains, eaSetK, eaAddK, eaSetTemp, eaAddTemp})
  s.storageV1ActionHasSnapshot.add hasSnapshot
  s.storageV1ActionSnapshotIds.add(if hasSnapshot and s.storageV1SnapshotIds.len > 0:
    s.storageV1SnapshotIds[^1] else: 0'u64)
  s.storageV1ActionHasKineticSet.add hasKinetic
  s.storageV1ActionKineticSetIds.add(if hasKinetic and s.storageV1KineticSetIds.len > 0:
    s.storageV1KineticSetIds[^1] else: 0'u64)
  s.storageV1ActionMessages.add actionResult.message
  if action == eaFeed:
    let fid = m.feedByName[actionResult.target]
    let doseMl = parseOptionalFloat(actionResult.requested)
    let beforeMl = actionResult.before
    let afterMl = actionResult.after
    let kmcAfter = m.V
    let kmcBefore = if afterMl > 0.0: kmcAfter * beforeMl / afterMl else: NaN
    s.storageV1FeedActionIds.add actionId
    s.storageV1FeedIds.add uint32(fid)
    s.storageV1FeedDoseMl.add doseMl
    s.storageV1FeedVolumeBeforeMl.add beforeMl
    s.storageV1FeedVolumeAfterMl.add afterMl
    s.storageV1FeedKmcVolumeBeforeL.add kmcBefore
    s.storageV1FeedKmcVolumeAfterL.add kmcAfter
  doAssert conditions.len == observedValues.len
  for i, condition in conditions:
    let conditionRecordId = uint64(s.storageV1ConditionRecordIds.len)
    s.storageV1ConditionRecordIds.add conditionRecordId
    s.storageV1ConditionActionIds.add actionId
    s.storageV1ConditionIndexes.add uint32(i)
    s.storageV1ConditionObservableIds.add conditionObservableId(m, condition)
    s.storageV1ConditionOperatorIds.add conditionOperatorId(condition)
    s.storageV1ConditionThresholds.add condition.threshold
    s.storageV1ConditionObservedValues.add observedValues[i]
    s.storageV1ConditionMet.add true

proc initStorageV1*(m: Model; s: var State; startedAt: string) =
  let root = m.outputDir
  for d in ["snapshots", "state", "chains", "sequences", "moments",
            "channel_events", "actions", "action_conditions", "feed_events", "monomer_balance", "species_balance", "kinetic_parameters", "diagnostics", ".work"]:
    createDir(root / d)
  createDir(root / "kinetic_parameters" / "sets")
  createDir(root / "kinetic_parameters" / "values")
  createDir(root / "diagnostics" / "channel_trace")
  if fileExists(root / "RESULTS_COMPLETE"):
    removeFile(root / "RESULTS_COMPLETE")
  copyFile(m.modelSourceFile, root / "input.model")
  writeJsonlAtomic(root / "schema.jsonl", schemaRecords(m))
  writePrettyJsonAtomic(root / "run_metadata.json", metadataNode(m, s, rsRunning, startedAt))
  writeTextAtomic(root / "diagnostics" / "run.log", "Slimmc run started at " & startedAt & "\n")
  captureKineticParameterSetV1(m, s, false)

proc captureChainsV1(m: Model; s: var State; sid: uint64) =
  type ChainV1Key = tuple[populationId: uint32, poolId: uint32, originId: uint32,
                          dp: uint64, leftEndId: uint32, rightEndId: uint32]
  var grouped = initTable[ChainV1Key, uint64]()
  for pid, pool in s.pools:
    let populationId = if m.pools[pid].kind == pkActive: 0'u32 else: 1'u32
    for ch in pool:
      doAssert ch.dp >= 1
      let key: ChainV1Key = (populationId, uint32(pid + 1), uint32(ord(ch.formedBy)),
        uint64(ch.dp), uint32(ch.eg1 + 2), uint32(ch.eg2 + 2))
      grouped[key] = grouped.getOrDefault(key, 0'u64) + 1'u64
  var keys: seq[ChainV1Key] = @[]
  for key in grouped.keys: keys.add key
  keys.sort(proc(a, b: ChainV1Key): int =
    result = cmp(a.populationId, b.populationId)
    if result == 0: result = cmp(a.poolId, b.poolId)
    if result == 0: result = cmp(a.dp, b.dp)
    if result == 0: result = cmp(a.leftEndId, b.leftEndId)
    if result == 0: result = cmp(a.rightEndId, b.rightEndId)
    if result == 0: result = cmp(a.originId, b.originId))
  # Pure small-molecule `rxn` models legitimately have no monomer and no
  # chain records. Avoid indexing m.species[-1] when an ordinary `save`
  # snapshot captures an empty chain population.
  if keys.len == 0:
    return
  doAssert m.monomerId >= 0 and m.monomerId < m.species.len
  let monomerMw = m.species[m.monomerId].mw
  for key in keys:
    let count = grouped[key]
    let eg1 = int(key.leftEndId) - 2
    let eg2 = int(key.rightEndId) - 2
    var mass = float64(key.dp) * monomerMw
    if m.massModel == "with_end_groups":
      mass += m.egMass[eg1] + m.egMass[eg2]
    doAssert mass > 0.0 and mass == mass and mass != Inf and mass != -Inf
    let mol = float64(count) / AvogadroConstantMolInv
    s.storageV1ChainRecordIds.add uint64(s.storageV1ChainRecordIds.len)
    s.storageV1ChainSnapshotIds.add sid
    s.storageV1ChainPopulationIds.add key.populationId
    s.storageV1ChainPoolIds.add key.poolId
    s.storageV1ChainOriginIds.add key.originId
    s.storageV1ChainDp.add key.dp
    s.storageV1ChainMolarMass.add mass
    s.storageV1ChainCounts.add count
    s.storageV1ChainMoles.add mol
    s.storageV1ChainConcentrations.add mol / m.V
    s.storageV1ChainLeftEndIds.add key.leftEndId
    s.storageV1ChainRightEndIds.add key.rightEndId
    s.storageV1ChainSequenceOffsets.add 0'u64
    s.storageV1ChainSequenceLengths.add 0'u64

proc captureMomentsV1(m: Model; s: var State; sid: uint64; rowStart: int) =
  let monomerMw = if m.monomerId >= 0: m.species[m.monomerId].mw else: 0.0
  for scopeId in 0'u32 .. 2'u32:
    for basisId in 0'u32 .. 1'u32:
      var chainCount = 0'u64
      var sumDp = 0.0
      var sumDp2 = 0.0
      var sumMass = 0.0
      var sumMass2 = 0.0
      var sumMass3 = 0.0
      for i in rowStart ..< s.storageV1ChainRecordIds.len:
        let populationId = s.storageV1ChainPopulationIds[i]
        if scopeId == 1'u32 and populationId != 0'u32: continue
        if scopeId == 2'u32 and populationId != 1'u32: continue
        let count = s.storageV1ChainCounts[i]
        let countF = float64(count)
        let dp = float64(s.storageV1ChainDp[i])
        let mass = if basisId == 0'u32: dp * monomerMw
                   else: s.storageV1ChainMolarMass[i]
        doAssert mass > 0.0 and mass == mass and mass != Inf and mass != -Inf
        chainCount += count
        sumDp += countF * dp
        sumDp2 += countF * dp * dp
        sumMass += countF * mass
        sumMass2 += countF * mass * mass
        sumMass3 += countF * mass * mass * mass
      let empty = chainCount == 0'u64
      let dpN = if empty: NaN else: sumDp / float64(chainCount)
      let dpW = if sumDp == 0.0: NaN else: sumDp2 / sumDp
      let mn = if empty: NaN else: sumMass / float64(chainCount)
      let mw = if sumMass == 0.0: NaN else: sumMass2 / sumMass
      let mz = if sumMass2 == 0.0: NaN else: sumMass3 / sumMass2
      let dispersity = if mn != mn or mn == 0.0: NaN else: mw / mn
      s.storageV1MomentSnapshotIds.add sid
      s.storageV1MomentPopulationScopeIds.add scopeId
      s.storageV1MomentMassBasisIds.add basisId
      s.storageV1MomentChainCounts.add chainCount
      s.storageV1MomentSumDp.add sumDp
      s.storageV1MomentSumDp2.add sumDp2
      s.storageV1MomentDpN.add dpN
      s.storageV1MomentDpW.add dpW
      s.storageV1MomentSumMass.add sumMass
      s.storageV1MomentSumMass2.add sumMass2
      s.storageV1MomentSumMass3.add sumMass3
      s.storageV1MomentMn.add mn
      s.storageV1MomentMw.add mw
      s.storageV1MomentMz.add mz
      s.storageV1MomentDispersity.add dispersity

proc captureStorageV1Snapshot*(m: Model; s: var State; reason: string;
                               isFinal: bool; hasChains: bool) =
  if s.storageV1SnapshotIds.len > 0 and
     s.storageV1KmcEvents[^1] == uint64(s.kmcEvent) and
     timeClose(s.storageV1Times[^1], s.t) and
     s.storageV1StateRevisions[^1] == s.stateRevision:
    let needsChainCapture = hasChains and not s.storageV1HasChains[^1]
    s.storageV1HasChains[^1] = s.storageV1HasChains[^1] or hasChains
    if needsChainCapture:
      let chainRowStart = s.storageV1ChainRecordIds.len
      captureChainsV1(m, s, s.storageV1SnapshotIds[^1])
      captureMomentsV1(m, s, s.storageV1SnapshotIds[^1], chainRowStart)
    if isFinal:
      s.storageV1IsFinal[^1] = true
      s.storageV1ReasonIds[^1] = reasonId(reason)
    return
  let sid = uint64(s.storageV1SnapshotIds.len)
  s.storageV1SnapshotIds.add sid
  s.storageV1StateRevisions.add s.stateRevision
  s.storageV1Times.add s.t
  s.storageV1KmcEvents.add uint64(s.kmcEvent)
  s.storageV1ReasonIds.add reasonId(reason)
  s.storageV1IsFinal.add isFinal
  s.storageV1HasChains.add hasChains
  s.storageV1HasSequences.add false
  s.storageV1ParameterSetIds.add uint64(max(0'i64, s.parameterStateId - 1))
  s.storageV1VolumeMl.add(if m.hasInitVolume: m.currentVolumeMl else: NaN)
  s.storageV1KmcVolumeL.add m.V
  var liveCount = 0'u64
  var deadCount = 0'u64
  for pid, chains in s.pools:
    if m.pools[pid].kind == pkActive:
      liveCount += uint64(chains.len)
    else:
      deadCount += uint64(chains.len)
  s.storageV1ChainCountLive.add liveCount
  s.storageV1ChainCountDead.add deadCount
  s.storageV1ChainCountTotal.add liveCount + deadCount
  if m.hasInitVolume:
    let physicalScale = (m.currentVolumeMl / 1000.0) / m.V
    for sp in 0 ..< m.species.len:
      let initial = m.species[sp].c0 * (m.initVolumeMl / 1000.0)
      let dosed = float64(s.speciesDosedN[sp]) / AvogadroConstantMolInv * physicalScale
      let total = float64(s.speciesExternalN[sp]) / AvogadroConstantMolInv * physicalScale
      let free = float64(s.n[sp]) / AvogadroConstantMolInv * physicalScale
      s.storageV1SpeciesBalanceSnapshotIds.add sid
      s.storageV1SpeciesBalanceEntityIds.add uint32(sp)
      s.storageV1SpeciesInitialMoles.add initial
      s.storageV1SpeciesDosedMoles.add dosed
      s.storageV1SpeciesTotalMoles.add total
      s.storageV1SpeciesFreeMoles.add free
      s.storageV1SpeciesConsumedMoles.add total - free
    if m.monomerId >= 0:
      let introduced = float64(s.mExpected) / AvogadroConstantMolInv * physicalScale
      let free = float64(s.n[m.monomerId]) / AvogadroConstantMolInv * physicalScale
      let incorporated = introduced - free
      s.storageV1MonomerBalanceSnapshotIds.add sid
      s.storageV1MonomerBalanceMonomerIds.add 0'u32
      s.storageV1MonomerInitialMoles.add m.species[m.monomerId].c0 * (m.initVolumeMl / 1000.0)
      s.storageV1MonomerIntroducedMoles.add introduced
      s.storageV1MonomerFreeMoles.add free
      s.storageV1MonomerIncorporatedMoles.add incorporated
      s.storageV1MonomerConversion.add(if introduced > 0.0: incorporated / introduced else: 0.0)
  var eid = 0'u32
  template addState(c: uint64) =
    s.storageV1StateSnapshotIds.add sid
    s.storageV1StateEntityIds.add eid
    s.storageV1StateCounts.add c
    let mol = float64(c) / AvogadroConstantMolInv
    s.storageV1StateMoles.add mol
    s.storageV1StateConcentrations.add mol / m.V
    inc eid
  for count in s.n: addState(uint64(count))
  var live = 0'u64
  var dead = 0'u64
  for pid, pool in s.pools:
    addState(uint64(pool.len))
    if m.pools[pid].kind == pkActive: live += uint64(pool.len)
    else: dead += uint64(pool.len)
  addState(live)
  addState(dead)
  if hasChains:
    let chainRowStart = s.storageV1ChainRecordIds.len
    captureChainsV1(m, s, sid)
    captureMomentsV1(m, s, sid, chainRowStart)

  var totalEvents = 0'u64
  for cid in 0 ..< m.channels.len:
    let events = uint64(s.channelFires[cid])
    let productive = if m.channels[cid].kind in {chRxnUni, chRxnBiDiff, chRxnBiSame}:
      uint64(s.channelSuccesses[cid])
    else:
      events
    let nonproductive = if m.channels[cid].kind in {chRxnUni, chRxnBiDiff, chRxnBiSame}:
      uint64(s.channelFailures[cid])
    else:
      0'u64
    doAssert events == productive + nonproductive
    s.storageV1ChannelSnapshotIds.add sid
    s.storageV1ChannelIds.add uint32(cid)
    s.storageV1ChannelEventCounts.add events
    s.storageV1ChannelProductiveCounts.add productive
    s.storageV1ChannelNonproductiveCounts.add nonproductive
    totalEvents += events
  doAssert totalEvents == uint64(s.kmcEvent)

proc publishStorageV1*(m: Model; s: State; status: RunStatus; startedAt,
                       finishedAt: string; wallTime: float64; exitCode = 0; terminationReason = "") =
  let root = m.outputDir
  writeNpyUint64(root / "snapshots" / "snapshot_id.npy", s.storageV1SnapshotIds)
  writeNpyFloat64(root / "snapshots" / "time.npy", s.storageV1Times)
  writeNpyUint64(root / "snapshots" / "kmc_event.npy", s.storageV1KmcEvents)
  writeNpyUint32(root / "snapshots" / "snapshot_reason_id.npy", s.storageV1ReasonIds)
  writeNpyBool(root / "snapshots" / "is_final.npy", s.storageV1IsFinal)
  writeNpyBool(root / "snapshots" / "has_chains.npy", s.storageV1HasChains)
  writeNpyBool(root / "snapshots" / "has_sequences.npy", s.storageV1HasSequences)
  writeNpyUint64(root / "snapshots" / "kinetic_parameter_set_id.npy", s.storageV1ParameterSetIds)
  writeNpyFloat64(root / "snapshots" / "volume_mL.npy", s.storageV1VolumeMl)
  writeNpyFloat64(root / "snapshots" / "kmc_volume_L.npy", s.storageV1KmcVolumeL)
  writeNpyUint64(root / "snapshots" / "chain_count_live.npy", s.storageV1ChainCountLive)
  writeNpyUint64(root / "snapshots" / "chain_count_dead.npy", s.storageV1ChainCountDead)
  writeNpyUint64(root / "snapshots" / "chain_count_total.npy", s.storageV1ChainCountTotal)
  writeNpyUint64(root / "state" / "snapshot_id.npy", s.storageV1StateSnapshotIds)
  writeNpyUint32(root / "state" / "entity_id.npy", s.storageV1StateEntityIds)
  writeNpyUint64(root / "state" / "count.npy", s.storageV1StateCounts)
  writeNpyFloat64(root / "state" / "moles.npy", s.storageV1StateMoles)
  writeNpyFloat64(root / "state" / "concentration.npy", s.storageV1StateConcentrations)
  writeNpyUint64(root / "chains" / "chain_record_id.npy", s.storageV1ChainRecordIds)
  writeNpyUint64(root / "chains" / "snapshot_id.npy", s.storageV1ChainSnapshotIds)
  writeNpyUint32(root / "chains" / "population_id.npy", s.storageV1ChainPopulationIds)
  writeNpyUint32(root / "chains" / "pool_id.npy", s.storageV1ChainPoolIds)
  writeNpyUint32(root / "chains" / "origin_id.npy", s.storageV1ChainOriginIds)
  writeNpyUint64(root / "chains" / "dp.npy", s.storageV1ChainDp)
  writeNpyFloat64(root / "chains" / "molar_mass.npy", s.storageV1ChainMolarMass)
  writeNpyUint64(root / "chains" / "count.npy", s.storageV1ChainCounts)
  writeNpyFloat64(root / "chains" / "moles.npy", s.storageV1ChainMoles)
  writeNpyFloat64(root / "chains" / "concentration.npy", s.storageV1ChainConcentrations)
  writeNpyUint32(root / "chains" / "left_end_id.npy", s.storageV1ChainLeftEndIds)
  writeNpyUint32(root / "chains" / "right_end_id.npy", s.storageV1ChainRightEndIds)
  writeNpyUint64(root / "chains" / "sequence_offset.npy", s.storageV1ChainSequenceOffsets)
  writeNpyUint64(root / "chains" / "sequence_length.npy", s.storageV1ChainSequenceLengths)
  writeNpyUint32(root / "sequences" / "symbols.npy", @[])
  writeNpyUint64(root / "moments" / "snapshot_id.npy", s.storageV1MomentSnapshotIds)
  writeNpyUint32(root / "moments" / "population_scope_id.npy", s.storageV1MomentPopulationScopeIds)
  writeNpyUint32(root / "moments" / "mass_basis_id.npy", s.storageV1MomentMassBasisIds)
  writeNpyUint64(root / "moments" / "chain_count.npy", s.storageV1MomentChainCounts)
  writeNpyFloat64(root / "moments" / "sum_dp.npy", s.storageV1MomentSumDp)
  writeNpyFloat64(root / "moments" / "sum_dp2.npy", s.storageV1MomentSumDp2)
  writeNpyFloat64(root / "moments" / "dp_n.npy", s.storageV1MomentDpN)
  writeNpyFloat64(root / "moments" / "dp_w.npy", s.storageV1MomentDpW)
  writeNpyFloat64(root / "moments" / "sum_molar_mass.npy", s.storageV1MomentSumMass)
  writeNpyFloat64(root / "moments" / "sum_molar_mass2.npy", s.storageV1MomentSumMass2)
  writeNpyFloat64(root / "moments" / "sum_molar_mass3.npy", s.storageV1MomentSumMass3)
  writeNpyFloat64(root / "moments" / "mn.npy", s.storageV1MomentMn)
  writeNpyFloat64(root / "moments" / "mw.npy", s.storageV1MomentMw)
  writeNpyFloat64(root / "moments" / "mz.npy", s.storageV1MomentMz)
  writeNpyFloat64(root / "moments" / "dispersity.npy", s.storageV1MomentDispersity)
  writeNpyUint64(root / "channel_events" / "snapshot_id.npy", s.storageV1ChannelSnapshotIds)
  writeNpyUint32(root / "channel_events" / "channel_id.npy", s.storageV1ChannelIds)
  writeNpyUint64(root / "channel_events" / "event_count.npy", s.storageV1ChannelEventCounts)
  writeNpyUint64(root / "channel_events" / "productive_event_count.npy", s.storageV1ChannelProductiveCounts)
  writeNpyUint64(root / "channel_events" / "nonproductive_event_count.npy", s.storageV1ChannelNonproductiveCounts)
  writeNpyUint64(root / "diagnostics" / "channel_trace" / "kmc_event.npy", s.storageTraceKmcEvents)
  writeNpyFloat64(root / "diagnostics" / "channel_trace" / "time.npy", s.storageTraceTimes)
  writeNpyFloat64(root / "diagnostics" / "channel_trace" / "dt.npy", s.storageTraceDt)
  writeNpyUint32(root / "diagnostics" / "channel_trace" / "channel_id.npy", s.storageTraceChannelIds)
  writeNpyFloat64(root / "diagnostics" / "channel_trace" / "rate.npy", s.storageTraceRates)
  writeNpyFloat64(root / "diagnostics" / "channel_trace" / "propensity.npy", s.storageTracePropensities)
  writeNpyFloat64(root / "diagnostics" / "channel_trace" / "total_propensity.npy", s.storageTraceTotalPropensities)
  writeNpyUint64(root / "actions" / "action_id.npy", s.storageV1ActionIds)
  writeNpyUint64(root / "actions" / "kmc_event.npy", s.storageV1ActionKmcEvents)
  writeNpyFloat64(root / "actions" / "time.npy", s.storageV1ActionTimes)
  writeNpyUint32(root / "actions" / "source_line.npy", s.storageV1ActionSourceLines)
  writeNpyUint32(root / "actions" / "action_type_id.npy", s.storageV1ActionTypeIds)
  writeNpyUint32(root / "actions" / "trigger_type_id.npy", s.storageV1ActionTriggerTypeIds)
  writeNpyFloat64(root / "actions" / "scheduled_time.npy", s.storageV1ActionScheduledTimes)
  writeNpyUint32(root / "actions" / "target_id.npy", s.storageV1ActionTargetIds)
  writeNpyFloat64(root / "actions" / "requested_value.npy", s.storageV1ActionRequestedValues)
  writeNpyFloat64(root / "actions" / "before_value.npy", s.storageV1ActionBeforeValues)
  writeNpyFloat64(root / "actions" / "after_value.npy", s.storageV1ActionAfterValues)
  writeNpyBool(root / "actions" / "state_changed.npy", s.storageV1ActionStateChanged)
  writeNpyBool(root / "actions" / "output_written.npy", s.storageV1ActionOutputWritten)
  writeNpyBool(root / "actions" / "has_snapshot.npy", s.storageV1ActionHasSnapshot)
  writeNpyUint64(root / "actions" / "snapshot_id.npy", s.storageV1ActionSnapshotIds)
  writeNpyBool(root / "actions" / "has_kinetic_parameter_set.npy", s.storageV1ActionHasKineticSet)
  writeNpyUint64(root / "actions" / "kinetic_parameter_set_id.npy", s.storageV1ActionKineticSetIds)
  writeNpyUint64(root / "feed_events" / "action_id.npy", s.storageV1FeedActionIds)
  writeNpyUint32(root / "feed_events" / "feed_id.npy", s.storageV1FeedIds)
  writeNpyFloat64(root / "feed_events" / "dose_mL.npy", s.storageV1FeedDoseMl)
  writeNpyFloat64(root / "feed_events" / "volume_before_mL.npy", s.storageV1FeedVolumeBeforeMl)
  writeNpyFloat64(root / "feed_events" / "volume_after_mL.npy", s.storageV1FeedVolumeAfterMl)
  writeNpyFloat64(root / "feed_events" / "kmc_volume_before_L.npy", s.storageV1FeedKmcVolumeBeforeL)
  writeNpyFloat64(root / "feed_events" / "kmc_volume_after_L.npy", s.storageV1FeedKmcVolumeAfterL)
  writeNpyUint64(root / "monomer_balance" / "snapshot_id.npy", s.storageV1MonomerBalanceSnapshotIds)
  writeNpyUint32(root / "monomer_balance" / "monomer_id.npy", s.storageV1MonomerBalanceMonomerIds)
  writeNpyFloat64(root / "monomer_balance" / "initial_moles.npy", s.storageV1MonomerInitialMoles)
  writeNpyFloat64(root / "monomer_balance" / "introduced_moles.npy", s.storageV1MonomerIntroducedMoles)
  writeNpyFloat64(root / "monomer_balance" / "free_moles.npy", s.storageV1MonomerFreeMoles)
  writeNpyFloat64(root / "monomer_balance" / "incorporated_moles.npy", s.storageV1MonomerIncorporatedMoles)
  writeNpyFloat64(root / "monomer_balance" / "conversion.npy", s.storageV1MonomerConversion)
  writeNpyUint64(root / "species_balance" / "snapshot_id.npy", s.storageV1SpeciesBalanceSnapshotIds)
  writeNpyUint32(root / "species_balance" / "entity_id.npy", s.storageV1SpeciesBalanceEntityIds)
  writeNpyFloat64(root / "species_balance" / "initial_moles.npy", s.storageV1SpeciesInitialMoles)
  writeNpyFloat64(root / "species_balance" / "dosed_moles.npy", s.storageV1SpeciesDosedMoles)
  writeNpyFloat64(root / "species_balance" / "total_moles.npy", s.storageV1SpeciesTotalMoles)
  writeNpyFloat64(root / "species_balance" / "free_moles.npy", s.storageV1SpeciesFreeMoles)
  writeNpyFloat64(root / "species_balance" / "consumed_moles.npy", s.storageV1SpeciesConsumedMoles)
  writeNpyUint64(root / "action_conditions" / "condition_record_id.npy", s.storageV1ConditionRecordIds)
  writeNpyUint64(root / "action_conditions" / "action_id.npy", s.storageV1ConditionActionIds)
  writeNpyUint32(root / "action_conditions" / "condition_index.npy", s.storageV1ConditionIndexes)
  writeNpyUint32(root / "action_conditions" / "observable_id.npy", s.storageV1ConditionObservableIds)
  writeNpyUint32(root / "action_conditions" / "operator_id.npy", s.storageV1ConditionOperatorIds)
  writeNpyFloat64(root / "action_conditions" / "threshold.npy", s.storageV1ConditionThresholds)
  writeNpyFloat64(root / "action_conditions" / "observed_value.npy", s.storageV1ConditionObservedValues)
  writeNpyBool(root / "action_conditions" / "condition_met.npy", s.storageV1ConditionMet)
  var messageRecords: seq[JsonNode] = @[]
  for i, message in s.storageV1ActionMessages:
    if message.len > 0:
      messageRecords.add %*{"action_id": $s.storageV1ActionIds[i], "message": message}
  writeJsonlAtomic(root / "actions" / "messages.jsonl", messageRecords)
  writeNpyUint64(root / "kinetic_parameters" / "sets" / "kinetic_parameter_set_id.npy", s.storageV1KineticSetIds)
  writeNpyUint64(root / "kinetic_parameters" / "sets" / "start_kmc_event.npy", s.storageV1KineticStartEvents)
  writeNpyFloat64(root / "kinetic_parameters" / "sets" / "start_time.npy", s.storageV1KineticStartTimes)
  writeNpyBool(root / "kinetic_parameters" / "sets" / "has_source_action.npy", s.storageV1KineticHasSourceAction)
  writeNpyUint64(root / "kinetic_parameters" / "sets" / "source_action_id.npy", s.storageV1KineticSourceActionIds)
  writeNpyUint64(root / "kinetic_parameters" / "values" / "kinetic_parameter_set_id.npy", s.storageV1KineticValueSetIds)
  writeNpyUint32(root / "kinetic_parameters" / "values" / "kinetic_parameter_id.npy", s.storageV1KineticParameterIds)
  writeNpyFloat64(root / "kinetic_parameters" / "values" / "value.npy", s.storageV1KineticValues)
  let validation = validateStorageV1(m, s, status)
  writeJsonlAtomic(root / "diagnostics" / "validation.jsonl", validation.records)
  let validationStatus = if validation.errorCount == 0: "passed" else: "failed"
  let effectiveStatus = if status == rsCompleted and validation.errorCount > 0: rsFailed else: status
  let effectiveExitCode = if effectiveStatus == rsFailed and exitCode == 0: 1 else: exitCode
  let storageManifest = writeStorageManifest(root)
  writePrettyJsonAtomic(root / "run_metadata.json",
    metadataNode(m, s, effectiveStatus, startedAt, finishedAt, wallTime, effectiveExitCode, terminationReason,
      validationStatus, validation.warningCount, validation.errorCount, storageManifest))
  if effectiveStatus == rsCompleted:
    writeTextAtomic(root / "RESULTS_COMPLETE", "slimmc-storage-v1\n")
    if dirExists(root / ".work"):
      removeDir(root / ".work")
  elif fileExists(root / "RESULTS_COMPLETE"):
    removeFile(root / "RESULTS_COMPLETE")
