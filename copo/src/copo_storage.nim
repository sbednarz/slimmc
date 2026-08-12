import std/[os, json, times, strutils, math, tables, algorithm]
import copo_types
import storage/copo_schema
import copo_sequence
import copo_stats
import ../../common/[npy_writer, result_json, results_types, sha256_file, storage_manifest, build_provenance]

proc utcNow(): string = now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")
proc reasonId(reason: string): uint32 =
  case reason
  of "initial": 0'u32
  of "scheduled": 1'u32
  of "action", "save", "save_chains": 2'u32
  of "final": 4'u32
  else: 3'u32

proc runId(m: Model): string =
  if m.modelStem.len > 0: m.modelStem else: "run_000001"


proc jsonArgs(args: seq[string]): JsonNode =
  result = newJArray()
  for arg in args: result.add %arg

proc copoConditionNode(m: Model; c: AtomicCondition): JsonNode =
  result = newJObject()
  result["observable"] = %conditionalObservableName(m, c)
  result["comparison"] = %comparisonName(c.comparison)
  result["threshold"] = %c.threshold

proc copoResolvedModelNode(m: Model): JsonNode =
  result = newJObject()
  result["schema"] = %"slimmc-resolved-model-v1"
  result["kinetic_model"] = %"copo"
  result["desc"] = %(if m.hasDescription: m.description else: "")

  result["parameters"] = %*{
    "kmc_volume_L": (if m.hasInitVolume and m.currentVolumeMl > 0.0: m.V * m.initVolumeMl / m.currentVolumeMl else: m.V),
    "init_volume_mL": (if m.hasInitVolume: m.initVolumeMl else: 0.0),
    "initial_temperature_K": m.T,
    "t_end_s": m.t_end,
    "max_events": m.max_steps,
    "when_check_events": m.whenCheckEvents,
    "seed": m.seed,
    "dp_max": m.dp_max,
    "sequence_mode": m.sequence_mode,
    "mass_model": (if m.mass_model == mmRepeatUnits: "repeat_units" else: "with_end_groups")
  }

  result["memory_policy"] = %*{
    "has_limit": m.memoryPolicy.hasLimit,
    "limit_bytes": m.memoryPolicy.limitBytes,
    "save_on_limit": m.memoryPolicy.snapshotOnLimit,
    "stop_on_limit": m.memoryPolicy.stopOnLimit
  }

  var monomers = newJArray()
  for item in m.monomers:
    monomers.add %*{"name": item.name, "initial_concentration_mol_L": item.c0, "molar_mass_g_mol": item.mw}
  result["monomers"] = monomers

  var species = newJArray()
  for item in m.species:
    species.add %*{"name": item.name, "initial_concentration_mol_L": item.c0}
  result["species"] = species

  var feeds = newJArray()
  for feed in m.feeds:
    var components = newJArray()
    for i, c in feed.monomerConcentrations:
      if c > 0.0:
        components.add %*{"name": m.monomers[i].name, "kind": "monomer", "concentration_mol_L": c}
    for i, c in feed.speciesConcentrations:
      if c > 0.0:
        components.add %*{"name": m.species[i].name, "kind": "species", "concentration_mol_L": c}
    feeds.add %*{"name": feed.name, "components": components}
  result["feeds"] = feeds

  var polymers = newJArray()
  for i, pool in m.pools:
    var node = newJObject()
    node["name"] = %pool.name
    node["population_activity"] = %(if pool.kind == pkActive: "live" else: "dead")
    if i < m.poolTerminalMer.len and m.poolTerminalMer[i] >= 0:
      node["terminal_monomer"] = %m.monomers[m.poolTerminalMer[i]].name
    else: node["terminal_monomer"] = newJNull()
    if i < m.poolPenultimateMer.len and m.poolPenultimateMer[i] >= 0:
      node["penultimate_monomer"] = %m.monomers[m.poolPenultimateMer[i]].name
    else: node["penultimate_monomer"] = newJNull()
    polymers.add node
  result["polymers"] = polymers

  var endgroups = newJArray()
  for item in m.endgroups:
    endgroups.add %*{"name": item.name, "molar_mass_g_mol": item.mw, "builtin": false}
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
    node["efficiency"] = %channel.efficiency
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
    for condition in action.conditions: conditions.add copoConditionNode(m, condition)
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
    "chain_dp_dtype": "int32",
    "chain_dp_max": $high(int32),
    "dead_pool": (if m.deadPoolId >= 0: m.pools[m.deadPoolId].name else: "")
  }

proc metadataNode(m: Model; s: State; status: string; finishedAt = ""; wallTime = 0.0; exitCode = 0; terminationReason = ""; validationStatus="not_run"; validationErrors=0; storageManifest: StorageManifest = StorageManifest()): JsonNode =
  result = newJObject()
  result["run_id"] = %runId(m)
  result["input_model_file"] = %"input.model"
  result["source_model_name"] = %extractFilename(m.modelFile)
  let inputPath = s.storageV1RunDir / "input.model"
  let schemaPath = s.storageV1RunDir / "schema.jsonl"
  if fileExists(inputPath): result["input_model_sha256"] = %sha256File(inputPath)
  result["storage"] = %StorageName
  result["storage_format_version"] = %StorageFormatVersion
  if fileExists(schemaPath): result["schema_sha256"] = %sha256File(schemaPath)
  result["seed"] = %($uint64(m.seed))
  result["engine"] = %AppName
  result["cli_version"] = %AppVersion
  result["engine_version"] = %AppVersion
  result["started_at_utc"] = %s.storageV1StartedAt
  if finishedAt.len > 0:
    result["finished_at_utc"] = %finishedAt
    result["wall_time_s"] = %wallTime
  result["run_status"] = %status
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
  result["kinetic_model"] = %"copo"
  result["n_monomers"] = %m.monomers.len
  result["dp_max"] = %($m.dp_max)
  result["sequence_mode"] = %m.sequence_mode
  result["engine_chain_dp_dtype"] = %"int32"
  result["engine_chain_dp_max"] = %($high(int32))
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
    "started_at_utc": s.storageV1StartedAt,
    "status": status,
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
  result["validation_warning_count"] = %0
  result["validation_error_count"] = %validationErrors
  result["channel_trace_enabled"] = %(s.storageTraceKmcEvents.len > 0 or s.channelTraceRowsWritten > 0)
  result["channel_trace_complete"] = %(not s.channelTraceTruncated)
  result["channel_trace_rows"] = %s.channelTraceRowsWritten
  result["channel_trace_truncated"] = %s.channelTraceTruncated
  var variables = newJArray()
  for v in m.variables:
    variables.add(%*{"kind": v.kind, "name": v.name, "value": v.value, "unit": v.unit})
  result["variables"] = variables
  result["model"] = copoResolvedModelNode(m)
  var storageInfo = %*{
    "name": StorageName,
    "format_version": StorageFormatVersion,
    "complete": status == "completed" and validationErrors == 0,
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


proc ensureDirs(root: string) =
  for d in [root, root/".work", root/"snapshots", root/"state", root/"channel_events", root/"channel_propensities", root/"actions", root/"action_conditions", root/"feed_events", root/"monomer_balance", root/"species_balance", root/"kinetic_parameters", root/"kinetic_parameters"/"sets", root/"kinetic_parameters"/"values", root/"chains", root/"chain_composition", root/"sequences", root/"moments", root/"memory", root/"microstructure_motifs", root/"block_statistics", root/"diagnostics"]: createDir(d)

proc captureKineticParameterSetV1*(m: Model; s: var State; hasSourceAction: bool; sourceActionId: uint64 = 0'u64) =
  let setId = uint64(s.storageV1KineticSetIds.len)
  s.storageV1KineticSetIds.add setId
  s.storageV1KineticStartEvents.add uint64(s.kmcEvent)
  s.storageV1KineticStartTimes.add s.t
  s.storageV1KineticHasSourceAction.add hasSourceAction
  s.storageV1KineticSourceActionIds.add(if hasSourceAction: sourceActionId else: 0'u64)
  var pid=0'u32
  template addv(v:float64)=
    s.storageV1KineticValueSetIds.add setId; s.storageV1KineticParameterIds.add pid; s.storageV1KineticValues.add v; inc pid
  addv(m.T)
  for rid,r in m.rates:
    addv(m.rateValue(rid))
    if r.declaredArrhenius:
      addv(r.Apre); addv(r.Ea)

proc parseMaybeFloat(x:string):float64 =
  if x.len==0: return NaN
  try: parseFloat(x)
  except ValueError: NaN

proc conditionObservableId(m:Model;a:AtomicCondition):uint32 =
  case a.observable
  of woTotalConversion: 1'u32
  of woMonomerConversion: uint32(2+a.targetId)
  of woSpeciesConc: uint32(2+m.monomers.len+a.targetId)
  of woMonomerConc: uint32(2+m.monomers.len+m.species.len+a.targetId)

proc conditionOperatorId(a:AtomicCondition):uint32 = (if a.comparison==coGreater:1'u32 else:2'u32)
proc triggerId(x:string):uint32 = (if x=="at":1'u32 elif x=="every":2'u32 elif x=="when":3'u32 else:0'u32)
proc actionTargetId(m:Model; action:ActionKind; args:seq[string]):uint32 =
  if action in {eaSetTemp,eaAddTemp}: return 1'u32
  if args.len==0: return 0'u32
  if action in {eaSetK,eaAddK} and m.rateByName.hasKey(args[0]): return uint32(2+m.rateByName[args[0]])
  let base=2+m.rates.len
  if m.monomerByName.hasKey(args[0]): return uint32(base+m.monomerByName[args[0]])
  if m.speciesByName.hasKey(args[0]): return uint32(base+m.monomers.len+m.speciesByName[args[0]])
  0'u32

proc captureActionV1*(m:Model;s:var State;lineNo:int;trigger:string;scheduledTime:string;
                      conditions:seq[AtomicCondition]; observed:seq[float64]; action:ActionKind;
                      args:seq[string]; beforeText,afterText:string; stateChanged,outputWritten:bool;
                      message:string="") =
  let aid=uint64(s.storageV1ActionIds.len)
  s.storageV1ActionIds.add aid; s.storageV1ActionKmcEvents.add uint64(s.kmcEvent); s.storageV1ActionTimes.add s.t
  s.storageV1ActionSourceLines.add uint32(lineNo); s.storageV1ActionTypeIds.add uint32(ord(action)); s.storageV1ActionTriggerTypeIds.add triggerId(trigger)
  s.storageV1ActionScheduledTimes.add parseMaybeFloat(scheduledTime); s.storageV1ActionTargetIds.add actionTargetId(m,action,args)
  s.storageV1ActionRequestedValues.add(if args.len>0: parseMaybeFloat(args[^1]) else: NaN)
  s.storageV1ActionBeforeValues.add parseMaybeFloat(beforeText); s.storageV1ActionAfterValues.add parseMaybeFloat(afterText)
  s.storageV1ActionStateChanged.add stateChanged; s.storageV1ActionOutputWritten.add outputWritten
  let hs=action in {eaSave,eaSaveChains,eaSetK,eaAddK,eaSetTemp,eaAddTemp}
  let hk=action in {eaSetK,eaAddK,eaSetTemp,eaAddTemp}
  s.storageV1ActionHasSnapshot.add hs; s.storageV1ActionSnapshotIds.add(if hs and s.storageV1SnapshotIds.len>0:s.storageV1SnapshotIds[^1] else:0'u64)
  s.storageV1ActionHasKineticSet.add hk; s.storageV1ActionKineticSetIds.add(if hk and s.storageV1KineticSetIds.len>0:s.storageV1KineticSetIds[^1] else:0'u64)
  s.storageV1ActionMessages.add message
  if action == eaFeed and args.len in [2, 3]:
    let fid = m.feedByName[args[0]]
    let rawDose = parseMaybeFloat(args[1])
    let doseMl = if args.len == 2: rawDose * 1000.0 else: rawDose
    let afterMl = m.currentVolumeMl
    let beforeMl = afterMl - doseMl
    let kmcAfter = m.V
    let kmcBefore = if afterMl > 0.0: kmcAfter * beforeMl / afterMl else: NaN
    s.storageV1FeedActionIds.add aid
    s.storageV1FeedIds.add uint32(fid)
    s.storageV1FeedDoseMl.add doseMl
    s.storageV1FeedVolumeBeforeMl.add beforeMl
    s.storageV1FeedVolumeAfterMl.add afterMl
    s.storageV1FeedKmcVolumeBeforeL.add kmcBefore
    s.storageV1FeedKmcVolumeAfterL.add kmcAfter
  for i,c in conditions:
    let cid=uint64(s.storageV1ConditionRecordIds.len)
    s.storageV1ConditionRecordIds.add cid; s.storageV1ConditionActionIds.add aid; s.storageV1ConditionIndexes.add uint32(i)
    s.storageV1ConditionObservableIds.add conditionObservableId(m,c); s.storageV1ConditionOperatorIds.add conditionOperatorId(c)
    s.storageV1ConditionThresholds.add c.threshold; s.storageV1ConditionObservedValues.add observed[i]; s.storageV1ConditionMet.add true

proc initStorageV1*(m: Model; s: var State) =
  s.storageV1RunDir = m.output_dir
  s.storageV1StartedAt = utcNow()
  ensureDirs(s.storageV1RunDir)
  if fileExists(s.storageV1RunDir / "RESULTS_COMPLETE"): removeFile(s.storageV1RunDir / "RESULTS_COMPLETE")
  if fileExists(m.modelFile): copyFile(m.modelFile, s.storageV1RunDir / "input.model")
  writeJsonlAtomic(s.storageV1RunDir / "schema.jsonl", schemaRecords(m))
  writePrettyJsonAtomic(s.storageV1RunDir / "run_metadata.json", metadataNode(m,s,"running") )
  writeFile(s.storageV1RunDir / "diagnostics" / "run.log", AppName & " Slimmc Storage v1 run log\n")
  captureKineticParameterSetV1(m,s,false)

proc addState(m: Model; s: var State; sid:uint64) =
  var eid=0'u32
  template addCount(n:int64) =
    let u=uint64(n)
    s.storageV1StateSnapshotIds.add sid; s.storageV1StateEntityIds.add eid; s.storageV1StateCounts.add u
    let mol=float64(u)/AvogadroConstantMolInv
    s.storageV1StateMoles.add mol; s.storageV1StateConcentrations.add mol/m.V; inc eid
  for n in s.monomerN: addCount(n)
  for n in s.speciesN: addCount(n)
  var live:int64=0
  for pool in s.livePools: live += int64(pool.len)
  var dead:int64=0
  for x in s.deadChains: dead += x.count
  addCount(live); addCount(dead)


proc endId(s: var State; name: string): uint32 =
  if name.len == 0: return 0'u32
  for i, known in s.storageV1ObservedEndGroups:
    if known == name: return uint32(i + 1)
  s.storageV1ObservedEndGroups.add name
  uint32(s.storageV1ObservedEndGroups.len)

proc addMomentRows(m: Model; s: var State; sid: uint64; firstRow: int) =
  let lastRow = s.storageV1ChainRecordIds.len
  for scope in 0'u32..2'u32:
    for basis in 0'u32..1'u32:
      var n=0'u64; var sd=0.0; var sd2=0.0; var sm=0.0; var sm2=0.0; var sm3=0.0
      for i in firstRow..<lastRow:
        let pop=s.storageV1ChainPopulationIds[i]
        if scope==1'u32 and pop!=0'u32: continue
        if scope==2'u32 and pop!=1'u32: continue
        let c=s.storageV1ChainCounts[i]; let d=float64(s.storageV1ChainDp[i])
        var mass=s.storageV1ChainMolarMass[i]
        if basis==0'u32:
          mass=0.0
          let base=i*m.monomers.len
          for j,mo in m.monomers: mass += float64(s.storageV1CompositionUnitCounts[base+j]) * mo.mw
        n += c; sd += float64(c)*d; sd2 += float64(c)*d*d
        sm += float64(c)*mass; sm2 += float64(c)*mass*mass; sm3 += float64(c)*mass*mass*mass
      s.storageV1MomentSnapshotIds.add sid; s.storageV1MomentPopulationScopeIds.add scope; s.storageV1MomentMassBasisIds.add basis
      s.storageV1MomentChainCounts.add n; s.storageV1MomentSumDp.add sd; s.storageV1MomentSumDp2.add sd2
      if n==0: s.storageV1MomentDpN.add NaN; s.storageV1MomentDpW.add NaN; s.storageV1MomentMn.add NaN; s.storageV1MomentMw.add NaN; s.storageV1MomentMz.add NaN; s.storageV1MomentDispersity.add NaN
      else:
        let dpn=sd/float64(n); let dpw=(if sd>0:sd2/sd else:NaN); let mn=sm/float64(n); let mw=(if sm>0:sm2/sm else:NaN); let mz=(if sm2>0:sm3/sm2 else:NaN)
        s.storageV1MomentDpN.add dpn; s.storageV1MomentDpW.add dpw; s.storageV1MomentMn.add mn; s.storageV1MomentMw.add mw; s.storageV1MomentMz.add mz; s.storageV1MomentDispersity.add(if mn>0:mw/mn else:NaN)
      s.storageV1MomentSumMass.add sm; s.storageV1MomentSumMass2.add sm2; s.storageV1MomentSumMass3.add sm3

proc addChainRows(m: Model; s: var State; sid: uint64) =
  let firstRow=s.storageV1ChainRecordIds.len
  template addCommon(pop:uint32; pool:uint32; origin:uint32; dpv:int32; massv:float64; countv:int64; left,right:string; first,pen,last:int; seqData:untyped; hasSeq:bool; nm:seq[int32]) =
    let rid=uint64(s.storageV1ChainRecordIds.len); let cnt=uint64(countv)
    s.storageV1ChainRecordIds.add rid; s.storageV1ChainSnapshotIds.add sid; s.storageV1ChainPopulationIds.add pop; s.storageV1ChainPoolIds.add pool; s.storageV1ChainOriginIds.add origin
    s.storageV1ChainDp.add uint64(dpv); s.storageV1ChainMolarMass.add massv; s.storageV1ChainCounts.add cnt
    let mol=float64(cnt)/AvogadroConstantMolInv; s.storageV1ChainMoles.add mol; s.storageV1ChainConcentrations.add mol/m.V
    s.storageV1ChainLeftEndIds.add endId(s,left); s.storageV1ChainRightEndIds.add endId(s,right)
    s.storageV1ChainHasFirst.add first>=0; s.storageV1ChainFirstIds.add(if first>=0:uint32(first) else:0'u32)
    s.storageV1ChainHasPenultimate.add pen>=0; s.storageV1ChainPenultimateIds.add(if pen>=0:uint32(pen) else:0'u32)
    s.storageV1ChainHasLast.add last>=0; s.storageV1ChainLastIds.add(if last>=0:uint32(last) else:0'u32)
    let off=uint64(s.storageV1SequenceSymbols.len); var slen=0'u64
    if hasSeq:
      for x in seqData: s.storageV1SequenceSymbols.add uint32(x); inc slen
    s.storageV1ChainHasSequence.add hasSeq; s.storageV1ChainSequenceOffsets.add(if hasSeq:off else:0'u64); s.storageV1ChainSequenceLengths.add(slen)
    for mi,x in nm: s.storageV1CompositionRecordIds.add rid; s.storageV1CompositionMonomerIds.add uint32(mi); s.storageV1CompositionUnitCounts.add uint64(x)
  for pi,p in s.livePools:
    for c in p:
      let hs=m.sequence_mode=="full"
      addCommon(0'u32,uint32(pi),uint32(ord(c.formedBy)),c.dp,c.mass,1,c.left_end,c.right_end,c.mers.first(),c.mers.prev(),c.mers.last(),c.mers.data,hs,c.nMer)
  for d in s.deadChains:
    var seqIds: seq[uint8] = @[]
    let hs=m.sequence_mode=="full" and d.sequenceStored
    if hs and d.sequenceText.len>0:
      for name in d.sequenceText.split('|'):
        if m.monomerByName.hasKey(name): seqIds.add uint8(m.monomerByName[name])
    addCommon(1'u32,uint32(m.deadPoolId),uint32(ord(d.formedBy)),d.dp,d.mass,d.count,d.left_end,d.right_end,d.firstMer,d.penultimateMer,d.lastMer,seqIds,hs,d.nMer)
  addMomentRows(m,s,sid,firstRow)

proc addMemoryRow(m: Model; s: var State; sid: uint64) =
  var live=0'u64; var deadRecords=uint64(s.deadChains.len); var dead=0'u64; var liveMers=0'u64; var deadMers=0'u64
  for p in s.livePools:
    live += uint64(p.len)
    for c in p: liveMers += uint64(c.dp)
  for d in s.deadChains:
    dead += uint64(d.count)
    if d.sequenceStored: deadMers += uint64(d.dp) * uint64(d.count)
  let liveSeq=liveMers; let liveObj=live*96'u64
  var blockCounters=0'u64
  for d in s.deadChains: blockCounters += uint64(d.blockCounts.len)
  let deadBytes=deadRecords*uint64(128+16*m.monomers.len)+blockCounters*24'u64
  s.storageV1MemorySnapshotIds.add sid; s.storageV1MemoryLiveChains.add live; s.storageV1MemoryDeadRecords.add deadRecords; s.storageV1MemoryDeadChains.add dead
  s.storageV1MemoryStoredLiveMers.add liveMers; s.storageV1MemoryStoredDeadMers.add deadMers; s.storageV1MemoryLiveSeqBytes.add liveSeq; s.storageV1MemoryLiveObjectBytes.add liveObj; s.storageV1MemoryDeadRecordBytes.add deadBytes; s.storageV1MemoryTotalEstBytes.add(liveSeq+liveObj+deadBytes)

proc addMicrostructureRows(m: Model; s: var State; sid: uint64) =
  ## Persist aggregate motif/block information independently of literal sequences.
  let dy = globalDyads(m, s)
  for motifId, value in dy:
    s.storageV1MotifSnapshotIds.add sid
    s.storageV1MotifOrders.add 2'u32
    s.storageV1MotifIds.add uint32(motifId)
    s.storageV1MotifCounts.add uint64(value)
  let tr = globalTriads(m, s)
  for motifId, value in tr:
    s.storageV1MotifSnapshotIds.add sid
    s.storageV1MotifOrders.add 3'u32
    s.storageV1MotifIds.add uint32(motifId)
    s.storageV1MotifCounts.add uint64(value)
  for bc in globalBlockCounts(m, s):
    if bc.count <= 0 or bc.length <= 0: continue
    s.storageV1BlockSnapshotIds.add sid
    s.storageV1BlockMonomerIds.add uint32(bc.monomerId)
    s.storageV1BlockLengths.add uint64(bc.length)
    s.storageV1BlockCounts.add uint64(bc.count)

proc captureStorageV1Snapshot*(m: Model; s: var State; reason:string; isFinal=false; hasChains=false; propensities: seq[float64] = @[]) =
  if s.storageV1SnapshotIds.len > 0:
    let i=s.storageV1SnapshotIds.high
    if s.storageV1KmcEvents[i]==uint64(s.kmcEvent) and timeClose(s.storageV1Times[i], s.t) and
       s.storageV1StateRevisions[i] == uint64(s.stateRevision):
      if isFinal: s.storageV1IsFinal[i]=true
      if hasChains and not s.storageV1HasChains[i]:
        s.storageV1HasChains[i]=true
        s.storageV1HasSequences[i]=(m.sequence_mode=="full")
        addChainRows(m,s,uint64(i)); addMemoryRow(m,s,uint64(i))
      return
  let sid=uint64(s.storageV1SnapshotIds.len)
  s.storageV1SnapshotIds.add sid; s.storageV1StateRevisions.add uint64(s.stateRevision); s.storageV1Times.add s.t; s.storageV1KmcEvents.add uint64(s.kmcEvent)
  s.storageV1ReasonIds.add reasonId(reason); s.storageV1IsFinal.add isFinal
  s.storageV1HasChains.add hasChains; s.storageV1HasSequences.add(hasChains and m.sequence_mode=="full")
  s.storageV1ParameterSetIds.add(if s.storageV1KineticSetIds.len>0: s.storageV1KineticSetIds[^1] else: 0'u64)
  s.storageV1VolumeMl.add(if m.hasInitVolume: m.currentVolumeMl else: NaN)
  s.storageV1KmcVolumeL.add m.V
  var liveCount = 0'u64
  for pool in s.livePools:
    liveCount += uint64(pool.len)
  var deadCount = 0'u64
  for chain in s.deadChains:
    if chain.count > 0:
      deadCount += uint64(chain.count)
  s.storageV1ChainCountLive.add liveCount
  s.storageV1ChainCountDead.add deadCount
  s.storageV1ChainCountTotal.add liveCount + deadCount
  if m.hasInitVolume:
    let physicalScale = (m.currentVolumeMl / 1000.0) / m.V
    var balanceId = 0'u32
    for mid in 0..<m.monomers.len:
      let initial = m.monomers[mid].c0 * (m.initVolumeMl / 1000.0)
      let dosed = float64(s.monomerDosedN[mid]) / AvogadroConstantMolInv * physicalScale
      let total = float64(s.monomerExternalN[mid]) / AvogadroConstantMolInv * physicalScale
      let free = float64(s.monomerN[mid]) / AvogadroConstantMolInv * physicalScale
      s.storageV1SpeciesBalanceSnapshotIds.add sid
      s.storageV1SpeciesBalanceEntityIds.add balanceId; inc balanceId
      s.storageV1SpeciesInitialMoles.add initial
      s.storageV1SpeciesDosedMoles.add dosed
      s.storageV1SpeciesTotalMoles.add total
      s.storageV1SpeciesFreeMoles.add free
      s.storageV1SpeciesConsumedMoles.add total - free
      let introduced = float64(s.monomerN0[mid]) / AvogadroConstantMolInv * physicalScale
      let incorporated = introduced - free
      s.storageV1MonomerBalanceSnapshotIds.add sid
      s.storageV1MonomerBalanceMonomerIds.add uint32(mid)
      s.storageV1MonomerInitialMoles.add initial
      s.storageV1MonomerIntroducedMoles.add introduced
      s.storageV1MonomerFreeMoles.add free
      s.storageV1MonomerIncorporatedMoles.add incorporated
      s.storageV1MonomerConversion.add(if introduced > 0.0: incorporated / introduced else: 0.0)
    for sp in 0..<m.species.len:
      let initial = m.species[sp].c0 * (m.initVolumeMl / 1000.0)
      let dosed = float64(s.speciesDosedN[sp]) / AvogadroConstantMolInv * physicalScale
      let total = float64(s.speciesExternalN[sp]) / AvogadroConstantMolInv * physicalScale
      let free = float64(s.speciesN[sp]) / AvogadroConstantMolInv * physicalScale
      s.storageV1SpeciesBalanceSnapshotIds.add sid
      s.storageV1SpeciesBalanceEntityIds.add balanceId; inc balanceId
      s.storageV1SpeciesInitialMoles.add initial
      s.storageV1SpeciesDosedMoles.add dosed
      s.storageV1SpeciesTotalMoles.add total
      s.storageV1SpeciesFreeMoles.add free
      s.storageV1SpeciesConsumedMoles.add total - free
  addState(m,s,sid)
  addMicrostructureRows(m,s,sid)
  var totalPropensity = 0.0
  if propensities.len notin [0, m.channels.len]:
    raise newException(ValueError, "snapshot propensity vector length mismatch")
  for cid in 0..<m.channels.len:
    let propensity = if propensities.len == m.channels.len: propensities[cid] else: NaN
    s.storageV1PropensitySnapshotIds.add sid
    s.storageV1PropensityChannelIds.add uint32(cid)
    s.storageV1PropensityValues.add propensity
    if propensity.classify notin {fcNan, fcInf, fcNegInf}: totalPropensity += propensity
    s.storageV1TotalPropensities.add 0.0 # filled below after total is known
    s.storageV1ChannelSnapshotIds.add sid; s.storageV1ChannelIds.add uint32(cid)
    let fires = uint64(s.channelFires[cid])
    s.storageV1ChannelEventCounts.add fires
    if m.channels[cid].kind in {chRxnUni,chRxnBiDiff,chRxnBiSame} and m.channels[cid].efficiency < 1.0 - Eps:
      s.storageV1ChannelProductiveCounts.add uint64(s.channelSuccesses[cid])
      s.storageV1ChannelNonproductiveCounts.add uint64(s.channelFailures[cid])
    else:
      s.storageV1ChannelProductiveCounts.add fires
      s.storageV1ChannelNonproductiveCounts.add 0'u64
  if m.channels.len > 0:
    for i in s.storageV1TotalPropensities.len - m.channels.len ..< s.storageV1TotalPropensities.len:
      s.storageV1TotalPropensities[i] = totalPropensity
  if hasChains: addChainRows(m,s,sid); addMemoryRow(m,s,sid)

proc writeColumns(s: State) =
  let r=s.storageV1RunDir
  writeNpyUint64(r/"snapshots"/"snapshot_id.npy",s.storageV1SnapshotIds)
  writeNpyFloat64(r/"snapshots"/"time.npy",s.storageV1Times)
  writeNpyUint64(r/"snapshots"/"kmc_event.npy",s.storageV1KmcEvents)
  writeNpyUint32(r/"snapshots"/"snapshot_reason_id.npy",s.storageV1ReasonIds)
  writeNpyBool(r/"snapshots"/"is_final.npy",s.storageV1IsFinal)
  writeNpyBool(r/"snapshots"/"has_chains.npy",s.storageV1HasChains)
  writeNpyBool(r/"snapshots"/"has_sequences.npy",s.storageV1HasSequences)
  writeNpyUint64(r/"snapshots"/"kinetic_parameter_set_id.npy",s.storageV1ParameterSetIds)
  writeNpyFloat64(r/"snapshots"/"volume_mL.npy",s.storageV1VolumeMl)
  writeNpyFloat64(r/"snapshots"/"kmc_volume_L.npy",s.storageV1KmcVolumeL)
  writeNpyUint64(r/"snapshots"/"chain_count_live.npy",s.storageV1ChainCountLive)
  writeNpyUint64(r/"snapshots"/"chain_count_dead.npy",s.storageV1ChainCountDead)
  writeNpyUint64(r/"snapshots"/"chain_count_total.npy",s.storageV1ChainCountTotal)
  writeNpyUint64(r/"state"/"snapshot_id.npy",s.storageV1StateSnapshotIds)
  writeNpyUint32(r/"state"/"entity_id.npy",s.storageV1StateEntityIds)
  writeNpyUint64(r/"state"/"count.npy",s.storageV1StateCounts)
  writeNpyFloat64(r/"state"/"moles.npy",s.storageV1StateMoles)
  writeNpyFloat64(r/"state"/"concentration.npy",s.storageV1StateConcentrations)
  writeNpyUint64(r/"channel_events"/"snapshot_id.npy",s.storageV1ChannelSnapshotIds)
  writeNpyUint32(r/"channel_events"/"channel_id.npy",s.storageV1ChannelIds)
  writeNpyUint64(r/"channel_events"/"event_count.npy",s.storageV1ChannelEventCounts)
  writeNpyUint64(r/"channel_events"/"productive_event_count.npy",s.storageV1ChannelProductiveCounts)
  writeNpyUint64(r/"channel_events"/"nonproductive_event_count.npy",s.storageV1ChannelNonproductiveCounts)
  writeNpyUint64(r/"channel_propensities"/"snapshot_id.npy",s.storageV1PropensitySnapshotIds)
  writeNpyUint32(r/"channel_propensities"/"channel_id.npy",s.storageV1PropensityChannelIds)
  writeNpyFloat64(r/"channel_propensities"/"propensity.npy",s.storageV1PropensityValues)
  writeNpyFloat64(r/"channel_propensities"/"total_propensity.npy",s.storageV1TotalPropensities)
  writeNpyUint64(r/"diagnostics"/"channel_trace"/"kmc_event.npy",s.storageTraceKmcEvents)
  writeNpyFloat64(r/"diagnostics"/"channel_trace"/"time.npy",s.storageTraceTimes)
  writeNpyFloat64(r/"diagnostics"/"channel_trace"/"dt.npy",s.storageTraceDt)
  writeNpyUint32(r/"diagnostics"/"channel_trace"/"channel_id.npy",s.storageTraceChannelIds)
  writeNpyFloat64(r/"diagnostics"/"channel_trace"/"rate.npy",s.storageTraceRates)
  writeNpyFloat64(r/"diagnostics"/"channel_trace"/"propensity.npy",s.storageTracePropensities)
  writeNpyFloat64(r/"diagnostics"/"channel_trace"/"total_propensity.npy",s.storageTraceTotalPropensities)
  writeNpyUint64(r/"actions"/"action_id.npy",s.storageV1ActionIds)
  writeNpyUint64(r/"actions"/"kmc_event.npy",s.storageV1ActionKmcEvents)
  writeNpyFloat64(r/"actions"/"time.npy",s.storageV1ActionTimes)
  writeNpyUint32(r/"actions"/"source_line.npy",s.storageV1ActionSourceLines)
  writeNpyUint32(r/"actions"/"action_type_id.npy",s.storageV1ActionTypeIds)
  writeNpyUint32(r/"actions"/"trigger_type_id.npy",s.storageV1ActionTriggerTypeIds)
  writeNpyFloat64(r/"actions"/"scheduled_time.npy",s.storageV1ActionScheduledTimes)
  writeNpyUint32(r/"actions"/"target_id.npy",s.storageV1ActionTargetIds)
  writeNpyFloat64(r/"actions"/"requested_value.npy",s.storageV1ActionRequestedValues)
  writeNpyFloat64(r/"actions"/"before_value.npy",s.storageV1ActionBeforeValues)
  writeNpyFloat64(r/"actions"/"after_value.npy",s.storageV1ActionAfterValues)
  writeNpyBool(r/"actions"/"state_changed.npy",s.storageV1ActionStateChanged)
  writeNpyBool(r/"actions"/"output_written.npy",s.storageV1ActionOutputWritten)
  writeNpyBool(r/"actions"/"has_snapshot.npy",s.storageV1ActionHasSnapshot)
  writeNpyUint64(r/"actions"/"snapshot_id.npy",s.storageV1ActionSnapshotIds)
  writeNpyBool(r/"actions"/"has_kinetic_parameter_set.npy",s.storageV1ActionHasKineticSet)
  writeNpyUint64(r/"actions"/"kinetic_parameter_set_id.npy",s.storageV1ActionKineticSetIds)
  writeNpyUint64(r/"feed_events"/"action_id.npy",s.storageV1FeedActionIds)
  writeNpyUint32(r/"feed_events"/"feed_id.npy",s.storageV1FeedIds)
  writeNpyFloat64(r/"feed_events"/"dose_mL.npy",s.storageV1FeedDoseMl)
  writeNpyFloat64(r/"feed_events"/"volume_before_mL.npy",s.storageV1FeedVolumeBeforeMl)
  writeNpyFloat64(r/"feed_events"/"volume_after_mL.npy",s.storageV1FeedVolumeAfterMl)
  writeNpyFloat64(r/"feed_events"/"kmc_volume_before_L.npy",s.storageV1FeedKmcVolumeBeforeL)
  writeNpyFloat64(r/"feed_events"/"kmc_volume_after_L.npy",s.storageV1FeedKmcVolumeAfterL)
  writeNpyUint64(r/"monomer_balance"/"snapshot_id.npy",s.storageV1MonomerBalanceSnapshotIds)
  writeNpyUint32(r/"monomer_balance"/"monomer_id.npy",s.storageV1MonomerBalanceMonomerIds)
  writeNpyFloat64(r/"monomer_balance"/"initial_moles.npy",s.storageV1MonomerInitialMoles)
  writeNpyFloat64(r/"monomer_balance"/"introduced_moles.npy",s.storageV1MonomerIntroducedMoles)
  writeNpyFloat64(r/"monomer_balance"/"free_moles.npy",s.storageV1MonomerFreeMoles)
  writeNpyFloat64(r/"monomer_balance"/"incorporated_moles.npy",s.storageV1MonomerIncorporatedMoles)
  writeNpyFloat64(r/"monomer_balance"/"conversion.npy",s.storageV1MonomerConversion)
  writeNpyUint64(r/"species_balance"/"snapshot_id.npy",s.storageV1SpeciesBalanceSnapshotIds)
  writeNpyUint32(r/"species_balance"/"entity_id.npy",s.storageV1SpeciesBalanceEntityIds)
  writeNpyFloat64(r/"species_balance"/"initial_moles.npy",s.storageV1SpeciesInitialMoles)
  writeNpyFloat64(r/"species_balance"/"dosed_moles.npy",s.storageV1SpeciesDosedMoles)
  writeNpyFloat64(r/"species_balance"/"total_moles.npy",s.storageV1SpeciesTotalMoles)
  writeNpyFloat64(r/"species_balance"/"free_moles.npy",s.storageV1SpeciesFreeMoles)
  writeNpyFloat64(r/"species_balance"/"consumed_moles.npy",s.storageV1SpeciesConsumedMoles)
  var msgs: seq[JsonNode] = @[]
  for i,x in s.storageV1ActionMessages:
    if x.len>0: msgs.add %*{"action_id": $(s.storageV1ActionIds[i]), "message": x}
  writeJsonlAtomic(r/"actions"/"messages.jsonl",msgs)
  writeNpyUint64(r/"action_conditions"/"condition_record_id.npy",s.storageV1ConditionRecordIds)
  writeNpyUint64(r/"action_conditions"/"action_id.npy",s.storageV1ConditionActionIds)
  writeNpyUint32(r/"action_conditions"/"condition_index.npy",s.storageV1ConditionIndexes)
  writeNpyUint32(r/"action_conditions"/"observable_id.npy",s.storageV1ConditionObservableIds)
  writeNpyUint32(r/"action_conditions"/"operator_id.npy",s.storageV1ConditionOperatorIds)
  writeNpyFloat64(r/"action_conditions"/"threshold.npy",s.storageV1ConditionThresholds)
  writeNpyFloat64(r/"action_conditions"/"observed_value.npy",s.storageV1ConditionObservedValues)
  writeNpyBool(r/"action_conditions"/"condition_met.npy",s.storageV1ConditionMet)
  writeNpyUint64(r/"kinetic_parameters"/"sets"/"kinetic_parameter_set_id.npy",s.storageV1KineticSetIds)
  writeNpyUint64(r/"kinetic_parameters"/"sets"/"start_kmc_event.npy",s.storageV1KineticStartEvents)
  writeNpyFloat64(r/"kinetic_parameters"/"sets"/"start_time.npy",s.storageV1KineticStartTimes)
  writeNpyBool(r/"kinetic_parameters"/"sets"/"has_source_action.npy",s.storageV1KineticHasSourceAction)
  writeNpyUint64(r/"kinetic_parameters"/"sets"/"source_action_id.npy",s.storageV1KineticSourceActionIds)
  writeNpyUint64(r/"kinetic_parameters"/"values"/"kinetic_parameter_set_id.npy",s.storageV1KineticValueSetIds)
  writeNpyUint32(r/"kinetic_parameters"/"values"/"kinetic_parameter_id.npy",s.storageV1KineticParameterIds)
  writeNpyFloat64(r/"kinetic_parameters"/"values"/"value.npy",s.storageV1KineticValues)


  writeNpyUint64(r/"chains"/"chain_record_id.npy",s.storageV1ChainRecordIds)
  writeNpyUint64(r/"chains"/"snapshot_id.npy",s.storageV1ChainSnapshotIds)
  writeNpyUint32(r/"chains"/"population_id.npy",s.storageV1ChainPopulationIds)
  writeNpyUint32(r/"chains"/"pool_id.npy",s.storageV1ChainPoolIds)
  writeNpyUint32(r/"chains"/"origin_id.npy",s.storageV1ChainOriginIds)
  writeNpyUint64(r/"chains"/"dp.npy",s.storageV1ChainDp)
  writeNpyFloat64(r/"chains"/"molar_mass.npy",s.storageV1ChainMolarMass)
  writeNpyUint64(r/"chains"/"count.npy",s.storageV1ChainCounts)
  writeNpyFloat64(r/"chains"/"moles.npy",s.storageV1ChainMoles)
  writeNpyFloat64(r/"chains"/"concentration.npy",s.storageV1ChainConcentrations)
  writeNpyUint32(r/"chains"/"left_end_id.npy",s.storageV1ChainLeftEndIds); writeNpyUint32(r/"chains"/"right_end_id.npy",s.storageV1ChainRightEndIds)
  writeNpyBool(r/"chains"/"has_first_monomer.npy",s.storageV1ChainHasFirst); writeNpyUint32(r/"chains"/"first_monomer_id.npy",s.storageV1ChainFirstIds)
  writeNpyBool(r/"chains"/"has_penultimate_monomer.npy",s.storageV1ChainHasPenultimate); writeNpyUint32(r/"chains"/"penultimate_monomer_id.npy",s.storageV1ChainPenultimateIds)
  writeNpyBool(r/"chains"/"has_last_monomer.npy",s.storageV1ChainHasLast); writeNpyUint32(r/"chains"/"last_monomer_id.npy",s.storageV1ChainLastIds)
  writeNpyBool(r/"chains"/"has_sequence.npy",s.storageV1ChainHasSequence); writeNpyUint64(r/"chains"/"sequence_offset.npy",s.storageV1ChainSequenceOffsets); writeNpyUint64(r/"chains"/"sequence_length.npy",s.storageV1ChainSequenceLengths)
  writeNpyUint64(r/"chain_composition"/"chain_record_id.npy",s.storageV1CompositionRecordIds); writeNpyUint32(r/"chain_composition"/"monomer_id.npy",s.storageV1CompositionMonomerIds); writeNpyUint64(r/"chain_composition"/"unit_count.npy",s.storageV1CompositionUnitCounts)
  writeNpyUint32(r/"sequences"/"symbols.npy",s.storageV1SequenceSymbols)
  writeNpyUint64(r/"moments"/"snapshot_id.npy",s.storageV1MomentSnapshotIds); writeNpyUint32(r/"moments"/"population_scope_id.npy",s.storageV1MomentPopulationScopeIds); writeNpyUint32(r/"moments"/"mass_basis_id.npy",s.storageV1MomentMassBasisIds); writeNpyUint64(r/"moments"/"chain_count.npy",s.storageV1MomentChainCounts)
  writeNpyFloat64(r/"moments"/"sum_dp.npy",s.storageV1MomentSumDp); writeNpyFloat64(r/"moments"/"sum_dp2.npy",s.storageV1MomentSumDp2); writeNpyFloat64(r/"moments"/"dp_n.npy",s.storageV1MomentDpN); writeNpyFloat64(r/"moments"/"dp_w.npy",s.storageV1MomentDpW)
  writeNpyFloat64(r/"moments"/"sum_molar_mass.npy",s.storageV1MomentSumMass); writeNpyFloat64(r/"moments"/"sum_molar_mass2.npy",s.storageV1MomentSumMass2); writeNpyFloat64(r/"moments"/"sum_molar_mass3.npy",s.storageV1MomentSumMass3); writeNpyFloat64(r/"moments"/"mn.npy",s.storageV1MomentMn); writeNpyFloat64(r/"moments"/"mw.npy",s.storageV1MomentMw); writeNpyFloat64(r/"moments"/"mz.npy",s.storageV1MomentMz); writeNpyFloat64(r/"moments"/"dispersity.npy",s.storageV1MomentDispersity)

  writeNpyUint64(r/"memory"/"snapshot_id.npy",s.storageV1MemorySnapshotIds); writeNpyUint64(r/"memory"/"live_chains.npy",s.storageV1MemoryLiveChains); writeNpyUint64(r/"memory"/"dead_records.npy",s.storageV1MemoryDeadRecords); writeNpyUint64(r/"memory"/"dead_chains.npy",s.storageV1MemoryDeadChains)
  writeNpyUint64(r/"memory"/"stored_live_mers.npy",s.storageV1MemoryStoredLiveMers); writeNpyUint64(r/"memory"/"stored_dead_mers.npy",s.storageV1MemoryStoredDeadMers); writeNpyUint64(r/"memory"/"live_seq_B.npy",s.storageV1MemoryLiveSeqBytes); writeNpyUint64(r/"memory"/"live_object_B.npy",s.storageV1MemoryLiveObjectBytes); writeNpyUint64(r/"memory"/"dead_record_B.npy",s.storageV1MemoryDeadRecordBytes); writeNpyUint64(r/"memory"/"total_est_B.npy",s.storageV1MemoryTotalEstBytes)
  writeNpyUint64(r/"microstructure_motifs"/"snapshot_id.npy",s.storageV1MotifSnapshotIds)
  writeNpyUint32(r/"microstructure_motifs"/"motif_order.npy",s.storageV1MotifOrders)
  writeNpyUint32(r/"microstructure_motifs"/"motif_id.npy",s.storageV1MotifIds)
  writeNpyUint64(r/"microstructure_motifs"/"count.npy",s.storageV1MotifCounts)
  writeNpyUint64(r/"block_statistics"/"snapshot_id.npy",s.storageV1BlockSnapshotIds)
  writeNpyUint32(r/"block_statistics"/"monomer_id.npy",s.storageV1BlockMonomerIds)
  writeNpyUint64(r/"block_statistics"/"block_length.npy",s.storageV1BlockLengths)
  writeNpyUint64(r/"block_statistics"/"block_count.npy",s.storageV1BlockCounts)

proc approxEqual(a,b: float64; rtol=1e-12): bool =
  if classify(a) in {fcNan, fcInf, fcNegInf} or classify(b) in {fcNan, fcInf, fcNegInf}: return false
  abs(a-b) <= rtol * max(abs(a), abs(b))

proc validateCore(m:Model;s:State): tuple[ok:bool, errors:int, checks:seq[JsonNode]] =
  result.ok=true
  template check(name:string; cond:bool) =
    let passed = cond
    result.checks.add %*{"check":name,"status":(if passed:"pass" else:"fail"),"severity":"error"}
    if not passed: result.ok=false; inc result.errors

  check("snapshots_present", s.storageV1SnapshotIds.len>0)
  var finals=0
  for x in s.storageV1IsFinal:
    if x: inc finals
  check("single_final_snapshot", finals==1)

  let ne=m.monomers.len+m.species.len+2
  check("state_dense_cardinality", s.storageV1StateCounts.len == s.storageV1SnapshotIds.len*ne)
  var stateUnitsOk=true
  for i,c in s.storageV1StateCounts:
    let mol=float64(c)/AvogadroConstantMolInv
    let sid = int(s.storageV1StateSnapshotIds[i])
    if not approxEqual(s.storageV1StateMoles[i],mol) or not approxEqual(s.storageV1StateConcentrations[i],mol/s.storageV1KmcVolumeL[sid]): stateUnitsOk=false
  check("state_count_unit_conversions",stateUnitsOk)

  var channelOk=true
  for i,sid in s.storageV1SnapshotIds:
    var total=0'u64
    for j in 0..<m.channels.len:
      let k=i*m.channels.len+j
      if s.storageV1ChannelEventCounts[k] != s.storageV1ChannelProductiveCounts[k]+s.storageV1ChannelNonproductiveCounts[k]: channelOk=false
      total += s.storageV1ChannelEventCounts[k]
    if total != s.storageV1KmcEvents[i]: channelOk=false
  check("channel_event_identities",channelOk)

  let na=s.storageV1ActionIds.len
  var actionLengthsOk=true
  for n in [s.storageV1ActionKmcEvents.len,s.storageV1ActionTimes.len,s.storageV1ActionSourceLines.len,s.storageV1ActionTypeIds.len,s.storageV1ActionTriggerTypeIds.len,s.storageV1ActionScheduledTimes.len,s.storageV1ActionTargetIds.len,s.storageV1ActionRequestedValues.len,s.storageV1ActionBeforeValues.len,s.storageV1ActionAfterValues.len,s.storageV1ActionStateChanged.len,s.storageV1ActionOutputWritten.len,s.storageV1ActionHasSnapshot.len,s.storageV1ActionSnapshotIds.len,s.storageV1ActionHasKineticSet.len,s.storageV1ActionKineticSetIds.len]:
    if n != na: actionLengthsOk=false
  check("actions_equal_column_lengths",actionLengthsOk)
  var actionIdsOk=true
  for i,id in s.storageV1ActionIds:
    if id != uint64(i): actionIdsOk=false
  check("action_ids_continuous",actionIdsOk)

  var setIdsOk=true
  for i,id in s.storageV1KineticSetIds:
    if id != uint64(i): setIdsOk=false
  check("kinetic_set_ids_continuous",setIdsOk)
  var np=1 + m.rates.len
  for r in m.rates:
    if r.declaredArrhenius: np += 2
  check("kinetic_parameter_sets_complete",s.storageV1KineticValues.len == s.storageV1KineticSetIds.len*np)
  check("action_conditions_equal_lengths",s.storageV1ConditionRecordIds.len == s.storageV1ConditionActionIds.len and s.storageV1ConditionRecordIds.len == s.storageV1ConditionIndexes.len)

  let nr=s.storageV1ChainRecordIds.len
  var chainLengthsOk=true
  for n in [s.storageV1ChainSnapshotIds.len,s.storageV1ChainPopulationIds.len,s.storageV1ChainPoolIds.len,s.storageV1ChainOriginIds.len,s.storageV1ChainDp.len,s.storageV1ChainMolarMass.len,s.storageV1ChainCounts.len,s.storageV1ChainMoles.len,s.storageV1ChainConcentrations.len,s.storageV1ChainLeftEndIds.len,s.storageV1ChainRightEndIds.len,s.storageV1ChainHasFirst.len,s.storageV1ChainFirstIds.len,s.storageV1ChainHasPenultimate.len,s.storageV1ChainPenultimateIds.len,s.storageV1ChainHasLast.len,s.storageV1ChainLastIds.len,s.storageV1ChainHasSequence.len,s.storageV1ChainSequenceOffsets.len,s.storageV1ChainSequenceLengths.len]:
    if n != nr: chainLengthsOk=false
  check("chains_equal_column_lengths",chainLengthsOk)

  var chainIdsOk=true; var chainUnitsOk=true; var compOk=true; var terminalsOk=true; var seqOk=true; var massOk=true
  check("chain_composition_dense_cardinality",s.storageV1CompositionUnitCounts.len == nr*m.monomers.len)
  for i in 0..<nr:
    if s.storageV1ChainRecordIds[i] != uint64(i): chainIdsOk=false
    let c=s.storageV1ChainCounts[i]; let mol=float64(c)/AvogadroConstantMolInv
    if c==0 or s.storageV1ChainDp[i]==0: chainUnitsOk=false
    let sid = int(s.storageV1ChainSnapshotIds[i])
    if sid < 0 or sid >= s.storageV1KmcVolumeL.len:
      chainUnitsOk=false
    else:
      let snapshotKmcVolumeL = s.storageV1KmcVolumeL[sid]
      if snapshotKmcVolumeL <= 0.0 or
         not approxEqual(s.storageV1ChainMoles[i],mol) or
         not approxEqual(s.storageV1ChainConcentrations[i],mol/snapshotKmcVolumeL):
        chainUnitsOk=false
    var sumdp=0'u64; var repeatMass=0.0
    let base=i*m.monomers.len
    if base+m.monomers.len <= s.storageV1CompositionUnitCounts.len:
      for j,mo in m.monomers:
        if s.storageV1CompositionRecordIds[base+j] != uint64(i) or s.storageV1CompositionMonomerIds[base+j] != uint32(j): compOk=false
        let u=s.storageV1CompositionUnitCounts[base+j]; sumdp += u; repeatMass += float64(u)*mo.mw
    else: compOk=false
    if sumdp != s.storageV1ChainDp[i]: compOk=false
    if s.storageV1ChainHasFirst[i] != (s.storageV1ChainDp[i]>=1): terminalsOk=false
    if s.storageV1ChainHasLast[i] != (s.storageV1ChainDp[i]>=1): terminalsOk=false
    if s.storageV1ChainHasPenultimate[i] != (s.storageV1ChainDp[i]>=2): terminalsOk=false
    if s.storageV1ChainHasFirst[i] and s.storageV1ChainFirstIds[i] >= uint32(m.monomers.len): terminalsOk=false
    if s.storageV1ChainHasPenultimate[i] and s.storageV1ChainPenultimateIds[i] >= uint32(m.monomers.len): terminalsOk=false
    if s.storageV1ChainHasLast[i] and s.storageV1ChainLastIds[i] >= uint32(m.monomers.len): terminalsOk=false
    let hs=s.storageV1ChainHasSequence[i]; let off=s.storageV1ChainSequenceOffsets[i]; let ln=s.storageV1ChainSequenceLengths[i]
    if hs:
      if ln != s.storageV1ChainDp[i] or off+ln > uint64(s.storageV1SequenceSymbols.len): seqOk=false
      else:
        var seqCounts=newSeq[uint64](m.monomers.len)
        for k in off..<off+ln:
          let sym=s.storageV1SequenceSymbols[int(k)]
          if sym>=uint32(m.monomers.len): seqOk=false
          else: inc seqCounts[int(sym)]
        for j in 0..<m.monomers.len:
          if seqCounts[j] != s.storageV1CompositionUnitCounts[base+j]: seqOk=false
        if ln>0 and s.storageV1SequenceSymbols[int(off)] != s.storageV1ChainFirstIds[i]: seqOk=false
        if ln>0 and s.storageV1SequenceSymbols[int(off+ln-1)] != s.storageV1ChainLastIds[i]: seqOk=false
        if ln>1 and s.storageV1SequenceSymbols[int(off+ln-2)] != s.storageV1ChainPenultimateIds[i]: seqOk=false
    else:
      if off!=0 or ln!=0: seqOk=false
    var expectedMass=repeatMass
    if m.mass_model == mmWithEndgroups:
      let le=s.storageV1ChainLeftEndIds[i]; let re=s.storageV1ChainRightEndIds[i]
      for eid in [le,re]:
        if eid > 0 and eid <= uint32(s.storageV1ObservedEndGroups.len):
          let name = s.storageV1ObservedEndGroups[int(eid)-1]
          if m.endgroupByName.hasKey(name):
            expectedMass += m.endgroups[m.endgroupByName[name]].mw
          # Undeclared mechanism end groups intentionally contribute 0.0,
          # matching the engine's compatibility semantics.
    if not approxEqual(s.storageV1ChainMolarMass[i],expectedMass): massOk=false
  check("chain_record_ids_continuous",chainIdsOk)
  check("chain_count_unit_conversions",chainUnitsOk)
  check("chain_composition_reconstructs_dp",compOk)
  check("chain_terminal_fields_consistent",terminalsOk)
  check("chain_sequences_consistent",seqOk)
  check("chain_molar_mass_reconstruction",massOk)

  var momentsOk = s.storageV1MomentSnapshotIds.len mod 6 == 0
  var mi=0
  while momentsOk and mi<s.storageV1MomentSnapshotIds.len:
    let sid=s.storageV1MomentSnapshotIds[mi]; let scope=s.storageV1MomentPopulationScopeIds[mi]; let basis=s.storageV1MomentMassBasisIds[mi]
    var n=0'u64; var sd=0.0; var sd2=0.0; var sm=0.0; var sm2=0.0; var sm3=0.0
    for i in 0..<nr:
      if s.storageV1ChainSnapshotIds[i]!=sid: continue
      let pop=s.storageV1ChainPopulationIds[i]
      if scope==1 and pop!=0: continue
      if scope==2 and pop!=1: continue
      let c=s.storageV1ChainCounts[i]; let d=float64(s.storageV1ChainDp[i]); var mass=s.storageV1ChainMolarMass[i]
      if basis==0:
        mass=0
        let base=i*m.monomers.len
        for j,mo in m.monomers: mass += float64(s.storageV1CompositionUnitCounts[base+j])*mo.mw
      n+=c; sd+=float64(c)*d; sd2+=float64(c)*d*d; sm+=float64(c)*mass; sm2+=float64(c)*mass*mass; sm3+=float64(c)*mass*mass*mass
    if s.storageV1MomentChainCounts[mi]!=n or not approxEqual(s.storageV1MomentSumDp[mi],sd) or not approxEqual(s.storageV1MomentSumDp2[mi],sd2) or not approxEqual(s.storageV1MomentSumMass[mi],sm) or not approxEqual(s.storageV1MomentSumMass2[mi],sm2) or not approxEqual(s.storageV1MomentSumMass3[mi],sm3): momentsOk=false
    inc mi
  check("moments_reconstruct_from_chains",momentsOk)

  var memOk = s.storageV1MemorySnapshotIds.len == s.storageV1MomentSnapshotIds.len div 6
  check("memory_rows_match_chain_snapshots",memOk)

proc finalizeStorageV1*(m:Model;s:var State;wallTime:float64;terminationReason="") =
  writeJsonlAtomic(s.storageV1RunDir / "schema.jsonl", schemaRecords(m, s.storageV1ObservedEndGroups))
  writeColumns(s)
  let v=validateCore(m,s)
  writeJsonlAtomic(s.storageV1RunDir/"diagnostics"/"validation.jsonl",v.checks)
  let status=if v.ok:"completed" else:"failed"
  let storageManifest = writeStorageManifest(s.storageV1RunDir)
  writePrettyJsonAtomic(s.storageV1RunDir / "run_metadata.json", metadataNode(m, s, status, utcNow(), wallTime, (if v.ok: 0 else: 1), terminationReason, (if v.ok: "passed" else: "failed"), v.errors, storageManifest))
  if v.ok:
    if dirExists(s.storageV1RunDir/".work"): removeDir(s.storageV1RunDir/".work")
    writeFile(s.storageV1RunDir/"RESULTS_COMPLETE","slimmc-storage-v1\n")


proc finalizeInterruptedStorageV1*(m: Model; s: var State; wallTime: float64; exitCode: int = 130) =
  ## Publish all currently captured canonical columns without marking the run complete.
  writeJsonlAtomic(s.storageV1RunDir / "schema.jsonl", schemaRecords(m, s.storageV1ObservedEndGroups))
  writeColumns(s)
  var lines: seq[JsonNode] = @[]
  lines.add %*{"check":"copo_core_storage_v1","status":"not_completed","severity":"info"}
  writeJsonlAtomic(s.storageV1RunDir/"diagnostics"/"validation.jsonl", lines)
  let storageManifest = writeStorageManifest(s.storageV1RunDir)
  writePrettyJsonAtomic(s.storageV1RunDir / "run_metadata.json",
    metadataNode(m, s, "interrupted", utcNow(), wallTime, exitCode, "user_interrupt", "not_completed", 0, storageManifest))
  if fileExists(s.storageV1RunDir / "RESULTS_COMPLETE"):
    removeFile(s.storageV1RunDir / "RESULTS_COMPLETE")

proc markFailedStorageV1*(m: Model; message: string; exitCode: int = 1) =
  ## Best-effort metadata repair for exceptions escaping the copo engine.
  let root = m.output_dir
  let path = root / "run_metadata.json"
  if not fileExists(path): return
  try:
    var node = parseJson(readFile(path))
    node["run_status"] = %"failed"
    node["exit_code"] = %exitCode
    node["finished_at_utc"] = %utcNow()
    node["validation_status"] = %"not_completed"
    node["failure_message"] = %message
    writePrettyJsonAtomic(path, node)
    if fileExists(root / "RESULTS_COMPLETE"): removeFile(root / "RESULTS_COMPLETE")
  except CatchableError:
    discard
