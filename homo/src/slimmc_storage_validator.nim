import std/[json, math, os, sets, sequtils]
import slimmc_types
import ../../common/results_types

type ValidationSummary* = object
  records*: seq[JsonNode]
  warningCount*: int
  errorCount*: int

proc addCheck(v: var ValidationSummary; name: string; ok: bool;
              severity = "error"; details = "") =
  var rec = newJObject()
  rec["check"] = %name
  rec["status"] = %(if ok: "pass" else: "fail")
  rec["severity"] = %severity
  if details.len > 0: rec["details"] = %details
  v.records.add rec
  if not ok:
    if severity == "warning": inc v.warningCount else: inc v.errorCount

proc approxEq(a, b: float64; rtol = 1e-12): bool =
  if classify(a) in {fcNan, fcInf, fcNegInf} or classify(b) in {fcNan, fcInf, fcNegInf}:
    return false
  abs(a - b) <= rtol * max(abs(a), abs(b))

proc continuous(ids: openArray[uint64]): bool =
  for i, x in ids:
    if x != uint64(i): return false
  true

proc validateStorageV1*(m: Model; s: State; status: RunStatus): ValidationSummary =
  let root = m.outputDir
  let requiredFiles = [
    "input.model", "schema.jsonl", "run_metadata.json",
    "snapshots/snapshot_id.npy", "snapshots/time.npy", "snapshots/kmc_event.npy",
    "state/snapshot_id.npy", "state/entity_id.npy", "state/count.npy",
    "chains/chain_record_id.npy", "moments/snapshot_id.npy",
    "channel_events/event_count.npy", "actions/action_id.npy",
    "action_conditions/condition_record_id.npy",
    "kinetic_parameters/sets/kinetic_parameter_set_id.npy",
    "kinetic_parameters/values/kinetic_parameter_id.npy",
    "diagnostics/run.log"
  ]
  var filesOk = true
  for p in requiredFiles:
    if not fileExists(root / p): filesOk = false
  result.addCheck("required_files_present", filesOk)

  let ns = s.storageV1SnapshotIds.len
  let snapshotColsOk = [s.storageV1Times.len, s.storageV1KmcEvents.len,
    s.storageV1ReasonIds.len, s.storageV1IsFinal.len, s.storageV1HasChains.len,
    s.storageV1HasSequences.len, s.storageV1ParameterSetIds.len].allIt(it == ns)
  result.addCheck("snapshot_column_lengths", snapshotColsOk)
  result.addCheck("snapshot_ids_continuous", continuous(s.storageV1SnapshotIds))
  var orderOk = true
  for i in 1 ..< ns:
    if s.storageV1Times[i] < s.storageV1Times[i-1] or
       s.storageV1KmcEvents[i] < s.storageV1KmcEvents[i-1]: orderOk = false
  result.addCheck("snapshot_order_nondecreasing", orderOk)
  var finals = 0
  for x in s.storageV1IsFinal:
    if x: inc finals
  result.addCheck("snapshot_finality", if status == rsCompleted: finals == 1 else: finals == 0,
    details = "final_count=" & $finals)

  let nEntities = m.species.len + m.pools.len + 2
  let expectedStateRows = ns * nEntities
  let stateColsOk = [s.storageV1StateSnapshotIds.len, s.storageV1StateEntityIds.len,
    s.storageV1StateCounts.len, s.storageV1StateMoles.len,
    s.storageV1StateConcentrations.len].allIt(it == expectedStateRows)
  result.addCheck("state_dense_shape", stateColsOk,
    details = "rows=" & $s.storageV1StateCounts.len & ", expected=" & $expectedStateRows)
  var stateOrderOk = stateColsOk
  var stateUnitsOk = stateColsOk
  if stateColsOk:
    for row in 0 ..< expectedStateRows:
      let sid = row div nEntities
      let eid = row mod nEntities
      if s.storageV1StateSnapshotIds[row] != uint64(sid) or
         s.storageV1StateEntityIds[row] != uint32(eid): stateOrderOk = false
      let mol = float64(s.storageV1StateCounts[row]) / AvogadroConstantMolInv
      let conc = mol / s.storageV1KmcVolumeL[sid]
      if not approxEq(s.storageV1StateMoles[row], mol) or
         not approxEq(s.storageV1StateConcentrations[row], conc): stateUnitsOk = false
  result.addCheck("state_dense_order", stateOrderOk)
  result.addCheck("state_count_unit_conversions", stateUnitsOk)

  let nc = m.channels.len
  let channelRows = ns * nc
  let channelColsOk = [s.storageV1ChannelSnapshotIds.len, s.storageV1ChannelIds.len,
    s.storageV1ChannelEventCounts.len, s.storageV1ChannelProductiveCounts.len,
    s.storageV1ChannelNonproductiveCounts.len].allIt(it == channelRows)
  result.addCheck("channel_event_dense_shape", channelColsOk)
  var channelIdentityOk = channelColsOk
  var channelSnapshotSumOk = channelColsOk
  if channelColsOk:
    for sid in 0 ..< ns:
      var total = 0'u64
      for cid in 0 ..< nc:
        let row = sid * nc + cid
        if s.storageV1ChannelSnapshotIds[row] != uint64(sid) or
           s.storageV1ChannelIds[row] != uint32(cid): channelIdentityOk = false
        if s.storageV1ChannelEventCounts[row] !=
           s.storageV1ChannelProductiveCounts[row] + s.storageV1ChannelNonproductiveCounts[row]:
          channelIdentityOk = false
        total += s.storageV1ChannelEventCounts[row]
      if total != s.storageV1KmcEvents[sid]: channelSnapshotSumOk = false
  result.addCheck("channel_event_identities", channelIdentityOk)
  result.addCheck("kmc_event_equals_channel_sum", channelSnapshotSumOk)

  result.addCheck("kinetic_set_ids_continuous", continuous(s.storageV1KineticSetIds))
  let nk = if s.storageV1KineticSetIds.len == 0: 0 else:
    s.storageV1KineticValues.len div s.storageV1KineticSetIds.len
  var kineticDenseOk = s.storageV1KineticSetIds.len > 0 and
    s.storageV1KineticValues.len == s.storageV1KineticSetIds.len * nk and
    s.storageV1KineticValueSetIds.len == s.storageV1KineticValues.len and
    s.storageV1KineticParameterIds.len == s.storageV1KineticValues.len
  if kineticDenseOk:
    for row in 0 ..< s.storageV1KineticValues.len:
      if s.storageV1KineticValueSetIds[row] != uint64(row div nk) or
         s.storageV1KineticParameterIds[row] != uint32(row mod nk) or
         classify(s.storageV1KineticValues[row]) in {fcNan, fcInf, fcNegInf}:
        kineticDenseOk = false
  result.addCheck("kinetic_parameter_sets_complete", kineticDenseOk)

  result.addCheck("action_ids_continuous", continuous(s.storageV1ActionIds))
  var actionMasksOk = true
  for i in 0 ..< s.storageV1ActionIds.len:
    if not s.storageV1ActionHasSnapshot[i] and s.storageV1ActionSnapshotIds[i] != 0: actionMasksOk = false
    if not s.storageV1ActionHasKineticSet[i] and s.storageV1ActionKineticSetIds[i] != 0: actionMasksOk = false
  result.addCheck("action_optional_masks", actionMasksOk)
  result.addCheck("condition_record_ids_continuous", continuous(s.storageV1ConditionRecordIds))
  var conditionsOk = true
  var lastAction = 0'u64
  var nextIndex = 0'u32
  for i in 0 ..< s.storageV1ConditionRecordIds.len:
    let aid = s.storageV1ConditionActionIds[i]
    if aid >= uint64(s.storageV1ActionIds.len): conditionsOk = false
    if i == 0 or aid != lastAction:
      nextIndex = 0
      lastAction = aid
    if s.storageV1ConditionIndexes[i] != nextIndex or not s.storageV1ConditionMet[i]: conditionsOk = false
    inc nextIndex
  result.addCheck("action_conditions_and_order", conditionsOk)

  result.addCheck("chain_record_ids_continuous", continuous(s.storageV1ChainRecordIds))
  var chainsOk = true
  var seen = initHashSet[string]()
  for i in 0 ..< s.storageV1ChainRecordIds.len:
    if s.storageV1ChainDp[i] < 1 or s.storageV1ChainCounts[i] < 1 or
       s.storageV1ChainMolarMass[i] <= 0 or
       not approxEq(s.storageV1ChainMoles[i], float64(s.storageV1ChainCounts[i]) / AvogadroConstantMolInv) or
       not approxEq(s.storageV1ChainConcentrations[i], s.storageV1ChainMoles[i] / s.storageV1KmcVolumeL[int(s.storageV1ChainSnapshotIds[i])]) or
       s.storageV1ChainSequenceLengths[i] != 0:
      chainsOk = false
    let key = $s.storageV1ChainSnapshotIds[i] & ":" & $s.storageV1ChainPopulationIds[i] & ":" &
      $s.storageV1ChainPoolIds[i] & ":" & $s.storageV1ChainOriginIds[i] & ":" &
      $s.storageV1ChainDp[i] & ":" & $s.storageV1ChainLeftEndIds[i] & ":" & $s.storageV1ChainRightEndIds[i]
    if key in seen: chainsOk = false else: seen.incl key
  result.addCheck("chains_valid_and_merged", chainsOk)

  let momentsColsOk = [s.storageV1MomentPopulationScopeIds.len, s.storageV1MomentMassBasisIds.len,
    s.storageV1MomentChainCounts.len, s.storageV1MomentSumDp.len, s.storageV1MomentSumDp2.len,
    s.storageV1MomentDpN.len, s.storageV1MomentDpW.len, s.storageV1MomentSumMass.len,
    s.storageV1MomentSumMass2.len, s.storageV1MomentSumMass3.len, s.storageV1MomentMn.len,
    s.storageV1MomentMw.len, s.storageV1MomentMz.len, s.storageV1MomentDispersity.len].allIt(
      it == s.storageV1MomentSnapshotIds.len)
  result.addCheck("moments_column_lengths", momentsColsOk)
  var momentsShapeOk = momentsColsOk
  if momentsColsOk:
    var snapshotsWithChains = 0
    for x in s.storageV1HasChains:
      if x: inc snapshotsWithChains
    momentsShapeOk = s.storageV1MomentSnapshotIds.len == snapshotsWithChains * 6
  result.addCheck("moments_six_rows_per_chain_snapshot", momentsShapeOk)
