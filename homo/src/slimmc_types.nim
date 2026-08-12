# slimmc_types.nim
# Shared types, constants, and small utilities for slimmc.

import strutils, math, tables, hashes
import ../../common/model_contract
import ../../common/safe_numeric

const
  AppName* = "slimmc"
  SlimmcVersion* {.strdefine.} = "5.0.0"
  AppVersion* = SlimmcVersion

  NA* = 6.02214076e23
  Rgas* = 8.31446261815324
  InfTime* = 1.0e300
  Eps* = 1.0e-12
  TimeRelEps* = 8.0e-15

proc timeTolerance*(a, b: float): float =
  ## Mixed absolute/relative tolerance for simulation-time comparisons.
  ## Eps preserves sub-second behavior; the relative term stays several
  ## floating-point ulps wide when t becomes very large.
  max(Eps, TimeRelEps * max(abs(a), abs(b)))

proc timeClose*(a, b: float): bool =
  abs(a - b) <= timeTolerance(a, b)

type
  SpeciesKind* = enum
    skSpecies, skMonomer

  PoolKind* = enum
    pkActive, pkDead

  FormationKind* = enum
    ## Provenance of the event that created the current chain record.
    ## Dead-chain values are written to chains.csv as formed_by.
    fbUnknown,
    fbInit,
    fbTransferM,
    fbTermC,
    fbTermD_H,
    fbTermD_U,
    fbTermX,
    fbTransferH

  RateKind* = enum
    rkFixed, rkArr

  ChannelKind* = enum
    chRxnUni,
    chRxnBiDiff,
    chRxnBiSame,
    chMacroInit,
    chMacroProp,
    chMacroDeprop,
    chMacroTermC,
    chMacroTermD,
    chMacroTermX,
    chMacroTransferH,
    chMacroTransferM

  ActionKind* = enum
    eaPrint,
    eaPrintInfo,
    eaSave,
    eaSaveChains,
    eaStop,
    eaSetK,
    eaAddK,
    eaSetTemp,
    eaAddTemp,
    eaSetC,
    eaAddC,
    eaFeed,
    eaPrintMemory

  ConditionalObservableKind* = enum
    woConversion,
    woSpeciesConc

  ComparisonKind* = enum
    coGreater,
    coLess

  SnapshotRecord* = object
    snapshotId*: int64
    kmcEvent*: int64
    t*: float
    level*: string
    isFinal*: bool
    reason*: string
    parameterStateId*: int64
    stateRevision*: int64

  SpeciesDef* = object
    name*: string
    kind*: SpeciesKind
    c0*: float
    mw*: float
    hasMw*: bool

  PoolDef* = object
    name*: string
    kind*: PoolKind

  VarDef* = object
    kind*: string
    name*: string
    value*: float
    unit*: string

  RateDef* = object
    name*: string
    kind*: RateKind
    kConst*: float
    Apre*: float
    Ea*: float
    declaredArrhenius*: bool

  RawLine* = object
    text*: string
    lineNo*: int

  RawTerm* = object
    name*: string
    stoich*: int

  Term* = object
    sp*: int
    stoich*: int

  Chain* = object
    eg1*: int
    dp*: int
    eg2*: int
    formedBy*: FormationKind

  StructKey* = object
    eg1*: int
    dp*: int
    eg2*: int
    formedBy*: FormationKind

  KmcChannel* = object
    name*: string
    expr*: string
    lineNo*: int
    kind*: ChannelKind
    kId*: int
    eff*: float

    lhs*: seq[Term]
    rhs*: seq[Term]

    sp1*: int
    sp2*: int
    pool1*: int
    pool2*: int
    poolOut*: int

  FeedDef* = object
    name*: string
    concentrations*: seq[float]

  ScheduledAction* = object
    startTime*: float
    nextTime*: float
    period*: float
    repeat*: bool
    remaining*: int64  # -1 unlimited, otherwise executions left
    active*: bool
    action*: ActionKind
    args*: seq[string]
    lineNo*: int

  AtomicCondition* = object
    observable*: ConditionalObservableKind
    speciesId*: int
    comparison*: ComparisonKind
    threshold*: float

  ConditionalAction* = object
    active*: bool
    conditions*: seq[AtomicCondition]
    # Convenience fields mirror the first condition for direct access.
    observable*: ConditionalObservableKind
    speciesId*: int
    comparison*: ComparisonKind
    threshold*: float
    action*: ActionKind
    args*: seq[string]
    lineNo*: int

  MemoryPolicy* = object
    ## Parity port from slimmc-copo v0.6: an opt-in runtime memory budget.
    ## `compact_dead`/`drop_dead_seq` are accepted as directive arguments for
    ## syntax parity but are no-ops here too -- slimmc's Chain is already a
    ## fixed-size record with no per-mer sequence to compact, matching
    ## copo's own "v0.5 always stores DEAD summaries" no-op for that token.
    limitBytes*: int64
    hasLimit*: bool
    snapshotOnLimit*: bool
    stopOnLimit*: bool

  MemoryEstimate* = object
    ## Conservative, transparent estimate of chain-storage memory. Not OS
    ## RSS. Simpler than slimmc-copo's equivalent because slimmc's Chain is
    ## a compact fixed-size record (eg1, dp, eg2, formedBy) with no stored
    ## per-mer sequence -- there is nothing analogous to copo's
    ## storedLiveMers/storedDeadMers to report.
    liveChains*: int64
    deadChains*: int64
    liveObjectBytes*: int64
    deadRecordBytes*: int64
    totalBytes*: int64

  Model* = object
    V*: float
    initVolumeMl*: float
    currentVolumeMl*: float
    hasInitVolume*: bool
    T*: float
    tEnd*: float
    maxEvents*: int64
    whenCheckEvents*: int64
    seed*: int64
    memoryPolicy*: MemoryPolicy

    # Source/model metadata and output paths for the v2.8 data contract.
    modelSourceFile*: string
    runId*: string
    outputDir*: string

    description*: string
    hasDescription*: bool
    variables*: seq[VarDef]

    massModel*: string
    dpMax*: int64
    sequenceMode*: string

    species*: seq[SpeciesDef]
    speciesByName*: Table[string, int]
    monomerId*: int

    feeds*: seq[FeedDef]
    feedByName*: Table[string, int]

    pools*: seq[PoolDef]
    poolByName*: Table[string, int]

    rates*: seq[RateDef]
    rateByName*: Table[string, int]

    egNames*: seq[string]
    egByName*: Table[string, int]
    egMass*: seq[float]
    egMassKnown*: seq[bool]
    egMassSource*: seq[string]
    egH*: int
    egU*: int
    egActive*: int

    channels*: seq[KmcChannel]
    scheduledActions*: seq[ScheduledAction]
    conditionalActions*: seq[ConditionalAction]
    rawReactions*: seq[RawLine]
    rawConditionalActions*: seq[RawLine]

  State* = object
    t*: float
    kmcEvent*: int64
    actionNo*: int64
    scheduledActionNo*: int64
    conditionalActionNo*: int64
    stateRevision*: int64
    parameterStateId*: int64
    startedAt*: string
    feedRemainders*: seq[seq[float]]
    speciesExternalN*: seq[int64]
    speciesDosedN*: seq[int64]
    # snapshotId identifies one logical snapshot of simulation state and is
    # once per logical snapshot, in saveSnapshot() (see slimmc_kmc.nim) --
    # never inside an individual file writer. `save` produces a light
    # the same exact state revision are coalesced, and a light snapshot can
    # be upgraded to full without duplicating Storage snapshot rows.
    snapshotId*: int64
    # Monotonic identifier of a compressed Storage chain record.
    chainRecordId*: int64
    snapshots*: seq[SnapshotRecord]
    # Runtime flags used by snapshot/action bookkeeping.
    anySnapshotWritten*: bool
    chainsWritten*: bool
    channelTraceRowsWritten*: int64  # item 67: trace_channels_rows
    channelTraceTruncated*: bool     # item 67: trace_channels_truncated (limit reached, sim continued)
    n*: seq[int64]
    pools*: seq[seq[Chain]]
    # Incremental eligibility counters keep SSA propensities O(1) in the
    # number of stored chains. They are runtime-only and are not persisted.
    poolPropagatableCounts*: seq[int64]
    poolDepropableCounts*: seq[int64]
    poolEligibilityTrackedLen*: seq[int64]
    channelFires*: seq[int64]
    channelSuccesses*: seq[int64]  # item 82: for f<1 rxn channels, fires = successes + failures (exact invariant)
    channelFailures*: seq[int64]
    sumA0Tau*: float    # item 86: sum of a0*tau across every SSA draw (should average to 1.0 under Exp(1))
    sumA0TauSq*: float  # item 86: sum of (a0*tau)^2, for the sample variance (should be ~1.0 under Exp(1))
    countA0Tau*: int64  # item 86: number of draws contributing to the above
    mExpected*: int64
    mBalance*: int64  # item 24: independent ledger for the single monomer, updated only by elementary rxns and set_c/add_c on the monomer species -- NOT by init/prop/deprop/transfer_m, which only move mass between free and incorporated pools (mass-conserving by construction). Should equal M_count + poly_M_count if bookkeeping is consistent.
    observedEg*: seq[bool]
    stopRequested*: bool
    stopLineNo*: int
    stopCheckSource*: string
    stopConditions*: seq[AtomicCondition]
    stopActualValues*: seq[float]

    # Slimmc Storage v1 stage C buffers (0-based snapshot ids, dense state).
    storageV1SnapshotIds*: seq[uint64]
    storageV1StateRevisions*: seq[int64]
    storageV1Times*: seq[float64]
    storageV1KmcEvents*: seq[uint64]
    storageV1ReasonIds*: seq[uint32]
    storageV1IsFinal*: seq[bool]
    storageV1HasChains*: seq[bool]
    storageV1HasSequences*: seq[bool]
    storageV1ParameterSetIds*: seq[uint64]
    storageV1VolumeMl*: seq[float64]
    storageV1KmcVolumeL*: seq[float64]
    storageV1ChainCountLive*: seq[uint64]
    storageV1ChainCountDead*: seq[uint64]
    storageV1ChainCountTotal*: seq[uint64]
    storageV1StateSnapshotIds*: seq[uint64]
    storageV1StateEntityIds*: seq[uint32]
    storageV1StateCounts*: seq[uint64]
    storageV1StateMoles*: seq[float64]
    storageV1StateConcentrations*: seq[float64]
    # Slimmc Storage v1 cumulative channel-event rows, dense by snapshot × channel.
    storageV1ChannelSnapshotIds*: seq[uint64]
    storageV1ChannelIds*: seq[uint32]
    storageV1ChannelEventCounts*: seq[uint64]
    storageV1ChannelProductiveCounts*: seq[uint64]
    storageV1ChannelNonproductiveCounts*: seq[uint64]

    # Slimmc Storage v1 compressed chain records.
    storageV1ChainRecordIds*: seq[uint64]
    storageV1ChainSnapshotIds*: seq[uint64]
    storageV1ChainPopulationIds*: seq[uint32]
    storageV1ChainPoolIds*: seq[uint32]
    storageV1ChainOriginIds*: seq[uint32]
    storageV1ChainDp*: seq[uint64]
    storageV1ChainMolarMass*: seq[float64]
    storageV1ChainCounts*: seq[uint64]
    storageV1ChainMoles*: seq[float64]
    storageV1ChainConcentrations*: seq[float64]
    storageV1ChainLeftEndIds*: seq[uint32]
    storageV1ChainRightEndIds*: seq[uint32]
    storageV1ChainSequenceOffsets*: seq[uint64]
    storageV1ChainSequenceLengths*: seq[uint64]

    # Slimmc Storage v1 derived moments, dense by chain snapshot × scope × mass basis.
    storageV1MomentSnapshotIds*: seq[uint64]
    storageV1MomentPopulationScopeIds*: seq[uint32]
    storageV1MomentMassBasisIds*: seq[uint32]
    storageV1MomentChainCounts*: seq[uint64]
    storageV1MomentSumDp*: seq[float64]
    storageV1MomentSumDp2*: seq[float64]
    storageV1MomentDpN*: seq[float64]
    storageV1MomentDpW*: seq[float64]
    storageV1MomentSumMass*: seq[float64]
    storageV1MomentSumMass2*: seq[float64]
    storageV1MomentSumMass3*: seq[float64]
    storageV1MomentMn*: seq[float64]
    storageV1MomentMw*: seq[float64]
    storageV1MomentMz*: seq[float64]
    storageV1MomentDispersity*: seq[float64]

    # Slimmc Storage v1 complete kinetic parameter sets.
    storageV1KineticSetIds*: seq[uint64]
    storageV1KineticStartEvents*: seq[uint64]
    storageV1KineticStartTimes*: seq[float64]
    storageV1KineticHasSourceAction*: seq[bool]
    storageV1KineticSourceActionIds*: seq[uint64]
    storageV1KineticValueSetIds*: seq[uint64]
    storageV1KineticParameterIds*: seq[uint32]
    storageV1KineticValues*: seq[float64]

    # Slimmc Storage v1 executed-action history.
    storageV1ActionIds*: seq[uint64]
    storageV1ActionKmcEvents*: seq[uint64]
    storageV1ActionTimes*: seq[float64]
    storageV1ActionSourceLines*: seq[uint32]
    storageV1ActionTypeIds*: seq[uint32]
    storageV1ActionTriggerTypeIds*: seq[uint32]
    storageV1ActionScheduledTimes*: seq[float64]
    storageV1ConditionRecordIds*: seq[uint64]
    storageV1ConditionActionIds*: seq[uint64]
    storageV1ConditionIndexes*: seq[uint32]
    storageV1ConditionObservableIds*: seq[uint32]
    storageV1ConditionOperatorIds*: seq[uint32]
    storageV1ConditionThresholds*: seq[float64]
    storageV1ConditionObservedValues*: seq[float64]
    storageV1ConditionMet*: seq[bool]
    storageV1ActionTargetIds*: seq[uint32]
    storageV1ActionRequestedValues*: seq[float64]
    storageV1ActionBeforeValues*: seq[float64]
    storageV1ActionAfterValues*: seq[float64]
    storageV1ActionStateChanged*: seq[bool]
    storageV1ActionOutputWritten*: seq[bool]
    storageV1ActionHasSnapshot*: seq[bool]
    storageV1ActionSnapshotIds*: seq[uint64]
    storageV1ActionHasKineticSet*: seq[bool]
    storageV1ActionKineticSetIds*: seq[uint64]
    storageV1ActionMessages*: seq[string]
    # Semibatch process history, one row per executed feed action.
    storageV1FeedActionIds*: seq[uint64]
    storageV1FeedIds*: seq[uint32]
    storageV1FeedDoseMl*: seq[float64]
    storageV1FeedVolumeBeforeMl*: seq[float64]
    storageV1FeedVolumeAfterMl*: seq[float64]
    storageV1FeedKmcVolumeBeforeL*: seq[float64]
    storageV1FeedKmcVolumeAfterL*: seq[float64]
    # Snapshot-dense monomer ledger.
    storageV1MonomerBalanceSnapshotIds*: seq[uint64]
    storageV1MonomerBalanceMonomerIds*: seq[uint32]
    storageV1MonomerInitialMoles*: seq[float64]
    storageV1MonomerIntroducedMoles*: seq[float64]
    storageV1MonomerFreeMoles*: seq[float64]
    storageV1MonomerIncorporatedMoles*: seq[float64]
    storageV1MonomerConversion*: seq[float64]
    # Snapshot-dense external material balance for all free species/monomers.
    storageV1SpeciesBalanceSnapshotIds*: seq[uint64]
    storageV1SpeciesBalanceEntityIds*: seq[uint32]
    storageV1SpeciesInitialMoles*: seq[float64]
    storageV1SpeciesDosedMoles*: seq[float64]
    storageV1SpeciesTotalMoles*: seq[float64]
    storageV1SpeciesFreeMoles*: seq[float64]
    storageV1SpeciesConsumedMoles*: seq[float64]

    storageTraceKmcEvents*: seq[uint64]
    storageTraceTimes*: seq[float64]
    storageTraceDt*: seq[float64]
    storageTraceChannelIds*: seq[uint32]
    storageTraceRates*: seq[float64]
    storageTracePropensities*: seq[float64]
    storageTraceTotalPropensities*: seq[float64]

  RunOptions* = object
    traceChannelsLimit*: int64
    debug*: bool

  MassMoments* = object
    nChains*: int64
    mn*: float
    mw*: float
    mz*: float
    pdi*: float

  DeadMoments* = object
    repeat*: MassMoments
    endgroups*: MassMoments
    endgroupsComplete*: bool
    dpn*: float
    dpw*: float

  MomentStats* = object
    ## Raw sums plus derived moments.  The raw fields make every value in
    ## Raw sums make every reported moment independently reconstructable.
    nChains*: int64
    sumDP*: float
    sumDP2*: float
    sumMass*: float
    sumMass2*: float
    sumMass3*: float
    dpn*: float
    dpw*: float
    mn*: float
    mw*: float
    mz*: float
    pdi*: float

  ActionResult* = object
    target*: string
    requested*: string
    before*: float
    after*: float
    message*: string
    outputWritten*: string
    hasNumeric*: bool
    stateChanged*: bool


proc sanitizeColumnToken*(raw: string): string =
  ## Stable, header-safe token for dynamically generated text table column names.
  ## It is deliberately derived from the declared reaction text.
  for c in raw:
    if (c >= 'A' and c <= 'Z') or (c >= 'a' and c <= 'z') or
       (c >= '0' and c <= '9') or c == '_':
      result.add c
    else:
      result.add '_'
  while result.contains("__"):
    result = result.replace("__", "_")
  if result.len == 0:
    result = "channel"

proc channelProgressIds*(m: Model): seq[string] =
  ## Deterministic identifiers in model/declaration order.
  ## The model file remains the authoritative reaction mapping.
  var used = initTable[string, int]()
  for ch in m.channels:
    let base = sanitizeColumnToken(ch.expr)
    let occurrence = used.getOrDefault(base, 0) + 1
    used[base] = occurrence
    if occurrence == 1:
      result.add base
    else:
      result.add base & "_" & $occurrence

proc deadPoolIds*(m: Model): seq[int] =
  ## Dead pools in declaration order for history-moment reporting.
  for pid, p in m.pools:
    if p.kind == pkDead:
      result.add pid

proc hash*(x: StructKey): Hash =
  result = x.eg1.hash !& x.dp.hash !& x.eg2.hash !& x.formedBy.hash
  result = !$result

proc formationKindName*(kind: FormationKind): string =
  ## Stable text written to the formed_by column in chains.csv.
  case kind
  of fbUnknown: "unknown"
  of fbInit: "init"
  of fbTransferM: "transfer_m"
  of fbTermC: "term_c"
  of fbTermD_H: "term_d_H"
  of fbTermD_U: "term_d_U"
  of fbTermX: "term_x"
  of fbTransferH: "transfer_h"

proc fail*(lineNo: int; msg: string) =
  ## Internal parser/validation failure. Do not terminate here: raise a
  ## catchable exception so parser checks are unit-testable and the CLI alone
  ## decides the process exit code. The "ERROR" prefix is added centrally by
  ## the CLI entry point (slimmc.nim), not here, so that runtime exceptions
  ## raised outside of fail() (e.g. in slimmc_kmc.nim) get the same prefix.
  if lineNo > 0:
    raise newException(ValueError, "line " & $lineNo & ": " & msg)
  else:
    raise newException(ValueError, msg)

proc warn*(lineNo: int; msg: string) =
  if lineNo > 0:
    stderr.writeLine("WARNING line " & $lineNo & ": " & msg)
  else:
    stderr.writeLine("WARNING: " & msg)

proc num*(x: float): string =
  ## Locale-independent float64 formatting with enough significant digits
  ## for safe round-trip through text table and other machine-readable outputs.
  formatFloat(x, ffDefault, 17)

proc numDec*(x: float; p: int = 4): string =
  formatFloat(x, ffDecimal, p)

proc conc*(n: int64; V: float): float =
  float(n) / (NA * V)

proc countFromConc*(c: float; V: float): int64 =
  checkedCountFromConc(c, NA, V, "concentration")

proc parseF*(s0: string; lineNo: int = 0; what: string = "number"): float =
  ## Parse a locale-independent floating-point number.
  let s = s0.strip()
  if s.len == 0:
    fail(lineNo, "empty " & what)

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
    fail(lineNo, "invalid " & what & ": " & s0)

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
      fail(lineNo, "invalid exponent in " & what & ": " & s0)

    expVal *= expSign

  if i != s.len:
    fail(lineNo, "invalid trailing characters in " & what & ": " & s0)

  result = sign * (intPart + fracPart / fracScale) * pow(10.0, float(expVal))

proc parseI64*(s0: string; lineNo: int = 0; what: string = "integer"): int64 =
  try:
    result = checkedParseInt64(s0, what)
  except ValueError as exc:
    fail(lineNo, exc.msg)

proc endGroupName*(speciesName: string): string =
  ## A terminal underscore marks a small radical in model syntax.
  ## The corresponding polymer end group omits that underscore.
  if speciesName.len > 1 and speciesName[^1] == '_':
    result = speciesName[0 ..< speciesName.len - 1]
  else:
    result = speciesName

proc transferMonomerEndGroupName*(monomerName: string): string =
  ## Provenance label for a chain born by transfer to the declared monomer.
  ## Example: M -> M_tr. This is not an additional repeat unit and is not a
  ## literal molecular formula for the chain end; the monomer-derived unit is
  ## already included in DP = 1 when the chain is created.
  result = endGroupName(monomerName) & "_tr"

proc ensureEg*(m: var Model; name: string): int =
  if m.egByName.hasKey(name):
    return m.egByName[name]

  result = m.egNames.len
  m.egNames.add name
  m.egByName[name] = result
  m.egMass.add 0.0
  m.egMassKnown.add false
  m.egMassSource.add "unknown"

proc setEgMass*(m: var Model; name: string; mass: float; source: string) =
  let id = m.ensureEg(name)
  m.egMass[id] = mass
  m.egMassKnown[id] = true
  m.egMassSource[id] = source

proc ensureBuiltinEg*(m: var Model; name: string; mass: float; source: string): int =
  result = m.ensureEg(name)
  if not m.egMassKnown[result]:
    m.egMass[result] = mass
    m.egMassKnown[result] = true
    m.egMassSource[result] = source

proc initModel*(): Model =
  result.V = 0.0
  result.initVolumeMl = 0.0
  result.currentVolumeMl = 0.0
  result.hasInitVolume = false
  result.T = DefaultTemperatureK
  result.tEnd = -1.0
  result.maxEvents = DefaultMaxSteps
  result.whenCheckEvents = DefaultWhenCheckEvents
  result.seed = DefaultSeed

  result.modelSourceFile = ""
  result.runId = ""
  result.outputDir = ""

  result.description = ""
  result.hasDescription = false
  result.variables = @[]
  result.massModel = DefaultMassModel
  result.dpMax = DefaultDpMax
  result.sequenceMode = DefaultSequenceMode

  result.species = @[]
  result.speciesByName = initTable[string, int]()
  result.monomerId = -1
  result.feeds = @[]
  result.feedByName = initTable[string, int]()

  result.pools = @[]
  result.poolByName = initTable[string, int]()

  result.rates = @[]
  result.rateByName = initTable[string, int]()

  result.egNames = @[]
  result.egByName = initTable[string, int]()
  result.egMass = @[]
  result.egMassKnown = @[]
  result.egMassSource = @[]

  result.channels = @[]
  result.scheduledActions = @[]
  result.conditionalActions = @[]
  result.rawReactions = @[]
  result.rawConditionalActions = @[]

  result.egH = result.ensureBuiltinEg("H", 1.008, "builtin")
  result.egU = result.ensureBuiltinEg("U", -1.008, "builtin")
  result.egActive = result.ensureBuiltinEg("ACTIVE", 0.0, "builtin")

proc rateValue*(m: Model; kId: int): float =
  let r = m.rates[kId]

  case r.kind
  of rkFixed:
    result = r.kConst
  of rkArr:
    result = r.Apre * exp(-r.Ea / (Rgas * m.T))

proc channelKindName*(k: ChannelKind): string =
  case k
  of chRxnUni: result = "rxn_uni"
  of chRxnBiDiff: result = "rxn_bi_diff"
  of chRxnBiSame: result = "rxn_bi_same"
  of chMacroInit: result = "macro_init"
  of chMacroProp: result = "macro_prop"
  of chMacroDeprop: result = "macro_deprop"
  of chMacroTermC: result = "macro_term_c"
  of chMacroTermD: result = "macro_term_d"
  of chMacroTermX: result = "macro_term_x"
  of chMacroTransferH: result = "macro_transfer_h"
  of chMacroTransferM: result = "macro_transfer_m"

proc actionKindName*(a: ActionKind): string =
  case a
  of eaPrint: result = "print"
  of eaPrintInfo: result = "print_info"
  of eaSave: result = "save"
  of eaSaveChains: result = "save_chains"
  of eaStop: result = "stop"
  of eaSetK: result = "set_k"
  of eaAddK: result = "add_k"
  of eaSetTemp: result = "set_temp"
  of eaAddTemp: result = "add_temp"
  of eaSetC: result = "set_c"
  of eaAddC: result = "add_c"
  of eaFeed: result = "feed"
  of eaPrintMemory: result = "print_memory"

proc comparisonName*(c: ComparisonKind): string =
  case c
  of coGreater: result = ">"
  of coLess: result = "<"

proc conditionalObservableName*(m: Model; a: AtomicCondition): string =
  case a.observable
  of woConversion: result = "X"
  of woSpeciesConc: result = "c:" & m.species[a.speciesId].name

proc conditionalObservableName*(m: Model; e: ConditionalAction): string =
  conditionalObservableName(m, AtomicCondition(observable: e.observable, speciesId: e.speciesId,
    comparison: e.comparison, threshold: e.threshold))

proc poolStats*(pool: seq[Chain]): tuple[n: int64, dpSum: int64, dpn: float, dpw: float] =
  let nChains = int64(pool.len)
  var sumDP: int64 = 0
  var sumDP2 = 0.0

  for ch in pool:
    sumDP += int64(ch.dp)
    sumDP2 += checkedSquareAsFloat(int64(ch.dp), "chain DP")

  var dpn = 0.0
  var dpw = 0.0

  if nChains > 0:
    dpn = float(sumDP) / float(nChains)

  if sumDP > 0:
    dpw = sumDP2 / float(sumDP)

  result = (n: nChains, dpSum: sumDP, dpn: dpn, dpw: dpw)

proc calcMTotal*(m: Model; s: State): int64 =
  if m.monomerId >= 0:
    result = s.n[m.monomerId]
  else:
    result = 0

  for pool in s.pools:
    for ch in pool:
      result += int64(ch.dp)
