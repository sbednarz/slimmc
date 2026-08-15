# copo_types.nim
# Shared model, channel and simulation-state types for slimmc-copo.

import math, tables

const
  AppName* = "slimmc-copo"
  SlimmcVersion* {.strdefine.} = "5.0.3"
  AppVersion* = SlimmcVersion
  NA* = 6.02214076e23
  Rgas* = 8.31446261815324
  InfTime* = 1.0e300
  Eps* = 1.0e-12
  TimeRelEps* = 8.0e-15

# v0.2 keeps the v0.1 terminal core and adds:
# - transfer to CTA/solvent-like species with optional reinitiation through normal init
# - terminal depropagation with sequence/pop bookkeeping
# - terminal/penultimate bookkeeping for active pools
# - live chains store full linear sequence
# - dead chains are compact summaries

proc timeTolerance*(a, b: float): float =
  ## Mixed absolute/relative tolerance for simulation-time comparisons.
  ## Eps preserves sub-second behavior; the relative term stays several
  ## floating-point ulps wide when t becomes very large.
  max(Eps, TimeRelEps * max(abs(a), abs(b)))

proc timeClose*(a, b: float): bool =
  abs(a - b) <= timeTolerance(a, b)

type
  PoolKind* = enum
    pkActive, pkDead

  FormationKind* = enum
    fbUnknown,
    fbInit,
    fbTermC,
    fbTermD_H,
    fbTermD_U,
    fbTermX,
    fbTransfer,
    fbTransferM

  SnapshotRecord* = object
    snapshotId*: int64
    kmcEvent*: int64
    t*: float
    level*: string
    isFinal*: bool
    reason*: string
    parameterStateId*: int64
    stateRevision*: int64

  RateKind* = enum
    rkFixed, rkArr

  SmallKind* = enum
    skSpecies,
    skMonomer

  SmallRef* = object
    kind*: SmallKind
    id*: int
    stoich*: int

  SmallDelta* = object
    kind*: SmallKind
    id*: int
    delta*: int

  ChannelKind* = enum
    chRxnUni,
    chRxnBiDiff,
    chRxnBiSame,
    chMacroInit,
    chMacroProp,
    chMacroTermC,
    chMacroTermD,
    chMacroTermX,
    chMacroTransfer,
    chMacroTransferM,
    chMacroDeprop

  MassModel* = enum
    mmRepeatUnits,
    mmWithEndgroups

  ActionKind* = enum
    eaPrint,
    eaPrintInfo,
    eaSave,
    eaSaveChains,
    eaStop,
    eaPrintMemory,
    eaSetK,
    eaAddK,
    eaSetTemp,
    eaAddTemp,
    eaSetC,
    eaAddC,
    eaFeed

  ConditionalObservableKind* = enum
    woTotalConversion,
    woMonomerConversion,
    woSpeciesConc,
    woMonomerConc

  ComparisonKind* = enum
    coGreater,
    coLess

  MonomerDef* = object
    name*: string
    c0*: float
    mw*: float

  SpeciesDef* = object
    name*: string
    c0*: float

  EndGroupDef* = object
    name*: string
    mw*: float

  PoolDef* = object
    name*: string
    kind*: PoolKind

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

  VarDef* = object
    kind*: string
    name*: string
    value*: float
    unit*: string

  LinearSequence* = object
    ## v0.1 uses one byte per mer for clarity. The type boundary is deliberate:
    ## later versions can replace this with packed 1/2-bit storage.
    data*: seq[uint8]

  BlockCount* = object
    ## Exact run/block-length counter for one monomer identity.
    ## Stored on dead summaries so block statistics survive sequenceText cleanup.
    monomerId*: int
    length*: int32
    count*: int64

  LiveChain* = object
    id*: int64
    left_end*: string
    mers*: LinearSequence
    right_end*: string
    nMer*: seq[int32]
    dp*: int32
    mass*: float
    last*: int
    prev*: int
    formedBy*: FormationKind

  DeadSummary* = object
    left_end*: string
    right_end*: string
    dp*: int32
    nMer*: seq[int32]
    mass*: float
    formedBy*: FormationKind
    firstMer*: int
    penultimateMer*: int
    lastMer*: int
    dyads*: seq[int64]
    triads*: seq[int64]
    blockCounts*: seq[BlockCount]
    sequenceStored*: bool
    sequenceText*: string
    count*: int64

  KmcChannel* = object
    name*: string
    lineNo*: int
    kind*: ChannelKind
    kId*: int

    # elementary rxn: small reactants (lhs, always consumed) and
    # products (rhs, formed only with probability `efficiency` -- items
    # 79/80/81). smallDeltas (net combined lhs+rhs) is kept only for
    # any remaining internal callers; the runtime uses smallReactants/
    # smallProducts directly so consumption and production can be
    # applied independently.
    smallReactants*: seq[SmallRef]
    smallProducts*: seq[SmallRef]
    smallDeltas*: seq[SmallDelta]
    efficiency*: float

    # init: speciesId + monomerId -> poolOut
    speciesId*: int
    monomerId*: int

    # prop/term/deprop/transfer: pool1/pool2 -> poolOut
    pool1*: int
    pool2*: int
    poolOut*: int

    # transfer: active pool + speciesId -> dead pool + speciesOutId
    speciesOutId*: int

  FeedDef* = object
    name*: string
    monomerConcentrations*: seq[float]
    speciesConcentrations*: seq[float]

  ScheduledAction* = object
    startTime*: float
    nextTime*: float
    period*: float
    repeat*: bool
    remaining*: int64
    active*: bool
    action*: ActionKind
    args*: seq[string]
    lineNo*: int

  AtomicCondition* = object
    observable*: ConditionalObservableKind
    targetId*: int
    comparison*: ComparisonKind
    threshold*: float

  ConditionalAction* = object
    active*: bool
    conditions*: seq[AtomicCondition]
    # Convenience fields mirror the first condition for direct access.
    observable*: ConditionalObservableKind
    targetId*: int
    comparison*: ComparisonKind
    threshold*: float
    action*: ActionKind
    args*: seq[string]
    lineNo*: int

  MemoryPolicy* = object
    limitBytes*: int64
    hasLimit*: bool
    snapshotOnLimit*: bool
    stopOnLimit*: bool

  Model* = object
    V*: float
    initVolumeMl*: float
    currentVolumeMl*: float
    hasInitVolume*: bool
    T*: float
    t_end*: float
    max_steps*: int64
    whenCheckEvents*: int64
    seed*: int64
    output_dir*: string
    output_dirWasSet*: bool
    modelFile*: string
    modelStem*: string
    startedAt*: string  # items 10/11/15: set once in runSimulation, read by every text table writer's metadata preamble
    description*: string
    hasDescription*: bool
    variables*: seq[VarDef]
    dp_max*: int64
    sequence_mode*: string
    mass_model*: MassModel

    monomers*: seq[MonomerDef]
    monomerByName*: Table[string, int]

    species*: seq[SpeciesDef]
    speciesByName*: Table[string, int]

    feeds*: seq[FeedDef]
    feedByName*: Table[string, int]

    endgroups*: seq[EndGroupDef]
    endgroupByName*: Table[string, int]

    pools*: seq[PoolDef]
    poolByName*: Table[string, int]
    deadPoolId*: int
    poolTerminalMer*: seq[int]
    poolPenultimateMer*: seq[int]

    rates*: seq[RateDef]
    rateByName*: Table[string, int]

    channels*: seq[KmcChannel]
    scheduledActions*: seq[ScheduledAction]
    conditionalActions*: seq[ConditionalAction]
    rawConditionalActions*: seq[RawLine]
    rawLines*: seq[RawLine]
    memoryPolicy*: MemoryPolicy

  State* = object
    t*: float
    kmcEvent*: int64
    actionNo*: int64
    stateRevision*: int64
    parameterStateId*: int64
    snapshotId*: int64
    # Monotonic identifier of a compressed Storage chain record.
    chainRecordId*: int64
    snapshots*: seq[SnapshotRecord]
    nextChainId*: int64
    # Set once any snapshot (save or save_chains) has actually been written.
    # `chainsWritten` additionally tracks whether save_chains specifically has
    # fired. Runtime bookkeeping uses both flags.
    # advertises a file the model never asked to be written -- see the parity
    # fix in runSimulation that removed the old unconditional start/end
    # snapshot (which used to guarantee every output file always existed
    # regardless of what the model requested).
    anySnapshotWritten*: bool
    chainsWritten*: bool
    channelTraceRowsWritten*: int64
    channelTraceTruncated*: bool
    feedMonomerRemainders*: seq[seq[float]]
    feedSpeciesRemainders*: seq[seq[float]]

    # Incremental live-chain eligibility accounting for the SSA hot path.
    # poolEligibilityTrackedLen is also the validity sentinel: states built or
    # mutated manually by tests/users fall back to exact scans until the engine
    # itself starts tracking that pool.
    poolEligibleCounts*: seq[int64]
    poolPropagatableCounts*: seq[int64]
    poolEligibilityTrackedLen*: seq[int64]
    channelDepropEligibleCounts*: seq[int64]

    # Slimmc Storage v1 core backend (stage K1).
    storageV1RunDir*: string
    storageV1StartedAt*: string
    storageV1SnapshotIds*: seq[uint64]
    storageV1StateRevisions*: seq[uint64]  # internal deduplication key; not serialized
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
    storageV1ChannelSnapshotIds*: seq[uint64]
    storageV1ChannelIds*: seq[uint32]
    storageV1ChannelEventCounts*: seq[uint64]
    storageV1ChannelProductiveCounts*: seq[uint64]
    storageV1ChannelNonproductiveCounts*: seq[uint64]
    storageV1PropensitySnapshotIds*: seq[uint64]
    storageV1PropensityChannelIds*: seq[uint32]
    storageV1PropensityValues*: seq[float64]
    storageV1TotalPropensities*: seq[float64]

    # Slimmc Storage v1 stage K2: complete kinetic sets and executed actions.
    storageV1KineticSetIds*: seq[uint64]
    storageV1KineticStartEvents*: seq[uint64]
    storageV1KineticStartTimes*: seq[float64]
    storageV1KineticHasSourceAction*: seq[bool]
    storageV1KineticSourceActionIds*: seq[uint64]
    storageV1KineticValueSetIds*: seq[uint64]
    storageV1KineticParameterIds*: seq[uint32]
    storageV1KineticValues*: seq[float64]
    storageV1ActionIds*: seq[uint64]
    storageV1ActionKmcEvents*: seq[uint64]
    storageV1ActionTimes*: seq[float64]
    storageV1ActionSourceLines*: seq[uint32]
    storageV1ActionTypeIds*: seq[uint32]
    storageV1ActionTriggerTypeIds*: seq[uint32]
    storageV1ActionScheduledTimes*: seq[float64]
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
    storageV1SpeciesBalanceSnapshotIds*: seq[uint64]
    storageV1SpeciesBalanceEntityIds*: seq[uint32]
    storageV1SpeciesInitialMoles*: seq[float64]
    storageV1SpeciesDosedMoles*: seq[float64]
    storageV1SpeciesTotalMoles*: seq[float64]
    storageV1SpeciesFreeMoles*: seq[float64]
    storageV1SpeciesConsumedMoles*: seq[float64]
    storageV1ConditionRecordIds*: seq[uint64]
    storageV1ConditionActionIds*: seq[uint64]
    storageV1ConditionIndexes*: seq[uint32]
    storageV1ConditionObservableIds*: seq[uint32]
    storageV1ConditionOperatorIds*: seq[uint32]
    storageV1ConditionThresholds*: seq[float64]
    storageV1ConditionObservedValues*: seq[float64]
    storageV1ConditionMet*: seq[bool]

    # Slimmc Storage v1 stage K3: chain-resolved copolymer data.
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
    storageV1ObservedEndGroups*: seq[string]
    storageV1ChainHasFirst*: seq[bool]
    storageV1ChainFirstIds*: seq[uint32]
    storageV1ChainHasPenultimate*: seq[bool]
    storageV1ChainPenultimateIds*: seq[uint32]
    storageV1ChainHasLast*: seq[bool]
    storageV1ChainLastIds*: seq[uint32]
    storageV1ChainHasSequence*: seq[bool]
    storageV1ChainSequenceOffsets*: seq[uint64]
    storageV1ChainSequenceLengths*: seq[uint64]
    storageV1CompositionRecordIds*: seq[uint64]
    storageV1CompositionMonomerIds*: seq[uint32]
    storageV1CompositionUnitCounts*: seq[uint64]
    storageV1SequenceSymbols*: seq[uint32]
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
    storageV1MemorySnapshotIds*: seq[uint64]
    storageV1MemoryLiveChains*: seq[uint64]
    storageV1MemoryDeadRecords*: seq[uint64]
    storageV1MemoryDeadChains*: seq[uint64]
    storageV1MemoryStoredLiveMers*: seq[uint64]
    storageV1MemoryStoredDeadMers*: seq[uint64]
    storageV1MemoryLiveSeqBytes*: seq[uint64]
    storageV1MemoryLiveObjectBytes*: seq[uint64]
    storageV1MemoryDeadRecordBytes*: seq[uint64]
    storageV1MemoryTotalEstBytes*: seq[uint64]

    # Storage microstructure parity: aggregate motifs and block histograms.
    storageV1MotifSnapshotIds*: seq[uint64]
    storageV1MotifOrders*: seq[uint32]
    storageV1MotifIds*: seq[uint32]
    storageV1MotifCounts*: seq[uint64]
    storageV1BlockSnapshotIds*: seq[uint64]
    storageV1BlockMonomerIds*: seq[uint32]
    storageV1BlockLengths*: seq[uint64]
    storageV1BlockCounts*: seq[uint64]

    speciesN*: seq[int64]
    monomerN*: seq[int64]
    monomerN0*: seq[int64]
    monomerBalance*: seq[int64]  # item 24: independent ledger, updated only by elementary rxns and set_c/add_c -- NOT by init/prop/deprop/transfer_m, which only move mass between free and incorporated pools (mass-conserving by construction). Should equal NAME_count + poly_NAME_count if bookkeeping is consistent.
    monomerExternalN*: seq[int64]
    monomerDosedN*: seq[int64]
    speciesExternalN*: seq[int64]
    speciesDosedN*: seq[int64]

    livePools*: seq[seq[LiveChain]]
    deadChains*: seq[DeadSummary]
    channelFires*: seq[int64]
    channelSuccesses*: seq[int64]  # item 82: for f<1 rxn channels, fires = successes + failures (exact invariant)
    channelFailures*: seq[int64]
    sumA0Tau*: float
    sumA0TauSq*: float
    countA0Tau*: int64
    stopRequested*: bool
    stopLineNo*: int
    stopCheckSource*: string
    stopConditions*: seq[AtomicCondition]
    stopActualValues*: seq[float]

    storageTraceKmcEvents*: seq[uint64]
    storageTraceTimes*: seq[float64]
    storageTraceDt*: seq[float64]
    storageTraceChannelIds*: seq[uint32]
    storageTraceRates*: seq[float64]
    storageTracePropensities*: seq[float64]
    storageTraceTotalPropensities*: seq[float64]

  RunOptions* = object
    checkOnly*: bool
    debug*: bool
    traceChannelsLimit*: int64

  MemoryEstimate* = object
    liveChains*: int64
    deadSummaries*: int64
    storedLiveMers*: int64
    storedDeadMers*: int64
    liveSeqBytes*: int64
    liveObjectBytes*: int64
    deadSummaryBytes*: int64
    totalBytes*: int64

  MomentStats* = object
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
    pdi*: float
    mz*: float

proc rateValue*(m: Model; kId: int): float =
  let r = m.rates[kId]
  case r.kind
  of rkFixed:
    result = r.kConst
  of rkArr:
    result = r.Apre * exp(-r.Ea / (Rgas * m.T))

proc formationName*(f: FormationKind): string =
  case f
  of fbUnknown: "unknown"
  of fbInit: "init"
  of fbTermC: "term_c"
  of fbTermD_H: "term_d_H"
  of fbTermD_U: "term_d_U"
  of fbTermX: "term_x"
  of fbTransfer: "transfer"
  of fbTransferM: "transfer_m"

proc poolKindName*(k: PoolKind): string =
  case k
  of pkActive: "active"
  of pkDead: "dead"

proc channelKindName*(k: ChannelKind): string =
  case k
  of chRxnUni: "rxn_uni"
  of chRxnBiDiff: "rxn_bi_diff"
  of chRxnBiSame: "rxn_bi_same"
  of chMacroInit: "macro_init"
  of chMacroProp: "macro_prop"
  of chMacroTermC: "macro_term_c"
  of chMacroTermD: "macro_term_d"
  of chMacroTermX: "macro_term_x"
  of chMacroTransfer: "macro_transfer"
  of chMacroTransferM: "macro_transfer_m"
  of chMacroDeprop: "macro_deprop"

proc actionKindName*(a: ActionKind): string =
  case a
  of eaPrint: result = "print"
  of eaPrintInfo: result = "print_info"
  of eaSave: result = "save"
  of eaSaveChains: result = "save_chains"
  of eaStop: result = "stop"
  of eaPrintMemory: result = "print_memory"
  of eaSetK: result = "set_k"
  of eaAddK: result = "add_k"
  of eaSetTemp: result = "set_temp"
  of eaAddTemp: result = "add_temp"
  of eaSetC: result = "set_c"
  of eaAddC: result = "add_c"
  of eaFeed: result = "feed"

proc comparisonName*(c: ComparisonKind): string =
  case c
  of coGreater: result = ">"
  of coLess: result = "<"

proc conditionalObservableName*(m: Model; a: AtomicCondition): string =
  case a.observable
  of woTotalConversion:
    result = "X"
  of woMonomerConversion:
    result = "X:" & m.monomers[a.targetId].name
  of woSpeciesConc:
    result = "c:" & m.species[a.targetId].name
  of woMonomerConc:
    result = "c:" & m.monomers[a.targetId].name

proc conditionalObservableName*(m: Model; e: ConditionalAction): string =
  conditionalObservableName(m, AtomicCondition(observable: e.observable, targetId: e.targetId,
    comparison: e.comparison, threshold: e.threshold))

proc mass_modelName*(mm: MassModel): string =
  case mm
  of mmRepeatUnits: "repeat_units"
  of mmWithEndgroups: "with_end_groups"
