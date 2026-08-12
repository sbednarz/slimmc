## Slimmc Storage schema builder.
## Schema construction is isolated and compiled without optimization because
## it is non-hot code and large JSON literals otherwise dominate build time.

{.push optimization: none.}

import std/[json, tables]
import ../copo_types
import ../../../common/[results_types, storage_schema]

proc schemaRecords*(m: Model; observedEndGroups: seq[string] = @[]): seq[JsonNode] {.codegenDecl: "__attribute__((optimize(\"O0\"))) $# $#$#".} =
  result.add storageRecord([("record_type", %("schema_header")), ("schema_name", %(StorageName)), ("schema_version", %(StorageFormatVersion)), ("byte_order", %("little")), ("column_storage", %("npy")), ("npy_format_version", %("1.0"))])
  for t in ["snapshots","state","channel_events","channel_propensities","actions","action_conditions","feed_events","monomer_balance","species_balance","kinetic_parameters/sets","kinetic_parameters/values","chains","chain_composition","sequences","moments","memory","microstructure_motifs","block_statistics","diagnostics/channel_trace"]:
    result.add storageRecord([("record_type", %("table")), ("name", %(t)), ("required", %(true))])
  template col(t,n,d,u: string) = result.add storageRecord([("record_type", %("column")), ("table", %(t)), ("name", %(n)), ("file", %(n & ".npy")), ("dtype", %(d)), ("unit", %(u)), ("required", %(true))])
  col("snapshots","snapshot_id","uint64","")
  col("snapshots","time","float64","s")
  col("snapshots","kmc_event","uint64","")
  col("snapshots","snapshot_reason_id","uint32","")
  col("snapshots","is_final","bool","")
  col("snapshots","has_chains","bool","")
  col("snapshots","has_sequences","bool","")
  col("snapshots","kinetic_parameter_set_id","uint64","")
  col("snapshots","volume_mL","float64","mL")
  col("snapshots","kmc_volume_L","float64","L")
  col("snapshots","chain_count_live","uint64","")
  col("snapshots","chain_count_dead","uint64","")
  col("snapshots","chain_count_total","uint64","")
  col("state","snapshot_id","uint64","")
  col("state","entity_id","uint32","")
  col("state","count","uint64","")
  col("state","moles","float64","mol")
  col("state","concentration","float64","mol/L")
  col("channel_events","snapshot_id","uint64","")
  col("channel_events","channel_id","uint32","")
  col("channel_events","event_count","uint64","")
  col("channel_events","productive_event_count","uint64","")
  col("channel_events","nonproductive_event_count","uint64","")
  col("diagnostics/channel_trace","kmc_event","uint64","")
  col("diagnostics/channel_trace","time","float64","s")
  col("diagnostics/channel_trace","dt","float64","s")
  col("diagnostics/channel_trace","channel_id","uint32","")
  col("diagnostics/channel_trace","rate","float64","engine_native")
  col("diagnostics/channel_trace","propensity","float64","s^-1")
  col("diagnostics/channel_trace","total_propensity","float64","s^-1")
  col("channel_propensities","snapshot_id","uint64","")
  col("channel_propensities","channel_id","uint32","")
  col("channel_propensities","propensity","float64","s^-1")
  col("channel_propensities","total_propensity","float64","s^-1")
  for (n,d) in [("action_id","uint64"),("kmc_event","uint64"),("time","float64"),("source_line","uint32"),("action_type_id","uint32"),("trigger_type_id","uint32"),("scheduled_time","float64"),("target_id","uint32"),("requested_value","float64"),("before_value","float64"),("after_value","float64"),("state_changed","bool"),("output_written","bool"),("has_snapshot","bool"),("snapshot_id","uint64"),("has_kinetic_parameter_set","bool"),("kinetic_parameter_set_id","uint64")]:
    col("actions",n,d,"")
  for (n,d) in [("condition_record_id","uint64"),("action_id","uint64"),("condition_index","uint32"),("observable_id","uint32"),("operator_id","uint32"),("threshold","float64"),("observed_value","float64"),("condition_met","bool")]:
    col("action_conditions",n,d,"")
  for (n,d) in [("kinetic_parameter_set_id","uint64"),("start_kmc_event","uint64"),("start_time","float64"),("has_source_action","bool"),("source_action_id","uint64")]:
    col("kinetic_parameters/sets",n,d,"")
  for (n,d) in [("kinetic_parameter_set_id","uint64"),("kinetic_parameter_id","uint32"),("value","float64")]:
    col("kinetic_parameters/values",n,d,"")
  for (n,d) in [("chain_record_id","uint64"),("snapshot_id","uint64"),("population_id","uint32"),("pool_id","uint32"),("origin_id","uint32"),("dp","uint64"),("molar_mass","float64"),("count","uint64"),("moles","float64"),("concentration","float64"),("left_end_id","uint32"),("right_end_id","uint32"),("has_first_monomer","bool"),("first_monomer_id","uint32"),("has_penultimate_monomer","bool"),("penultimate_monomer_id","uint32"),("has_last_monomer","bool"),("last_monomer_id","uint32"),("has_sequence","bool"),("sequence_offset","uint64"),("sequence_length","uint64")]:
    col("chains",n,d,"")
  for (n,d) in [("chain_record_id","uint64"),("monomer_id","uint32"),("unit_count","uint64")]: col("chain_composition",n,d,"")
  col("sequences","symbols","uint32","")
  for (n,d) in [("snapshot_id","uint64"),("population_scope_id","uint32"),("mass_basis_id","uint32"),("chain_count","uint64"),("sum_dp","float64"),("sum_dp2","float64"),("dp_n","float64"),("dp_w","float64"),("sum_molar_mass","float64"),("sum_molar_mass2","float64"),("sum_molar_mass3","float64"),("mn","float64"),("mw","float64"),("mz","float64"),("dispersity","float64")]: col("moments",n,d,"")
  for (n,d) in [("snapshot_id","uint64"),("live_chains","uint64"),("dead_records","uint64"),("dead_chains","uint64"),("stored_live_mers","uint64"),("stored_dead_mers","uint64"),("live_seq_B","uint64"),("live_object_B","uint64"),("dead_record_B","uint64"),("total_est_B","uint64")]: col("memory",n,d,"")
  for (n,d) in [("snapshot_id","uint64"),("motif_order","uint32"),("motif_id","uint32"),("count","uint64")]: col("microstructure_motifs",n,d,"")
  for (n,d) in [("snapshot_id","uint64"),("monomer_id","uint32"),("block_length","uint64"),("block_count","uint64")]: col("block_statistics",n,d,"")
  for id,name in ["initial","scheduled","action","manual","final"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("snapshot_reasons")), ("id", %(id)), ("name", %(name))])
  var eid=0
  for monomerId, mo in m.monomers:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %("monomer_" & mo.name)), ("kind", %("monomer"))]); inc eid
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("monomers")), ("id", %(monomerId)), ("name", %(mo.name)), ("molar_mass_increment", %(mo.mw))])
  for sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %("species_" & sp.name)), ("kind", %("species"))]); inc eid
  for name in ["live_chains","dead_chains"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %(name)), ("kind", %("aggregate"))]); inc eid
  for id,ch in m.channels:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("channels")), ("id", %(id)), ("name", %(ch.name)), ("kind", %(channelKindName(ch.kind)))])
  for id,name in ["unknown","at","every","when"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_triggers")), ("id", %(id)), ("name", %(name))])
  for a in ActionKind:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_types")), ("id", %(ord(a))), ("name", %(actionKindName(a)))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("condition_operators")), ("id", %(1)), ("name", %(">"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("condition_operators")), ("id", %(2)), ("name", %("<"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(0)), ("name", %("temperature_K")), ("kind", %("temperature")), ("unit", %("K"))])
  var kpid=1
  for r in m.rates:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kpid)), ("name", %(r.name)), ("kind", %("rate_constant")), ("unit", %("model_defined"))]); inc kpid
    if r.declaredArrhenius:
      result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kpid)), ("name", %(r.name & "_A")), ("kind", %("arrhenius_A")), ("unit", %("model_defined"))]); inc kpid
      result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kpid)), ("name", %(r.name & "_Ea")), ("kind", %("arrhenius_Ea")), ("unit", %("J/mol"))]); inc kpid
  for id,name in ["live","dead"]: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_populations")), ("id", %(id)), ("name", %(name))])
  for id,p in m.pools: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_pools")), ("id", %(id)), ("name", %(p.name))])
  for f in FormationKind: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_origins")), ("id", %(ord(f))), ("name", %(formationName(f)))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_end_types")), ("id", %(0)), ("name", %("not_applicable"))])
  for id,name in observedEndGroups:
    var rec = storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_end_types")), ("id", %(id+1)), ("name", %(name))])
    if m.endgroupByName.hasKey(name):
      rec["molar_mass_contribution"] = %m.endgroups[m.endgroupByName[name]].mw
    result.add rec
  for id,mo in m.monomers: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("sequence_symbols")), ("id", %(id)), ("name", %(mo.name)), ("kind", %("monomer")), ("monomer_id", %(id))])
  col("feed_events","action_id","uint64","")
  col("feed_events","feed_id","uint32","")
  col("feed_events","dose_mL","float64","mL")
  col("feed_events","volume_before_mL","float64","mL")
  col("feed_events","volume_after_mL","float64","mL")
  col("feed_events","kmc_volume_before_L","float64","L")
  col("feed_events","kmc_volume_after_L","float64","L")
  col("monomer_balance","snapshot_id","uint64","")
  col("monomer_balance","monomer_id","uint32","")
  col("monomer_balance","initial_moles","float64","mol")
  col("monomer_balance","introduced_moles","float64","mol")
  col("monomer_balance","free_moles","float64","mol")
  col("monomer_balance","incorporated_moles","float64","mol")
  col("monomer_balance","conversion","float64","")
  col("species_balance","snapshot_id","uint64","")
  col("species_balance","entity_id","uint32","")
  col("species_balance","initial_moles","float64","mol")
  col("species_balance","dosed_moles","float64","mol")
  col("species_balance","total_moles","float64","mol")
  col("species_balance","free_moles","float64","mol")
  col("species_balance","consumed_moles","float64","mol")
  var balanceId = 0
  for mo in m.monomers:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("balance_entities")), ("id", %(balanceId)), ("name", %(mo.name)), ("kind", %("monomer"))]); inc balanceId
  for sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("balance_entities")), ("id", %(balanceId)), ("name", %(sp.name)), ("kind", %("species"))]); inc balanceId
  for fid, feed in m.feeds:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("feeds")), ("id", %(fid)), ("name", %(feed.name))])
  for id,name in ["all","live","dead"]: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("population_scope")), ("id", %(id)), ("name", %(name))])
  for id,name in ["repeat_units","with_end_groups"]: result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("mass_bases")), ("id", %(id)), ("name", %(name))])
  var motifId = 0
  for a in 0 ..< m.monomers.len:
    for b in 0 ..< m.monomers.len:
      result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("microstructure_dyads")), ("id", %(motifId)), ("name", %(m.monomers[a].name & "|" & m.monomers[b].name)), ("order", %(2))]); inc motifId
  motifId = 0
  for a in 0 ..< m.monomers.len:
    for b in 0 ..< m.monomers.len:
      for c in 0 ..< m.monomers.len:
        result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("microstructure_triads")), ("id", %(motifId)), ("name", %(m.monomers[a].name & "|" & m.monomers[b].name & "|" & m.monomers[c].name)), ("order", %(3))]); inc motifId
  result.add storageRecord([("record_type", %("rule")), ("name", %("state_dense")), ("description", %("Each snapshot contains every state entity in ascending entity_id order."))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("channel_event_sum")), ("description", %("For each snapshot, kmc_event equals the sum of cumulative event_count over channels."))])


{.pop.}
