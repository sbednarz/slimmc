## Slimmc Storage schema builder.
## Schema construction is isolated and compiled without optimization because
## it is non-hot code and large JSON literals otherwise dominate build time.

{.push optimization: none.}

import std/json
import ../slimmc_types
import ../../../common/[results_types, storage_schema]

proc schemaRecords*(m: Model): seq[JsonNode] {.codegenDecl: "__attribute__((optimize(\"O0\"))) $# $#$#".} =
  result.add storageRecord([("record_type", %("schema_header")), ("schema_name", %(StorageName)), ("schema_version", %(StorageFormatVersion)), ("byte_order", %("little")), ("column_storage", %("npy")), ("npy_format_version", %("1.0"))])
  for table in ["snapshots", "state", "chains", "sequences", "moments", "channel_events", "actions", "action_conditions", "feed_events", "monomer_balance", "species_balance", "kinetic_parameters/sets", "kinetic_parameters/values", "diagnostics/channel_trace"]:
    result.add storageRecord([("record_type", %("table")), ("name", %(table)), ("required", %(true))])
  template col(t,n,f,d,u: string) =
    result.add storageRecord([("record_type", %("column")), ("table", %(t)), ("name", %(n)), ("file", %(f)), ("dtype", %(d)), ("unit", %(u)), ("required", %(true))])
  col("snapshots","snapshot_id","snapshot_id.npy","uint64","")
  col("snapshots","time","time.npy","float64","s")
  col("snapshots","kmc_event","kmc_event.npy","uint64","")
  col("snapshots","snapshot_reason_id","snapshot_reason_id.npy","uint32","")
  col("snapshots","is_final","is_final.npy","bool","")
  col("snapshots","has_chains","has_chains.npy","bool","")
  col("snapshots","has_sequences","has_sequences.npy","bool","")
  col("snapshots","kinetic_parameter_set_id","kinetic_parameter_set_id.npy","uint64","")
  col("snapshots","volume_mL","volume_mL.npy","float64","mL")
  col("snapshots","kmc_volume_L","kmc_volume_L.npy","float64","L")
  col("snapshots","chain_count_live","chain_count_live.npy","uint64","")
  col("snapshots","chain_count_dead","chain_count_dead.npy","uint64","")
  col("snapshots","chain_count_total","chain_count_total.npy","uint64","")
  col("state","snapshot_id","snapshot_id.npy","uint64","")
  col("state","entity_id","entity_id.npy","uint32","")
  col("state","count","count.npy","uint64","")
  col("state","moles","moles.npy","float64","mol")
  col("state","concentration","concentration.npy","float64","mol/L")
  col("chains","chain_record_id","chain_record_id.npy","uint64","")
  col("chains","snapshot_id","snapshot_id.npy","uint64","")
  col("chains","population_id","population_id.npy","uint32","")
  col("chains","pool_id","pool_id.npy","uint32","")
  col("chains","origin_id","origin_id.npy","uint32","")
  col("chains","dp","dp.npy","uint64","")
  col("chains","molar_mass","molar_mass.npy","float64","g/mol")
  col("chains","count","count.npy","uint64","")
  col("chains","moles","moles.npy","float64","mol")
  col("chains","concentration","concentration.npy","float64","mol/L")
  col("chains","left_end_id","left_end_id.npy","uint32","")
  col("chains","right_end_id","right_end_id.npy","uint32","")
  col("chains","sequence_offset","sequence_offset.npy","uint64","")
  col("chains","sequence_length","sequence_length.npy","uint64","")
  col("sequences","symbols","symbols.npy","uint32","")
  col("moments","snapshot_id","snapshot_id.npy","uint64","")
  col("moments","population_scope_id","population_scope_id.npy","uint32","")
  col("moments","mass_basis_id","mass_basis_id.npy","uint32","")
  col("moments","chain_count","chain_count.npy","uint64","")
  col("moments","sum_dp","sum_dp.npy","float64","")
  col("moments","sum_dp2","sum_dp2.npy","float64","")
  col("moments","dp_n","dp_n.npy","float64","")
  col("moments","dp_w","dp_w.npy","float64","")
  col("moments","sum_molar_mass","sum_molar_mass.npy","float64","g/mol")
  col("moments","sum_molar_mass2","sum_molar_mass2.npy","float64","(g/mol)^2")
  col("moments","sum_molar_mass3","sum_molar_mass3.npy","float64","(g/mol)^3")
  col("moments","mn","mn.npy","float64","g/mol")
  col("moments","mw","mw.npy","float64","g/mol")
  col("moments","mz","mz.npy","float64","g/mol")
  col("moments","dispersity","dispersity.npy","float64","")
  col("channel_events","snapshot_id","snapshot_id.npy","uint64","")
  col("channel_events","channel_id","channel_id.npy","uint32","")
  col("channel_events","event_count","event_count.npy","uint64","")
  col("channel_events","productive_event_count","productive_event_count.npy","uint64","")
  col("channel_events","nonproductive_event_count","nonproductive_event_count.npy","uint64","")
  col("diagnostics/channel_trace","kmc_event","kmc_event.npy","uint64","")
  col("diagnostics/channel_trace","time","time.npy","float64","s")
  col("diagnostics/channel_trace","dt","dt.npy","float64","s")
  col("diagnostics/channel_trace","channel_id","channel_id.npy","uint32","")
  col("diagnostics/channel_trace","rate","rate.npy","float64","engine_native")
  col("diagnostics/channel_trace","propensity","propensity.npy","float64","s^-1")
  col("diagnostics/channel_trace","total_propensity","total_propensity.npy","float64","s^-1")
  col("kinetic_parameters/sets","kinetic_parameter_set_id","kinetic_parameter_set_id.npy","uint64","")
  col("kinetic_parameters/sets","start_kmc_event","start_kmc_event.npy","uint64","")
  col("kinetic_parameters/sets","start_time","start_time.npy","float64","s")
  col("kinetic_parameters/sets","has_source_action","has_source_action.npy","bool","")
  col("kinetic_parameters/sets","source_action_id","source_action_id.npy","uint64","")
  col("kinetic_parameters/values","kinetic_parameter_set_id","kinetic_parameter_set_id.npy","uint64","")
  col("kinetic_parameters/values","kinetic_parameter_id","kinetic_parameter_id.npy","uint32","")
  col("kinetic_parameters/values","value","value.npy","float64","")
  col("actions","action_id","action_id.npy","uint64","")
  col("actions","kmc_event","kmc_event.npy","uint64","")
  col("actions","time","time.npy","float64","s")
  col("actions","source_line","source_line.npy","uint32","")
  col("actions","action_type_id","action_type_id.npy","uint32","")
  col("actions","trigger_type_id","trigger_type_id.npy","uint32","")
  col("actions","scheduled_time","scheduled_time.npy","float64","s")
  col("actions","target_id","target_id.npy","uint32","")
  col("actions","requested_value","requested_value.npy","float64","")
  col("actions","before_value","before_value.npy","float64","")
  col("actions","after_value","after_value.npy","float64","")
  col("actions","state_changed","state_changed.npy","bool","")
  col("actions","output_written","output_written.npy","bool","")
  col("actions","has_snapshot","has_snapshot.npy","bool","")
  col("actions","snapshot_id","snapshot_id.npy","uint64","")
  col("actions","has_kinetic_parameter_set","has_kinetic_parameter_set.npy","bool","")
  col("actions","kinetic_parameter_set_id","kinetic_parameter_set_id.npy","uint64","")
  col("action_conditions","condition_record_id","condition_record_id.npy","uint64","")
  col("action_conditions","action_id","action_id.npy","uint64","")
  col("action_conditions","condition_index","condition_index.npy","uint32","")
  col("action_conditions","observable_id","observable_id.npy","uint32","")
  col("action_conditions","operator_id","operator_id.npy","uint32","")
  col("action_conditions","threshold","threshold.npy","float64","")
  col("action_conditions","observed_value","observed_value.npy","float64","")
  col("action_conditions","condition_met","condition_met.npy","bool","")
  col("feed_events","action_id","action_id.npy","uint64","")
  col("feed_events","feed_id","feed_id.npy","uint32","")
  col("feed_events","dose_mL","dose_mL.npy","float64","mL")
  col("feed_events","volume_before_mL","volume_before_mL.npy","float64","mL")
  col("feed_events","volume_after_mL","volume_after_mL.npy","float64","mL")
  col("feed_events","kmc_volume_before_L","kmc_volume_before_L.npy","float64","L")
  col("feed_events","kmc_volume_after_L","kmc_volume_after_L.npy","float64","L")
  col("monomer_balance","snapshot_id","snapshot_id.npy","uint64","")
  col("monomer_balance","monomer_id","monomer_id.npy","uint32","")
  col("monomer_balance","initial_moles","initial_moles.npy","float64","mol")
  col("monomer_balance","introduced_moles","introduced_moles.npy","float64","mol")
  col("monomer_balance","free_moles","free_moles.npy","float64","mol")
  col("monomer_balance","incorporated_moles","incorporated_moles.npy","float64","mol")
  col("monomer_balance","conversion","conversion.npy","float64","")
  col("species_balance","snapshot_id","snapshot_id.npy","uint64","")
  col("species_balance","entity_id","entity_id.npy","uint32","")
  col("species_balance","initial_moles","initial_moles.npy","float64","mol")
  col("species_balance","dosed_moles","dosed_moles.npy","float64","mol")
  col("species_balance","total_moles","total_moles.npy","float64","mol")
  col("species_balance","free_moles","free_moles.npy","float64","mol")
  col("species_balance","consumed_moles","consumed_moles.npy","float64","mol")
  for bid, sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("balance_entities")), ("id", %(bid)), ("name", %(sp.name)), ("kind", %((if sp.kind == skMonomer: "monomer" else: "species")))])
  for id, name in ["all", "live", "dead"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("population_scope")), ("id", %(id)), ("name", %(name))])
  for id, name in ["repeat_units", "with_end_groups"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("mass_bases")), ("id", %(id)), ("name", %(name))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("moments_reconstructable")), ("description", %("Every moments row is reconstructed from canonical chains rows for the same snapshot, population scope and mass basis."))])
  for id, name in ["print", "print_info", "save", "save_chains", "stop",
                   "set_k", "add_k", "set_temp", "add_temp", "set_c", "add_c", "feed", "print_memory"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_types")), ("id", %(id)), ("name", %(name))])
  for id, name in ["not_applicable", "at", "every", "when"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_trigger_types")), ("id", %(id)), ("name", %(name))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_observables")), ("id", %(0)), ("name", %("not_applicable"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_observables")), ("id", %(1)), ("name", %("conversion")), ("kind", %("conversion"))])
  for i, sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_observables")), ("id", %(2 + i)), ("name", %("concentration_" & sp.name)), ("kind", %("species_concentration")), ("species_id", %(i))])
  for id, name in ["not_applicable", "greater_than", "less_than"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_operators")), ("id", %(id)), ("name", %(name))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_targets")), ("id", %(0)), ("name", %("not_applicable"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_targets")), ("id", %(1)), ("name", %("temperature")), ("kind", %("temperature"))])
  var targetId = 2
  for rate in m.rates:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_targets")), ("id", %(targetId)), ("name", %(rate.name)), ("kind", %("rate_constant"))])
    inc targetId
  for sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("action_targets")), ("id", %(targetId)), ("name", %(sp.name)), ("kind", %("species"))])
    inc targetId
  for fid, feed in m.feeds:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("feeds")), ("id", %(fid)), ("name", %(feed.name))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("actions_optional_masks")), ("description", %("Absent optional foreign keys use has_*=false and identifier zero; absent numeric values are NaN."))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("action_conditions_and")), ("description", %("Rows with the same action_id are ordered by condition_index and form one AND conjunction. Scheduled actions have no condition rows."))])

  var kineticParameterId = 0
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kineticParameterId)), ("name", %("temperature")), ("kind", %("temperature")), ("unit", %("K"))])
  inc kineticParameterId
  for rate in m.rates:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kineticParameterId)), ("name", %(rate.name)), ("kind", %("rate_constant")), ("unit", %("engine_native"))])
    inc kineticParameterId
    if rate.declaredArrhenius:
      result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kineticParameterId)), ("name", %(rate.name & "__A")), ("kind", %("arrhenius_A")), ("unit", %("engine_native"))])
      inc kineticParameterId
      result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("kinetic_parameter_definitions")), ("id", %(kineticParameterId)), ("name", %(rate.name & "__Ea")), ("kind", %("arrhenius_Ea")), ("unit", %("J/mol"))])
      inc kineticParameterId
  result.add storageRecord([("record_type", %("rule")), ("name", %("kinetic_parameter_sets_complete")), ("description", %("Every set contains exactly one finite value for every kinetic_parameter_id."))])

  for cid, ch in m.channels:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("channels")), ("id", %(cid)), ("name", %(channelProgressIds(m)[cid])), ("kind", %(channelKindName(ch.kind))), ("description", %(ch.expr))])
  for id, name in ["initial", "scheduled", "action", "manual", "final"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("snapshot_reasons")), ("id", %(id)), ("name", %(name))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_populations")), ("id", %(0)), ("name", %("live"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_populations")), ("id", %(1)), ("name", %("dead"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_pools")), ("id", %(0)), ("name", %("not_applicable"))])
  for i, pool in m.pools:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_pools")), ("id", %(i + 1)), ("name", %(pool.name)), ("kind", %((if pool.kind == pkActive: "active" else: "dead")))])
  for id, name in ["unknown", "init", "transfer_m", "term_c", "term_d_H",
                   "term_d_U", "term_x", "transfer_h"]:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_origins")), ("id", %(id)), ("name", %(name))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_end_types")), ("id", %(0)), ("name", %("not_applicable"))])
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_end_types")), ("id", %(1)), ("name", %("unknown"))])
  for i, name in m.egNames:
    var rec = storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("chain_end_types")), ("id", %(i + 2)), ("name", %(name))])
    if i < m.egMassKnown.len and m.egMassKnown[i]:
      rec["molar_mass_contribution"] = %m.egMass[i]
    result.add rec
  result.add storageRecord([("record_type", %("rule")), ("name", %("chains_compressed_unique")), ("description", %("Within one snapshot, identical population/pool/origin/DP/end records are merged by count."))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("chains_snapshot_blocks")), ("description", %("Chain records of one snapshot form one contiguous deterministically ordered block."))])
  var eid = 0
  var monomerId = 0
  for sp in m.species:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %(sp.name)), ("kind", %((if sp.kind == skMonomer: "monomer" else: "species")))])
    if sp.kind == skMonomer:
      var rec = storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("monomers")), ("id", %(monomerId)), ("name", %(sp.name))])
      if sp.hasMw:
        rec["molar_mass_increment"] = %sp.mw
      result.add rec
      inc monomerId
    inc eid
  for pool in m.pools:
    result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %(pool.name)), ("kind", %("chain_pool"))])
    inc eid
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %("live_chains")), ("kind", %("chain_population"))]); inc eid
  result.add storageRecord([("record_type", %("dictionary_entry")), ("dictionary", %("state_entities")), ("id", %(eid)), ("name", %("dead_chains")), ("kind", %("chain_population"))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("state_dense")), ("description", %("Each snapshot contains all state entities in ascending entity_id order."))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("channel_events_dense")), ("description", %("Each snapshot contains all channels in ascending channel_id order."))])
  result.add storageRecord([("record_type", %("rule")), ("name", %("channel_events_match_kmc_event")), ("description", %("At each snapshot, kmc_event equals the sum of cumulative event_count over all channels."))])


{.pop.}
