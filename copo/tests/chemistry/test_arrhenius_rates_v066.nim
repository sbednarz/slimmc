import unittest, os, math, strutils, json, tables
import copo_types
import copo_parser
import copo_kmc
import storage/copo_schema

proc writeArrModel066(name, body: string): string =
  result = getTempDir() / name
  writeFile(result, body)

suite "Arrhenius rates v0.6.13":
  test "parser accepts rate NAME arr Apre Ea and evaluates at model temperature":
    let path = writeArrModel066("slimmc_copo_arr_rate_v066.model", """
param kmc_volume 1.0e-18
param temperature 300.0
param t_end 1.0
param max_steps 10
var rate kp_arr L_mol_s
monomer A 0.0 100.0
monomer B 0.0 128.0
species R 0.0
polymer PA active
polymer PB active
polymer D dead
rate ki const 1.0
rate kp_arr arr 1.0e5 1000.0
macro init R + A -> PA ki
macro prop PA + B -> PB kp_arr
""")
    let m = parseModel(path)
    let kid = m.rateByName["kp_arr"]
    check m.rates[kid].kind == rkArr
    check m.rates[kid].Apre == 1.0e5
    check m.rates[kid].Ea == 1000.0
    let expected = 1.0e5 * exp(-1000.0 / (Rgas * 300.0))
    check abs(m.rateValue(kid) - expected) < 1.0e-8 * expected
    check m.variables.len == 1
    check abs(m.variables[0].value - expected) < 1.0e-8 * expected

  test "Arrhenius rate value follows temperature in propensity calculation":
    var m = Model()
    m.V = 1.0
    m.T = 300.0
    m.monomers = @[MonomerDef(name: "A", c0: 0.0, mw: 100.0)]
    m.pools = @[PoolDef(name: "PA", kind: pkActive)]
    m.poolTerminalMer = @[0]
    m.poolPenultimateMer = @[-1]
    m.rates = @[RateDef(name: "karr", kind: rkArr, Apre: 1.0e5, Ea: 5000.0)]
    m.channels = @[KmcChannel(name: "deprop_PA_A", kind: chMacroDeprop, kId: 0, pool1: 0, poolOut: 0, monomerId: 0)]

    var s = State()
    s.livePools = newSeq[seq[LiveChain]](1)
    var c = LiveChain(id: 1, left_end: "R", right_end: "ACTIVE", nMer: @[int32(2)], dp: 2, mass: 200.0, last: 0, prev: 0, formedBy: fbInit)
    c.mers.data = @[uint8(0), uint8(0)]
    s.livePools[0].add c

    let p300 = computePropensities(m, s)[0]
    m.T = 600.0
    let p600 = computePropensities(m, s)[0]
    check p600 > p300
    check abs(p300 - (1.0e5 * exp(-5000.0 / (Rgas * 300.0)))) < 1.0e-8 * p300
    check abs(p600 - (1.0e5 * exp(-5000.0 / (Rgas * 600.0)))) < 1.0e-8 * p600

  test "Storage schema advertises Arrhenius parameter definitions":
    let path = writeArrModel066("slimmc_copo_arr_runinfo_v066.model", """
param kmc_volume 1.0e-18
param temperature 310.0
param t_end 1.0
param max_steps 10
monomer A 0.0 100.0
monomer B 0.0 128.0
species R 0.0
polymer PA active
polymer PB active
polymer D dead
rate ki 1.0
rate kp_arr arr 1.0e5 1000.0
macro init R + A -> PA ki
macro prop PA + B -> PB kp_arr
""")
    let m = parseModel(path)
    let txt = $schemaRecords(m)
    check txt.contains("kinetic_parameter_definitions")
    check txt.contains("kp_arr")
    check txt.contains("arrhenius_A")
    check txt.contains("arrhenius_Ea")
    check txt.contains("J/mol")

  test "Storage schema publishes all monomer repeat-unit masses":
    var m = Model()
    m.monomers = @[
      MonomerDef(name: "A", c0: 0.0, mw: 100.0),
      MonomerDef(name: "B", c0: 0.0, mw: 128.0)
    ]
    let txt = $schemaRecords(m)
    check txt.count("\"dictionary\":\"monomers\"") == 2
    check txt.contains("\"name\":\"A\",\"molar_mass_increment\":100.0")
    check txt.contains("\"name\":\"B\",\"molar_mass_increment\":128.0")

  test "negative Arrhenius pre-exponential factor is rejected":
    let path = writeArrModel066("slimmc_copo_bad_arr_rate_v066.model", """
param kmc_volume 1.0e-18
param t_end 1.0
param max_steps 10
monomer A 0.0 100.0
monomer B 0.0 128.0
species R 0.0
polymer PA active
polymer PB active
polymer D dead
rate bad arr -1.0 1000.0
""")
    expect ValueError:
      discard parseModel(path)
