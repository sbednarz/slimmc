import unittest, os, random
import copo_types
import copo_parser
import copo_kmc

proc writeTempModel(name, body: string): string =
  result = getTempDir() / name
  writeFile(result, body)

suite "mass model v0.5":
  test "repeat_units remains the default even when endgroups are declared":
    let path = writeTempModel("slimmc_copo_mass_repeat_units.model", """
param kmc_volume 1.0e-18
param t_end 1.0
param max_steps 10
monomer A 0.0 100.0
monomer B 0.0 128.0
species R 0.0
endgroup R 10.0
endgroup ACTIVE 0.0
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
macro init R + A -> PA ki_a
""")
    var m = parseModel(path)
    var s = initState(m)
    s.speciesN[0] = 1
    s.monomerN[0] = 1
    var rng = initRand(501)
    applyChannel(m, s, 0, rng)
    check m.mass_model == mmRepeatUnits
    check s.livePools[0][0].mass == 100.0

  test "with_end_groups adds initiator and transfer endgroup masses":
    let path = writeTempModel("slimmc_copo_mass_transfer.model", """
param kmc_volume 1.0e-18
param t_end 1.0
param max_steps 10
param mass_model with_end_groups
monomer A 0.0 100.0
monomer B 0.0 128.0
species R 0.0
species CTA 0.0
species Rcta 0.0
endgroup R 10.0
endgroup ACTIVE 0.0
endgroup CTA 50.0
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate kp_ab 1.0
rate ktr 1.0
macro init R + A -> PA ki_a
macro prop PA + B -> PB kp_ab
macro transfer PB + CTA -> D + Rcta ktr
""")
    var m = parseModel(path)
    var s = initState(m)
    var rng = initRand(502)
    s.speciesN[0] = 1
    s.monomerN[0] = 1
    applyChannel(m, s, 0, rng)
    check s.livePools[0][0].mass == 110.0
    s.monomerN[1] = 1
    applyChannel(m, s, 1, rng)
    check s.livePools[1][0].mass == 238.0
    s.speciesN[1] = 1
    applyChannel(m, s, 2, rng)
    check s.deadChains[0].right_end == "CTA"
    check s.deadChains[0].mass == 288.0

  test "with_end_groups handles combination end masses":
    let path = writeTempModel("slimmc_copo_mass_term_c.model", """
param kmc_volume 1.0e-18
param t_end 1.0
param max_steps 10
param mass_model with_end_groups
monomer A 0.0 100.0
monomer B 0.0 128.0
species R1 0.0
species R2 0.0
endgroup R1 10.0
endgroup R2 20.0
endgroup ACTIVE 0.0
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate ktc 1.0
macro init R1 + A -> PA ki_a
macro init R2 + B -> PB ki_b
macro term_c PA + PB -> D ktc
""")
    var m = parseModel(path)
    var s = initState(m)
    var rng = initRand(503)
    s.speciesN[0] = 1
    s.monomerN[0] = 1
    applyChannel(m, s, 0, rng)
    s.speciesN[1] = 1
    s.monomerN[1] = 1
    applyChannel(m, s, 1, rng)
    applyChannel(m, s, 2, rng)
    check s.deadChains[0].left_end == "R1"
    check s.deadChains[0].right_end == "R2"
    check s.deadChains[0].mass == 258.0
