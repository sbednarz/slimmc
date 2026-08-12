import unittest, os, tables
import copo_types
import copo_parser

suite "parser v0.6 penultimate propagation syntax":
  test "parses penultimate active pools through ordinary macro prop channels":
    let path = getTempDir() / "slimmc_copo_parser_v06_penultimate.model"
    writeFile(path, """
desc "parser v0.6 penultimate propagation"
var rate kp_ab_a L_mol_s
param kmc_volume 1.0e-19
param t_end 1.0
param max_steps 10
param seed 1
monomer A 0.1 100.0
monomer B 0.1 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer PAA active
polymer PAB active
polymer PBA active
polymer PBB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa_a 10.0
rate kp_aa_b 20.0
rate kp_ab_a 30.0
rate kp_ab_b 40.0
rate kp_ba_a 50.0
rate kp_ba_b 60.0
rate kp_bb_a 70.0
rate kp_bb_b 80.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PAA kp_aa_a
macro prop PA + B -> PAB kp_aa_b
macro prop PB + A -> PBA kp_ba_a
macro prop PB + B -> PBB kp_ba_b
macro prop PAA + A -> PAA kp_aa_a
macro prop PAA + B -> PAB kp_aa_b
macro prop PAB + A -> PBA kp_ab_a
macro prop PAB + B -> PBB kp_ab_b
macro prop PBA + A -> PAA kp_ba_a
macro prop PBA + B -> PAB kp_ba_b
macro prop PBB + A -> PBA kp_bb_a
macro prop PBB + B -> PBB kp_bb_b
""")
    let m = parseModel(path)
    check m.description == "parser v0.6 penultimate propagation"
    check m.variables.len == 1
    check m.variables[0].kind == "rate"
    check m.variables[0].name == "kp_ab_a"
    check m.variables[0].value == 30.0
    check m.channels.len == 14

    let paa = m.poolByName["PAA"]
    let pab = m.poolByName["PAB"]
    let pba = m.poolByName["PBA"]
    let pbb = m.poolByName["PBB"]
    check m.poolTerminalMer[paa] == 0
    check m.poolPenultimateMer[paa] == 0
    check m.poolTerminalMer[pab] == 1
    check m.poolPenultimateMer[pab] == 0
    check m.poolTerminalMer[pba] == 0
    check m.poolPenultimateMer[pba] == 1
    check m.poolTerminalMer[pbb] == 1
    check m.poolPenultimateMer[pbb] == 1
