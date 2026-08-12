import unittest, os, tables, strutils
import copo_types
import copo_parser

suite "parser v0.2 syntax":
  test "parses desc var transfer and deprop without a case directive":
    let path = getTempDir() / "slimmc_copo_parser_v02.model"
    writeFile(path, """
desc "parser v0.2 chemistry syntax"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-19
param t_end 1.0
param max_steps 10
param seed 1
param sequence_mode full
monomer A 0.1 100.0
monomer B 0.1 128.0
species R 1.0e-4
species CTA 1.0e-3
species Rcta 0.0
polymer PA active
polymer PB active
polymer D dead
rate kp_ab 120.0
rate ki_a 1.0
rate ki_b 1.0
rate ktr_a 1.0
rate kdeprop_ba 1.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro transfer PA + CTA -> D + Rcta ktr_a
macro deprop PA -> PB + A kdeprop_ba
""")
    let m = parseModel(path)
    check m.description == "parser v0.2 chemistry syntax"
    check m.variables.len == 1
    check m.variables[0].kind == "rate"
    check m.variables[0].name == "kp_ab"
    check m.variables[0].value == 120.0
    check m.variables[0].unit == "L_mol_s"
    check m.output_dir.endsWith("results" / "slimmc_copo_parser_v02")
    check m.poolTerminalMer[m.poolByName["PA"]] == 0
    check m.poolTerminalMer[m.poolByName["PB"]] == 1
    check m.channels[^2].kind == chMacroTransfer
    check m.channels[^1].kind == chMacroDeprop

  test "multiple var declarations resolve independent comparison axes":
    let path = getTempDir() / "slimmc_copo_parser_v02_two_var.model"
    writeFile(path, """
desc "two var axes"
var monomer A mol_L
var param temperature K
param kmc_volume 1.0e-19
param temperature 343.15
param t_end 1.0
monomer A 0.1 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
macro init R + A -> PA ki_a
""")
    let m = parseModel(path)
    check m.variables.len == 2
    check m.variables[0].name == "A"
    check m.variables[0].value == 0.1
    check m.variables[1].name == "temperature"
    check m.variables[1].value == 343.15

  test "duplicate var target name is rejected":
    let path = getTempDir() / "slimmc_copo_parser_v02_duplicate_var.model"
    writeFile(path, """
var monomer A mol_L
var param A K
param kmc_volume 1.0e-19
param t_end 1.0
monomer A 0.1 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
macro init R + A -> PA ki_a
""")
    expect(ValueError):
      discard parseModel(path)

  test "case directive is intentionally not part of model syntax":
    let path = getTempDir() / "slimmc_copo_parser_v02_case_directive.model"
    writeFile(path, """
case COP_SHOULD_NOT_PARSE
desc "case should fail"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-19
monomer A 0.1 100.0
monomer B 0.1 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
macro init R + A -> PA ki_a
""")
    expect(ValueError):
      discard parseModel(path)

  test "leading underscore is accepted by the shared identifier contract":
    let path = getTempDir() / "slimmc_copo_parser_v02_leading_underscore.model"
    writeFile(path, """
param kmc_volume 1.0e-19
param t_end 1.0
monomer A 0.1 100.0
monomer B 0.2 128.0
species _R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
macro init _R + A -> PA ki_a
""")
    let m = parseModel(path)
    check m.speciesByName.hasKey("_R")

  test "core model names are globally unique ignoring ASCII case":
    let path = getTempDir() / "slimmc_copo_parser_v02_case_collision.model"
    writeFile(path, """
param kmc_volume 1.0e-19
param t_end 1.0
monomer A 0.1 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate r 1.0
macro init R + A -> PA r
""")
    expect(ValueError):
      discard parseModel(path)
