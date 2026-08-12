import unittest, os, strutils
import copo_types
import copo_parser

proc baseModel(extraParams = ""; extraMacros = ""): string =
  result = """
desc "parser sequence-mode syntax"
param kmc_volume 1.0e-19
param t_end 1.0
param max_steps 10
param seed 1
""" & extraParams & """
monomer A 0.1 100.0
monomer B 0.1 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kd 0.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
""" & extraMacros

suite "parser sequence_mode and dp_max":
  test "parses dp_max as numeric var target":
    let path = getTempDir() / "slimmc_copo_parser_dpmax.model"
    writeFile(path, baseModel("""
param dp_max 4
param sequence_mode composition
var param dp_max DP
"""))
    let m = parseModel(path)
    check m.dp_max == 4
    check m.sequence_mode == "composition"
    check m.variables.len == 1
    check m.variables[0].name == "dp_max"
    check m.variables[0].value == 4.0
    check m.output_dir.endsWith("results" / "slimmc_copo_parser_dpmax")

  test "composition is accepted without depropagation":
    let path = getTempDir() / "slimmc_copo_parser_composition.model"
    writeFile(path, baseModel("param sequence_mode composition\n"))
    check parseModel(path).sequence_mode == "composition"

  test "full is accepted with depropagation":
    let path = getTempDir() / "slimmc_copo_parser_full_deprop.model"
    writeFile(path, baseModel("param sequence_mode full\n", "macro deprop PA -> PB + A kd\n"))
    check parseModel(path).sequence_mode == "full"

  test "composition is rejected with any declared depropagation channel":
    let path = getTempDir() / "slimmc_copo_parser_composition_deprop.model"
    writeFile(path, baseModel("param sequence_mode composition\n", "macro deprop PA -> PB + A kd\n"))
    expect ValueError:
      discard parseModel(path)
