import unittest, os, strutils
import slimmc_types
import slimmc_parser

proc writeTempModel(name, text: string): string =
  result = getTempDir() / name
  writeFile(result, text.strip() & "\n")

let base = """
desc "parser test"
param kmc_volume 1.0e-18
param temperature 298.15
param t_end 0.01
param max_steps 10
param seed 1
monomer M 1.0 100.0
species R 0.0
species Q 1.0
polymer P active
polymer D dead
rate ki 1.0
rate kp 1.0
rate ktr 1.0
"""

suite "parser validation v2.8.5":
  test "species MW is rejected with ValueError":
    let path = writeTempModel("slimmc_parser_species_mw.model", base.replace("species Q 1.0", "species Q 1.0 123.4") & "rxn Q -> R ktr")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)

  test "generic transfer alias is accepted":
    let path = writeTempModel("slimmc_parser_transfer_alias.model", base & """
endgroup R 68.0
macro init R + M -> P ki
macro prop P + M -> P kp
macro transfer P + Q -> D + R ktr
""")
    let m = parseModel(path)
    check m.channels.len == 3
    removeFile(path)

  test "unknown statement raises ValueError instead of quitting":
    let path = writeTempModel("slimmc_parser_unknown_statement.model", base & "nonsense Q -> R ktr")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)

  test "single and double quotes both work for desc (item 36)":
    check tokenize("desc \"double quoted\"", 1) == @["desc", "double quoted"]
    check tokenize("desc 'single quoted'", 1) == @["desc", "single quoted"]

  test "escape sequences \\n \\r \\t \\\\ \\\" \\' all work inside quoted strings (item 38)":
    let toks = tokenize("desc \"a\\nb\\rc\\td\\\\e\\\"f\\'g\"", 1)
    check toks.len == 2
    check toks[1] == "a\nb\rc\td\\e\"f\'g"

  test "A+B->C (no spaces) parses to the same channel count as A + B -> C (item 33/34)":
    let spaced = base & """
endgroup R 68.0
macro init R + M -> P ki
macro prop P + M -> P kp
"""
    let nospace = spaced.replace(" + ", "+").replace(" -> ", "->")
    let pathSpaced = writeTempModel("slimmc_parser_ws_spaced.model", spaced)
    let pathNospace = writeTempModel("slimmc_parser_ws_nospace.model", nospace)
    let mSpaced = parseModel(pathSpaced)
    let mNospace = parseModel(pathNospace)
    check mSpaced.channels.len == mNospace.channels.len
    check mSpaced.channels.len == 2
    removeFile(pathSpaced)
    removeFile(pathNospace)

  test "explicit positive exponent (1.0e+4) still tokenizes as one numeric token":
    check tokenize("rate ki const 2.0e+4", 1) == @["rate", "ki", "const", "2.0e+4"]

  test "keywords stay case-sensitive (item 39): 'Save' is a valid species name":
    let path = writeTempModel("slimmc_parser_case.model", base.replace("species Q 1.0", "species Save 1.0") & "rxn Save -> R ktr")
    let m = parseModel(path)
    check m.channels.len == 1
    removeFile(path)

  test "at_memory rejects unknown tokens instead of silently ignoring them (item 52)":
    let path = writeTempModel("slimmc_parser_atmemory_unknown.model",
      base & "at_memory 10MB bogus_token\nrxn Q -> R ktr\n")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)

  test "at_memory no longer accepts inert compact_dead/drop_dead_seq (item 53)":
    for tok in ["compact_dead", "drop_dead_seq"]:
      let path = writeTempModel("slimmc_parser_atmemory_inert.model",
        base & "at_memory 10MB " & tok & "\nrxn Q -> R ktr\n")
      expect(ValueError):
        discard parseModel(path)
      removeFile(path)

  test "bare print (no message) is a parser error; print_info is the canonical progress action (items 56/58)":
    let pathBarePrint = writeTempModel("slimmc_parser_bare_print.model", base & "rxn Q -> R ktr\nat 0.0 print\n")
    expect(ValueError):
      discard parseModel(pathBarePrint)
    removeFile(pathBarePrint)

    let pathPrintInfo = writeTempModel("slimmc_parser_print_info.model", base & "rxn Q -> R ktr\nat 0.0 print_info\n")
    let m = parseModel(pathPrintInfo)
    check m.scheduledActions.len == 1
    check m.scheduledActions[0].action == eaPrintInfo
    removeFile(pathPrintInfo)

  test "print with a message still works and requires exactly one argument (item 57)":
    let pathOneMsg = writeTempModel("slimmc_parser_print_msg.model", base & "rxn Q -> R ktr\nat 0.0 print \"hello\"\n")
    let m = parseModel(pathOneMsg)
    check m.scheduledActions[0].action == eaPrint
    check m.scheduledActions[0].args == @["hello"]
    removeFile(pathOneMsg)

  test "print_info takes no arguments":
    let path = writeTempModel("slimmc_parser_print_info_args.model", base & "rxn Q -> R ktr\nat 0.0 print_info \"oops\"\n")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)

  test "output_dir must be quoted and use identifier-like path segments":
    let unquotedPath = writeTempModel("slimmc_parser_outputdir_unquoted.model",
      base.replace("param seed 1", "param seed 1\nparam output_dir results/run_01") & "rxn Q -> R ktr\n")
    expect(ValueError):
      discard parseModel(unquotedPath)
    removeFile(unquotedPath)

    let invalidPath = writeTempModel("slimmc_parser_outputdir_invalid.model",
      base.replace("param seed 1", "param seed 1\nparam output_dir \"results/run-01\"") & "rxn Q -> R ktr\n")
    expect(ValueError):
      discard parseModel(invalidPath)
    removeFile(invalidPath)

    let quotedPath = writeTempModel("slimmc_parser_outputdir_quoted.model",
      base.replace("param seed 1", "param seed 1\nparam output_dir \"results/run_01\"") & "rxn Q -> R ktr\n")
    let m = parseModel(quotedPath)
    check m.outputDir.endsWith("results/run_01")
    removeFile(quotedPath)

  test "at_memory uses save or stop and rejects snapshot":
    let validPath = writeTempModel("slimmc_parser_atmemory_save.model",
      base & "at_memory 10MB save stop\nrxn Q -> R ktr\n")
    let m = parseModel(validPath)
    check m.memoryPolicy.snapshotOnLimit
    check m.memoryPolicy.stopOnLimit
    removeFile(validPath)

    let oldPath = writeTempModel("slimmc_parser_atmemory_snapshot.model",
      base & "at_memory 10MB snapshot\nrxn Q -> R ktr\n")
    expect(ValueError):
      discard parseModel(oldPath)
    removeFile(oldPath)


  test "multiple var declarations resolve independent comparison axes":
    let path = writeTempModel("slimmc_parser_multiple_var.model", base & """
var monomer M mol_L
var param temperature K
""")
    let m = parseModel(path)
    check m.variables.len == 2
    check m.variables[0].name == "M"
    check m.variables[0].value == 1.0
    check m.variables[1].name == "temperature"
    check m.variables[1].value == 298.15
    removeFile(path)

  test "duplicate var target name is rejected":
    let path = writeTempModel("slimmc_parser_duplicate_var.model", base & """
var monomer M mol_L
var param M K
""")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)

  test "leading underscore is accepted by the shared identifier contract":
    let path = writeTempModel("slimmc_parser_leading_underscore.model",
      base.replace("species Q 1.0", "species _Q 1.0") & "rxn _Q -> R ktr")
    let m = parseModel(path)
    check m.channels.len == 1
    removeFile(path)

  test "core model names are globally unique ignoring ASCII case":
    let path = writeTempModel("slimmc_parser_case_collision.model",
      base & "rate R const 1.0\n")
    expect(ValueError):
      discard parseModel(path)
    removeFile(path)
