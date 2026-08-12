import unittest, os
import copo_types
import copo_parser

suite "parser whitespace and operator-spacing equivalence (items 33/34)":
  test "A+B->C (no spaces) parses to the same channels as A + B -> C":
    let spaced = """
desc "spaced"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
param max_steps 10
monomer A 0.20 100.12
monomer B 0.20 128.17
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 2.0e4
rate ki_b 2.0e4
rate kp_aa 300.0
rate kp_ab 120.0
rate kp_ba 250.0
rate kp_bb 80.0
rate ktc 1.0e7
rate ktd 1.0e7
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let nospace = spaced.replace(" + ", "+").replace(" -> ", "->")
    let pathSpaced = getTempDir() / "slimmc_copo_ws_spaced.model"
    let pathNospace = getTempDir() / "slimmc_copo_ws_nospace.model"
    writeFile(pathSpaced, spaced)
    writeFile(pathNospace, nospace)
    let mSpaced = parseModel(pathSpaced)
    let mNospace = parseModel(pathNospace)
    check mSpaced.channels.len == mNospace.channels.len
    check mSpaced.channels.len == 12
    for i in 0 ..< mSpaced.channels.len:
      check mSpaced.channels[i].name == mNospace.channels[i].name
      check mSpaced.channels[i].kind == mNospace.channels[i].kind

  test "tab-separated tokens parse identically to space-separated":
    let spaced = "macro init R + A -> PA ki_a"
    let tabbed = "macro\tinit\tR\t+\tA\t->\tPA\tki_a"
    check tokenize(spaced, 1) == tokenize(tabbed, 1)

  test "explicit positive exponent (1.0e+4) still tokenizes as one numeric token, not split on '+'":
    let toks = tokenize("rate ki_a 2.0e+4", 1)
    check toks == @["rate", "ki_a", "2.0e+4"]

  test "'+' between identifiers still splits into its own token even without spaces":
    let toks = tokenize("R+A->PA", 1)
    check toks == @["R", "+", "A", "->", "PA"]

  test "keywords stay case-sensitive (item 39): 'Save' is a valid species name, exact-case 'save' is still reserved":
    let ok = """
desc "case ok"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species Save 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
macro init Save + A -> PA ki_a
macro init Save + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let bad = ok.replace("species Save", "species save").replace("Save +", "save +")
    let pathOk = getTempDir() / "slimmc_copo_case_ok.model"
    let pathBad = getTempDir() / "slimmc_copo_case_bad.model"
    writeFile(pathOk, ok)
    writeFile(pathBad, bad)
    discard parseModel(pathOk)  # must not raise
    expect ValueError:
      discard parseModel(pathBad)

  test "single and double quotes both work for desc (item 36)":
    check tokenize("desc \"double quoted\"", 1) == @["desc", "double quoted"]
    check tokenize("desc 'single quoted'", 1) == @["desc", "single quoted"]

  test "escape sequences \\n \\r \\t \\\\ \\\" \\' all work inside quoted strings (item 38)":
    let toks = tokenize("desc \"a\\nb\\rc\\td\\\\e\\\"f\\'g\"", 1)
    check toks.len == 2
    check toks[1] == "a\nb\rc\td\\e\"f\'g"

  test "# inside a single-quoted string is not treated as a comment start":
    let toks = tokenize("desc 'has a # inside'", 1)
    check toks == @["desc", "has a # inside"]

  test "default max_steps=10_000_000_000 and when_check_events=1 (item 45, matches homo)":
    let text = """
desc "defaults check"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let path = getTempDir() / "slimmc_copo_defaults.model"
    writeFile(path, text)
    let m = parseModel(path)
    check m.max_steps == 10_000_000_000
    check m.whenCheckEvents == 1

  test "mass_model accepts only canonical values, case-sensitive, no legacy aliases (item 51)":
    check parseMassModelValue("repeat_units", 1) == mmRepeatUnits
    check parseMassModelValue("with_end_groups", 1) == mmWithEndgroups
    for bad in ["repeatunits", "repeat-unit", "RepeatUnits", "with_endgroups", "endgroups", "end_groups", "WITH_END_GROUPS"]:
      expect ValueError:
        discard parseMassModelValue(bad, 1)

  test "monomer/species/polymer reject excess tokens instead of silently ignoring them (item 48)":
    proc buildModel(monomerAExtra, speciesExtra: string): string =
      "desc \"excess tokens\"\nvar rate kp_ab L_mol_s\nparam kmc_volume 1.0e-18\nparam t_end 0.0\n" &
      "monomer A 0.2 100.0" & monomerAExtra & "\nmonomer B 0.2 128.0\n" &
      "species R 1.0e-4" & speciesExtra & "\n" &
      "polymer PA active\npolymer PB active\npolymer D dead\n" &
      "rate ki_a 1.0\nrate ki_b 1.0\nrate kp_aa 1.0\nrate kp_ab 1.0\nrate kp_ba 1.0\nrate kp_bb 1.0\nrate ktc 1.0\nrate ktd 1.0\n" &
      "macro init R + A -> PA ki_a\nmacro init R + B -> PB ki_b\n" &
      "macro prop PA + A -> PA kp_aa\nmacro prop PA + B -> PB kp_ab\nmacro prop PB + A -> PA kp_ba\nmacro prop PB + B -> PB kp_bb\n" &
      "macro term_c PA + PA -> D ktc\nmacro term_c PA + PB -> D ktc\nmacro term_c PB + PB -> D ktc\n" &
      "macro term_d PA + PA -> D ktd\nmacro term_d PA + PB -> D ktd\nmacro term_d PB + PB -> D ktd\n"

    # sanity: the well-formed model must parse fine
    let goodPath = getTempDir() / "slimmc_copo_excess_good.model"
    writeFile(goodPath, buildModel("", ""))
    discard parseModel(goodPath)

    let badMonomerPath = getTempDir() / "slimmc_copo_excess_monomer.model"
    writeFile(badMonomerPath, buildModel(" extra", ""))
    expect ValueError:
      discard parseModel(badMonomerPath)

    let badSpeciesPath = getTempDir() / "slimmc_copo_excess_species.model"
    writeFile(badSpeciesPath, buildModel("", " extra"))
    expect ValueError:
      discard parseModel(badSpeciesPath)

  test "rate const/arr and polymer active/dead matching is case-sensitive (item 39)":
    check parseMassModelValue("repeat_units", 1) == mmRepeatUnits  # sanity: function still importable/working
    let bad = """
desc "case sensitivity"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA Active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let path = getTempDir() / "slimmc_copo_polymer_case.model"
    writeFile(path, bad)
    expect ValueError:
      discard parseModel(path)

  test "at_memory rejects unknown tokens instead of silently ignoring them (item 52)":
    let text = """
desc "at_memory unknown"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
at_memory 10MB bogus_token
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let path = getTempDir() / "slimmc_copo_atmemory_unknown.model"
    writeFile(path, text)
    expect ValueError:
      discard parseModel(path)

  test "at_memory no longer accepts inert compact_dead/drop_dead_seq (item 53)":
    let textTemplate = """
desc "at_memory inert"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
at_memory 10MB $1
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    for tok in ["compact_dead", "drop_dead_seq"]:
      let path = getTempDir() / "slimmc_copo_atmemory_inert.model"
      writeFile(path, textTemplate.replace("$1", tok))
      expect ValueError:
        discard parseModel(path)

  test "'snapshot'/'chains'/'memory' action aliases are no longer accepted (item 55 -- canonical spellings only)":
    for badAction in ["snapshot", "chains", "memory"]:
      expect ValueError:
        discard parseActionKind(badAction, 1)
    # canonical spellings still work
    check parseActionKind("save", 1) == eaSave
    check parseActionKind("save_chains", 1) == eaSaveChains
    check parseActionKind("print_memory", 1) == eaPrintMemory

  test "action names are case-sensitive (item 39): 'SAVE'/'Save' are rejected, only exact-case 'save' works":
    check parseActionKind("save", 1) == eaSave
    for badCase in ["SAVE", "Save", "sAVE"]:
      expect ValueError:
        discard parseActionKind(badCase, 1)

  test "transfer_h is no longer accepted as a macro keyword (item 55 -- 'transfer' is canonical)":
    let bad = """
desc "transfer_h removed"
var rate kp_ab L_mol_s
param kmc_volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
species CTA 1.0e-4
species Rcta 0.0
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktr_a 1.0
rate ktc 1.0
rate ktd 1.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro transfer_h PA + CTA -> D + Rcta ktr_a
macro term_c PA + PA -> D ktc
macro term_d PA + PA -> D + D ktd
"""
    let path = getTempDir() / "slimmc_copo_transfer_h_removed.model"
    writeFile(path, bad)
    expect ValueError:
      discard parseModel(path)

  test "volume and t_end are mandatory (item 44) -- missing either raises, matching homo's error text":
    proc buildModel(includeVolume, includeTEnd: bool): string =
      result = "desc \"mandatory params\"\nvar rate kp_ab L_mol_s\n"
      if includeVolume: result &= "param kmc_volume 1.0e-18\n"
      if includeTEnd: result &= "param t_end 1.0\n"
      result &= "monomer A 0.2 100.0\nmonomer B 0.2 128.0\nspecies R 1.0e-4\n"
      result &= "polymer PA active\npolymer PB active\npolymer D dead\n"
      result &= "rate ki_a 1.0\nrate ki_b 1.0\nrate kp_aa 1.0\nrate kp_ab 1.0\nrate kp_ba 1.0\nrate kp_bb 1.0\nrate ktc 1.0\nrate ktd 1.0\n"
      result &= "macro init R + A -> PA ki_a\nmacro init R + B -> PB ki_b\n"
      result &= "macro prop PA + A -> PA kp_aa\nmacro prop PA + B -> PB kp_ab\nmacro prop PB + A -> PA kp_ba\nmacro prop PB + B -> PB kp_bb\n"
      result &= "macro term_c PA + PA -> D ktc\nmacro term_c PA + PB -> D ktc\nmacro term_c PB + PB -> D ktc\n"
      result &= "macro term_d PA + PA -> D ktd\nmacro term_d PA + PB -> D ktd\nmacro term_d PB + PB -> D ktd\n"

    let goodPath = getTempDir() / "slimmc_copo_mandatory_good.model"
    writeFile(goodPath, buildModel(true, true))
    discard parseModel(goodPath)  # must not raise

    let noVolumePath = getTempDir() / "slimmc_copo_mandatory_novolume.model"
    writeFile(noVolumePath, buildModel(false, true))
    expect ValueError:
      discard parseModel(noVolumePath)

    let noTEndPath = getTempDir() / "slimmc_copo_mandatory_notend.model"
    writeFile(noTEndPath, buildModel(true, false))
    expect ValueError:
      discard parseModel(noTEndPath)

  test "bare print (no message) is a parser error; print_info is canonical progress action (items 56/58)":
    proc buildModel(actionLine: string): string =
      "desc \"print_info test\"\nvar rate kp_ab L_mol_s\nparam kmc_volume 1.0e-18\nparam t_end 0.0\n" &
      "monomer A 0.2 100.0\nmonomer B 0.2 128.0\nspecies R 1.0e-4\n" &
      "polymer PA active\npolymer PB active\npolymer D dead\n" &
      "rate ki_a 1.0\nrate ki_b 1.0\nrate kp_aa 1.0\nrate kp_ab 1.0\nrate kp_ba 1.0\nrate kp_bb 1.0\nrate ktc 1.0\nrate ktd 1.0\n" &
      "macro init R + A -> PA ki_a\nmacro init R + B -> PB ki_b\n" &
      "macro prop PA + A -> PA kp_aa\nmacro prop PA + B -> PB kp_ab\nmacro prop PB + A -> PA kp_ba\nmacro prop PB + B -> PB kp_bb\n" &
      "macro term_c PA + PA -> D ktc\nmacro term_c PA + PB -> D ktc\nmacro term_c PB + PB -> D ktc\n" &
      "macro term_d PA + PA -> D ktd\nmacro term_d PA + PB -> D ktd\nmacro term_d PB + PB -> D ktd\n" &
      actionLine & "\n"

    let barePath = getTempDir() / "slimmc_copo_bare_print.model"
    writeFile(barePath, buildModel("at 0.0 print"))
    expect ValueError:
      discard parseModel(barePath)

    let infoPath = getTempDir() / "slimmc_copo_print_info.model"
    writeFile(infoPath, buildModel("at 0.0 print_info"))
    let m = parseModel(infoPath)
    check m.scheduledActions.len == 1
    check m.scheduledActions[0].action == eaPrintInfo

    let msgPath = getTempDir() / "slimmc_copo_print_msg.model"
    writeFile(msgPath, buildModel("at 0.0 print \"hello\""))
    let m2 = parseModel(msgPath)
    check m2.scheduledActions[0].action == eaPrint
    check m2.scheduledActions[0].args == @["hello"]

    let infoArgsPath = getTempDir() / "slimmc_copo_print_info_args.model"
    writeFile(infoArgsPath, buildModel("at 0.0 print_info \"oops\""))
    expect ValueError:
      discard parseModel(infoArgsPath)

  test "top-level statement keywords are case-sensitive (item 39): 'PARAM' is rejected":
    let bad = """
desc "keyword case"
var rate kp_ab L_mol_s
PARAM volume 1.0e-18
param t_end 0.0
monomer A 0.2 100.0
monomer B 0.2 128.0
species R 1.0e-4
polymer PA active
polymer PB active
polymer D dead
rate ki_a 1.0
rate ki_b 1.0
rate kp_aa 1.0
rate kp_ab 1.0
rate kp_ba 1.0
rate kp_bb 1.0
rate ktc 1.0
rate ktd 1.0
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc
macro term_d PA + PA -> D ktd
macro term_d PA + PB -> D ktd
macro term_d PB + PB -> D ktd
"""
    let path = getTempDir() / "slimmc_copo_keyword_case.model"
    writeFile(path, bad)
    expect ValueError:
      discard parseModel(path)

  test "output_dir must be quoted and use identifier-like path segments":
    proc buildModel(outputDirLine: string): string =
      "desc \"output_dir contract\"\nvar rate kp_ab L_mol_s\nparam kmc_volume 1.0e-18\nparam t_end 0.0\n" &
      outputDirLine & "\n" &
      "monomer A 0.2 100.0\nmonomer B 0.2 128.0\nspecies R 1.0e-4\n" &
      "polymer PA active\npolymer PB active\npolymer D dead\n" &
      "rate ki_a 1.0\nrate ki_b 1.0\nrate kp_aa 1.0\nrate kp_ab 1.0\nrate kp_ba 1.0\nrate kp_bb 1.0\nrate ktc 1.0\nrate ktd 1.0\n" &
      "macro init R + A -> PA ki_a\nmacro init R + B -> PB ki_b\n" &
      "macro prop PA + A -> PA kp_aa\nmacro prop PA + B -> PB kp_ab\nmacro prop PB + A -> PA kp_ba\nmacro prop PB + B -> PB kp_bb\n" &
      "macro term_c PA + PA -> D ktc\nmacro term_c PA + PB -> D ktc\nmacro term_c PB + PB -> D ktc\n" &
      "macro term_d PA + PA -> D ktd\nmacro term_d PA + PB -> D ktd\nmacro term_d PB + PB -> D ktd\n"

    let unquotedPath = getTempDir() / "slimmc_copo_outputdir_unquoted.model"
    writeFile(unquotedPath, buildModel("param output_dir results/run_01"))
    expect ValueError:
      discard parseModel(unquotedPath)

    let invalidPath = getTempDir() / "slimmc_copo_outputdir_invalid.model"
    writeFile(invalidPath, buildModel("param output_dir \"results/run-01\""))
    expect ValueError:
      discard parseModel(invalidPath)

    let quotedPath = getTempDir() / "slimmc_copo_outputdir_quoted.model"
    writeFile(quotedPath, buildModel("param output_dir \"results/run_01\""))
    let m = parseModel(quotedPath)
    check m.output_dir.endsWith("results/run_01")

  test "at_memory uses save or stop and rejects snapshot":
    proc buildAtMemory(action: string): string =
      "desc \"at_memory contract\"\nvar rate kp_ab L_mol_s\nparam kmc_volume 1.0e-18\nparam t_end 0.0\n" &
      "monomer A 0.2 100.0\nmonomer B 0.2 128.0\nspecies R 1.0e-4\n" &
      "polymer PA active\npolymer PB active\npolymer D dead\n" &
      "rate ki_a 1.0\nrate ki_b 1.0\nrate kp_aa 1.0\nrate kp_ab 1.0\nrate kp_ba 1.0\nrate kp_bb 1.0\nrate ktc 1.0\nrate ktd 1.0\n" &
      "at_memory 10MB " & action & "\n" &
      "macro init R + A -> PA ki_a\nmacro init R + B -> PB ki_b\n" &
      "macro prop PA + A -> PA kp_aa\nmacro prop PA + B -> PB kp_ab\nmacro prop PB + A -> PA kp_ba\nmacro prop PB + B -> PB kp_bb\n" &
      "macro term_c PA + PA -> D ktc\nmacro term_c PA + PB -> D ktc\nmacro term_c PB + PB -> D ktc\n" &
      "macro term_d PA + PA -> D ktd\nmacro term_d PA + PB -> D ktd\nmacro term_d PB + PB -> D ktd\n"

    let validPath = getTempDir() / "slimmc_copo_atmemory_save.model"
    writeFile(validPath, buildAtMemory("save stop"))
    let m = parseModel(validPath)
    check m.memoryPolicy.snapshotOnLimit
    check m.memoryPolicy.stopOnLimit

    let oldPath = getTempDir() / "slimmc_copo_atmemory_snapshot.model"
    writeFile(oldPath, buildAtMemory("snapshot"))
    expect ValueError:
      discard parseModel(oldPath)
