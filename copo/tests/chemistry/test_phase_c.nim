import unittest, os, strutils
import copo_types
import copo_parser

proc writeTempModelPhaseC(name, text: string): string =
  result = getTempDir() / name
  writeFile(result, text.strip() & "\n")

let base = """
desc "phase C parser test"
param kmc_volume 1.0e-18
param temperature 298.15
param t_end 1.0
param max_steps 1000
param when_check_events 5
param seed 123
monomer M 0.1 100.0
monomer B 0.0 128.0
species Q 0.2
species R 0.0
polymer D dead
rate k const 1.0
rxn Q -> R k
"""

suite "copo scheduling/actions/conditions phase C":
  test "H11 parses at, every-from-zero, and shifted from/step schedules":
    let path = writeTempModelPhaseC("copo_phase_c_schedule.model", base & """
at 0.0 save
at 0.25 set_k k 2.0
every 0.10 save_chains
from 0.30 step 0.20 save
""")
    let m = parseModel(path)
    check m.scheduledActions.len == 4
    check m.scheduledActions[0].nextTime == 0.0
    check not m.scheduledActions[0].repeat
    check m.scheduledActions[1].action == eaSetK
    check m.scheduledActions[1].args == @["k", "2.0"]
    check m.scheduledActions[2].repeat
    check m.scheduledActions[2].nextTime == 0.0
    check m.scheduledActions[2].period == 0.10
    check m.scheduledActions[3].repeat
    check m.scheduledActions[3].nextTime == 0.30
    check m.scheduledActions[3].period == 0.20
    removeFile(path)

  test "H11 rejects non-positive every period and scheduled stop":
    for line in ["every 0 save", "every -1 save", "from -1 step 1 save", "from 0 step 0 save", "from 0 save", "at 0.2 stop"]:
      let path = writeTempModelPhaseC("copo_phase_c_bad_schedule.model", base & line & "\n")
      expect(ValueError):
        discard parseModel(path)
      removeFile(path)

  test "H12 parses set/add rate and temperature actions":
    let path = writeTempModelPhaseC("copo_phase_c_rates.model", base & """
at 0.1 set_k k 2.5
at 0.2 add_k k -0.5
at 0.3 set_temp 320
at 0.4 add_temp 5
""")
    let m = parseModel(path)
    check m.scheduledActions.len == 4
    check m.scheduledActions[0].action == eaSetK
    check m.scheduledActions[1].action == eaAddK
    check m.scheduledActions[2].action == eaSetTemp
    check m.scheduledActions[3].action == eaAddTemp
    removeFile(path)

  test "H13 parses set_c and signed add_c":
    let path = writeTempModelPhaseC("copo_phase_c_conc.model", base & """
at 0.1 set_c Q 0.3
at 0.2 add_c Q 0.1
at 0.3 add_c Q -0.05
""")
    let m = parseModel(path)
    check m.scheduledActions.len == 3
    check m.scheduledActions[0].action == eaSetC
    check m.scheduledActions[1].action == eaAddC
    check m.scheduledActions[2].args == @["Q", "-0.05"]
    removeFile(path)

  test "H14 parses AND conditions and keeps separate when lines as independent actions":
    let path = writeTempModelPhaseC("copo_phase_c_when.model", base & """
when X > 0.1 and c Q < 0.15 set_k k 0
when c R > 0.01 print "R appeared"
when c R > 0.02 stop
""")
    let m = parseModel(path)
    check m.whenCheckEvents == 5
    check m.conditionalActions.len == 3
    check m.conditionalActions[0].conditions.len == 2
    check m.conditionalActions[0].action == eaSetK
    check m.conditionalActions[1].conditions.len == 1
    check m.conditionalActions[2].action == eaStop
    removeFile(path)

  test "H14 rejects unsupported comparison and incomplete condition":
    for line in ["when c Q >= 0.1 stop", "when c Q > stop", "when Q > 0.1 stop"]:
      let path = writeTempModelPhaseC("copo_phase_c_bad_when.model", base & line & "\n")
      expect(ValueError):
        discard parseModel(path)
      removeFile(path)


suite "copo simulation-time tolerance":
  test "mixed absolute and relative tolerance scales only at large times":
    check timeTolerance(1.0e-6, 1.0e-6) == Eps
    check timeTolerance(1.0e8, 1.0e8) > Eps
    check timeTolerance(1.0e8, 1.0e8) < 1.0e-5
    check timeClose(1.0e8, 1.0e8 + 1.0e-7)
    check not timeClose(1.0e8, 1.0e8 + 1.0e-4)
