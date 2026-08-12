import unittest
import ../model_contract

suite "shared model contract":
  test "identifier grammar is shared and Python-compatible":
    for name in ["A", "_A", "kp_IA_BMA", "run_001"]:
      check isValidModelIdentifier(name)
    for name in ["", "1A", "A-B", "A.B", "A+B", "A B"]:
      check not isValidModelIdentifier(name)

  test "model paths are quoted separately but use identifier-like segments":
    for path in ["results", "results/run_01", "/net/project/results", "C:/results/run_01"]:
      check isValidModelPath(path)
    for path in ["", ".", "../results", "results/run-01", "results/run.01", "results/run 01", "results#1", "results//run_01"]:
      check not isValidModelPath(path)

  test "language keywords are case-sensitive":
    check isReservedModelIdentifier("save")
    check not isReservedModelIdentifier("Save")

  test "case-insensitive collision key is stable":
    check modelNameKey("kp_IA") == modelNameKey("KP_ia")

  test "common defaults are explicit and portable":
    check DefaultTemperatureK == 298.15
    check DefaultMaxSteps == 10_000_000_000'i64
    check DefaultWhenCheckEvents == 1'i64
    check DefaultSeed == 12345'i64
    check DefaultDpMax == int64(high(int32))
    check DefaultSequenceMode == "composition"
    check DefaultMassModel == "repeat_units"
