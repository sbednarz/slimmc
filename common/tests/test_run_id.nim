import std/[unittest]
import ../run_id

suite "run_id contract":
  test "accepts canonical identifiers":
    for value in ["run", "run_001", "_test", "DBI80_deprop", "A1"]:
      check isValidRunId(value)

  test "rejects non-canonical identifiers":
    for value in ["", "1run", "run-001", "run.001", "run+001", "run=001",
                  "run 001", "run?", "run*", "[run]", "ąrun"]:
      check not isValidRunId(value)

  test "derives and validates model stem":
    check validateModelRunId("/tmp/fig5a_exp100.model") == "fig5a_exp100"
    expect ValueError:
      discard validateModelRunId("/tmp/[fig5a]?.model")
