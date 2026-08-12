import ../../src/copo_stats
import ../../../common/safe_numeric

suite "checked numeric safety v0.3.0 CLI hardening":
  test "shared int64 parser rejects overflow":
    expect ValueError:
      discard checkedParseInt64("9223372036854775808", "integer parameter")
    check checkedParseInt64("-9223372036854775808") == low(int64)

  test "concentration conversion rejects non-finite and out-of-range populations":
    expect ValueError:
      discard countFromConc(Inf, 1.0)
    expect ValueError:
      discard countFromConc(1.0e300, 1.0)

  test "int32 chain arithmetic rejects overflow":
    expect ValueError:
      discard checkedAddInt32(high(int32), 1'i32, "chain DP")
