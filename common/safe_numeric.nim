# Shared checked numeric operations for the homo and copo Nim engines.

import math, strutils

const MaxExactIntegerFloat64* = 9007199254740992.0 # 2^53

proc requireFinite*(x: float64; what: string): float64 =
  if x.classify in {fcNan, fcInf, fcNegInf}:
    raise newException(ValueError, what & " must be finite")
  result = x

proc requirePositiveFinite*(x: float64; what: string): float64 =
  discard requireFinite(x, what)
  if x <= 0.0:
    raise newException(ValueError, what & " must be > 0")
  result = x

proc requireNonNegativeFinite*(x: float64; what: string): float64 =
  discard requireFinite(x, what)
  if x < 0.0:
    raise newException(ValueError, what & " must be >= 0")
  result = x

proc checkedParseInt64*(s0: string; what: string = "integer"): int64 =
  let s = s0.strip()
  if s.len == 0:
    raise newException(ValueError, "empty " & what)
  var i = 0
  var negative = false
  if s[i] == '+':
    inc i
  elif s[i] == '-':
    negative = true
    inc i
  if i >= s.len or not s[i].isDigit:
    raise newException(ValueError, "invalid " & what & ": " & s0)
  var magnitude: uint64 = 0
  let limit = if negative: uint64(high(int64)) + 1'u64 else: uint64(high(int64))
  while i < s.len and s[i].isDigit:
    let digit = uint64(ord(s[i]) - ord('0'))
    if magnitude > (limit - digit) div 10'u64:
      raise newException(ValueError, what & " exceeds int64 range: " & s0)
    magnitude = magnitude * 10'u64 + digit
    inc i
  if i != s.len:
    raise newException(ValueError, "invalid trailing characters in " & what & ": " & s0)
  if negative:
    if magnitude == uint64(high(int64)) + 1'u64:
      result = low(int64)
    else:
      result = -int64(magnitude)
  else:
    result = int64(magnitude)

proc checkedAddInt64*(a, b: int64; what: string): int64 =
  if (b > 0 and a > high(int64) - b) or (b < 0 and a < low(int64) - b):
    raise newException(ValueError, what & " exceeds int64 range")
  result = a + b

proc checkedSubInt64*(a, b: int64; what: string): int64 =
  if b == low(int64):
    if a >= 0: raise newException(ValueError, what & " exceeds int64 range")
    return a - b
  result = checkedAddInt64(a, -b, what)

proc checkedAddInt32*(a, b: int32; what: string): int32 =
  let v = int64(a) + int64(b)
  if v < int64(low(int32)) or v > int64(high(int32)):
    raise newException(ValueError, what & " exceeds int32 range")
  result = int32(v)

proc checkedCountFromConc*(c, avogadro, volume: float64; what: string; allowNegative = false): int64 =
  discard requireFinite(c, what & " concentration")
  discard requirePositiveFinite(volume, what & " volume")
  discard requirePositiveFinite(avogadro, "Avogadro constant")
  if not allowNegative and c < 0.0:
    raise newException(ValueError, what & " concentration must be >= 0")
  let exact = c * avogadro * volume
  if exact.classify in {fcNan, fcInf, fcNegInf}:
    raise newException(ValueError, what & " molecule count is non-finite")
  if exact < float64(low(int64)) or exact > float64(high(int64)):
    raise newException(ValueError, what & " molecule count exceeds int64 range")
  result = int64(round(exact))

proc checkedSquareAsFloat*(x: int64; what: string): float64 =
  let d = float64(x)
  result = d * d
  if result.classify in {fcNan, fcInf, fcNegInf}:
    raise newException(ValueError, what & " square is non-finite")
