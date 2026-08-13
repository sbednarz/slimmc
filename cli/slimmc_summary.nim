import std/[os, strutils, json, math]

const AppName = "slimmc-summary"
const AppVersion {.strdefine.} = "5.0.1"

type NpyData = object
  descr: string
  count: int
  offset: int
  bytes: string

proc fail(msg: string) =
  stderr.writeLine("ERROR: " & msg)
  quit(1)

proc parseNpy(path: string): NpyData =
  let b = readFile(path)
  if b.len < 10 or b[0..5] != "\x93NUMPY":
    raise newException(ValueError, "invalid NPY file: " & path)
  let major = ord(b[6])
  var hlen, off: int
  if major == 1:
    hlen = ord(b[8]) or (ord(b[9]) shl 8)
    off = 10 + hlen
  elif major in [2,3]:
    if b.len < 12: raise newException(ValueError, "truncated NPY header")
    hlen = ord(b[8]) or (ord(b[9]) shl 8) or (ord(b[10]) shl 16) or (ord(b[11]) shl 24)
    off = 12 + hlen
  else:
    raise newException(ValueError, "unsupported NPY version")
  let hs = b[(off-hlen)..<off]
  proc extractAfter(key: string): string =
    let p = hs.find(key)
    if p < 0: return ""
    let q1 = hs.find('\'', p + key.len)
    let q2 = hs.find('\'', q1 + 1)
    if q1 < 0 or q2 < 0: return ""
    hs[q1+1..<q2]
  let descr = extractAfter("'descr':")
  let sp = hs.find("'shape':")
  var n = 0
  if sp >= 0:
    let lp = hs.find('(', sp)
    let comma = hs.find(',', lp)
    if lp >= 0 and comma > lp:
      try: n = parseInt(hs[lp+1..<comma].strip())
      except ValueError: discard
  result = NpyData(descr: descr, count: n, offset: off, bytes: b)

proc u32le(s: string; p: int): uint32 =
  uint32(ord(s[p])) or (uint32(ord(s[p+1])) shl 8) or (uint32(ord(s[p+2])) shl 16) or (uint32(ord(s[p+3])) shl 24)
proc u64le(s: string; p: int): uint64 =
  uint64(u32le(s,p)) or (uint64(u32le(s,p+4)) shl 32)
proc lastFloat(path: string): float =
  let a = parseNpy(path)
  if a.count == 0: return NaN
  case a.descr
  of "<f8", "|f8": cast[float64](u64le(a.bytes, a.offset + (a.count-1)*8))
  of "<f4", "|f4": float(cast[float32](u32le(a.bytes, a.offset + (a.count-1)*4)))
  of "<i8", "|i8": float(cast[int64](u64le(a.bytes, a.offset + (a.count-1)*8)))
  of "<u8", "|u8": float(u64le(a.bytes, a.offset + (a.count-1)*8))
  of "<i4", "|i4": float(cast[int32](u32le(a.bytes, a.offset + (a.count-1)*4)))
  of "<u4", "|u4": float(u32le(a.bytes, a.offset + (a.count-1)*4))
  else: NaN
proc maxFloat(path: string): float =
  let a = parseNpy(path)
  result = NaN
  for i in 0..<a.count:
    var v: float
    case a.descr
    of "<f8", "|f8": v = cast[float64](u64le(a.bytes, a.offset+i*8))
    of "<f4", "|f4": v = float(cast[float32](u32le(a.bytes, a.offset+i*4)))
    of "<i8", "|i8": v = float(cast[int64](u64le(a.bytes, a.offset+i*8)))
    of "<u8", "|u8": v = float(u64le(a.bytes, a.offset+i*8))
    of "<i4", "|i4": v = float(cast[int32](u32le(a.bytes, a.offset+i*4)))
    of "<u4", "|u4": v = float(u32le(a.bytes, a.offset+i*4))
    else: return NaN
    if result.isNaN or v > result: result = v
proc npyCount(path: string): int = parseNpy(path).count
proc lastMaybe(root, rel: string): float =
  let p = root / rel
  if fileExists(p): lastFloat(p) else: NaN
proc countMaybe(root, rel: string): int =
  let p = root / rel
  if fileExists(p): npyCount(p) else: 0
proc dirSize(path: string): int64 =
  for p in walkDirRec(path):
    if fileExists(p): result += getFileSize(p)
proc fmt(v: float): string =
  if v.isNaN: "-" else: formatFloat(v, ffDefault, 6)
proc jnum(v: float): JsonNode =
  if v.isNaN: newJNull() else: %v

proc validationCounts(root: string): tuple[status:string,warnings:int,errors:int] =
  let p = root / "diagnostics" / "validation.jsonl"
  if not fileExists(p): return ("-",0,0)
  for line in lines(p):
    if line.strip.len == 0: continue
    try:
      let j = parseJson(line)
      let st = j{"status"}.getStr("").toLowerAscii
      let sev = j{"severity"}.getStr("").toLowerAscii
      if st in ["fail","failed","error"]:
        if sev == "warning": inc result.warnings else: inc result.errors
      elif st in ["warning","warn"]: inc result.warnings
    except CatchableError: discard
  result.status = if result.errors > 0: "FAIL" else: "PASS"

proc main() =
  let args = commandLineParams()
  if args.len == 1 and args[0] in ["-v", "--version"]:
    echo AppName & " " & AppVersion
    quit(0)
  if args.len == 0 or args[0] in ["-h","--help"]:
    echo AppName & " " & AppVersion
    echo "usage: slimmc-summary RUN [--format txt|json] [-o FILE]"
    quit(if args.len == 0: 1 else: 0)
  var root = args[0]
  var format = "txt"
  var output = ""
  var i = 1
  while i < args.len:
    case args[i]
    of "--format":
      inc i
      if i >= args.len: fail("missing --format value")
      format = args[i]
    of "-o", "--output":
      inc i
      if i >= args.len: fail("missing output path")
      output = args[i]
    else: fail("unknown option: " & args[i])
    inc i
  if not dirExists(root): fail("run directory not found: " & root)
  let mp = root / "run_metadata.json"
  if not fileExists(mp): fail("run_metadata.json not found")
  let md = parseJson(readFile(mp))
  let status = md{"run_status"}.getStr("unknown")
  let engine = md{"engine"}.getStr(md{"engine_name"}.getStr("-"))
  let engineVersion = md{"engine_version"}.getStr("-")
  let cliVersion = md{"cli_version"}.getStr("-")
  let seed =
    if md.hasKey("seed"): md["seed"].getStr($md["seed"])
    elif md.hasKey("rng_seed"): md["rng_seed"].getStr($md["rng_seed"])
    else: "-"
  let vc = validationCounts(root)
  let finalTime = lastMaybe(root, "snapshots/time.npy")
  let finalEvent = lastMaybe(root, "snapshots/kmc_event.npy")
  let snapshots = countMaybe(root, "snapshots/snapshot_id.npy")
  let mn = lastMaybe(root, "moments/mn.npy")
  let mw = lastMaybe(root, "moments/mw.npy")
  let mz = lastMaybe(root, "moments/mz.npy")
  let disp = lastMaybe(root, "moments/dispersity.npy")
  let dpn = lastMaybe(root, "moments/dpn.npy")
  let dpw = lastMaybe(root, "moments/dpw.npy")
  let peakMemory = if fileExists(root / "memory" / "total_est_B.npy"): maxFloat(root / "memory" / "total_est_B.npy") else: NaN
  let sizeB = dirSize(root)
  var outputText: string
  if format == "json":
    var j = newJObject()
    j["path"] = %extractFilename(root)
    j["engine"] = %engine; j["engine_version"] = %engineVersion; j["cli_version"] = %cliVersion
    j["status"] = %status; j["validation"] = %vc.status; j["validation_warnings"] = %vc.warnings; j["validation_errors"] = %vc.errors
    j["seed"] = %seed; j["final_time"] = jnum(finalTime); j["final_event"] = jnum(finalEvent); j["snapshots"] = %snapshots
    j["dpn"] = jnum(dpn); j["dpw"] = jnum(dpw); j["mn"] = jnum(mn); j["mw"] = jnum(mw); j["mz"] = jnum(mz); j["dispersity"] = jnum(disp)
    j["peak_memory_B"] = jnum(peakMemory); j["results_size_B"] = %sizeB
    outputText = pretty(j)
  elif format == "txt":
    outputText = "Slimmc run summary\n\n" &
      "Run:                 " & extractFilename(root) & "/\n" &
      "Engine:              " & engine & "\n" &
      "Engine version:      " & engineVersion & "\n" &
      "CLI version:         " & cliVersion & "\n" &
      "Status:              " & status & "\n" &
      "Validation:          " & vc.status & "\n" &
      "Warnings:            " & $vc.warnings & "\n" &
      "Errors:              " & $vc.errors & "\n" &
      "Seed:                " & seed & "\n\n" &
      "Final time:          " & fmt(finalTime) & "\n" &
      "Final KMC event:     " & fmt(finalEvent) & "\n" &
      "Snapshots:           " & $snapshots & "\n\n" &
      "Final moments:\n" &
      "  DPN                " & fmt(dpn) & "\n" &
      "  DPW                " & fmt(dpw) & "\n" &
      "  Mn                 " & fmt(mn) & "\n" &
      "  Mw                 " & fmt(mw) & "\n" &
      "  Mz                 " & fmt(mz) & "\n" &
      "  Dispersity         " & fmt(disp) & "\n\n" &
      "Peak memory:         " & fmt(peakMemory) & " B\n" &
      "Results size:        " & $sizeB & " B"
  else: fail("unsupported format: " & format)
  if output.len > 0: writeFile(output, outputText & "\n") else: echo outputText

when isMainModule: main()
