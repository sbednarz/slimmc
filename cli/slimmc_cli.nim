## slimmc -- single dispatcher binary between the homo (homopolymer)
## and copo (copolymer) engines.
##
## This links both engines directly into one executable and calls the
## right one's runMain() as a plain in-process proc call, based on the
## model file's declared monomer count (see dispatch_logic.nim). There
## is exactly one file to build, ship, and put on PATH -- no sibling
## binary is ever exec'd, and no "binary not found, build the family
## first" failure mode is possible, because there is nothing left to
## find at runtime.
##
## Both homo and copo declare identically-named top-level symbols
## (AppName, AppVersion, RunOptions, parseModel, printCheck,
## runSimulation, ...) -- a real but fully and cleanly resolvable
## collision: `import ... as` aliases plus fully-qualified
## `alias.symbol` calls compile without ambiguity. No changes were
## needed inside homo's or copo's own modules for this. The actual
## per-engine change needed was making each engine's CLI logic callable
## as `runMain(args: seq[string])` instead of only reading process argv
## directly -- confined entirely to each engine's single entry-point
## file (slimmc_homo.nim, slimmc_copo.nim); none of the ~8500 lines of
## actual simulation/parser/IO logic in either engine needed touching.
##
## homo/ and copo/ still each build their own standalone CLI binary
## (`make build` in either directory) for solo engine development --
## but nothing in this package's tests, checks, or validation workflows
## invokes those directly any more. Every script and Makefile target
## that needs to run a model goes through this one binary, which
## dispatches internally.
##
## Making every script go through this one dispatcher surfaced a real
## dispatch-rule gap: 20 of homo's own 44 engine validation cases are
## pure-kinetics models with zero `monomer` declarations (no polymer
## chain growth at all) -- these previously worked only because
## run_validation.py called homo's own binary directly, never through
## any dispatcher. The rule (dispatch_logic.nim) now treats 0 monomers
## as a legitimate homo case rather than a dispatcher-level error.

import std/[os, strutils]
import "../homo/src/slimmc_homo" as homoEngine
import "../copo/src/slimmc_copo" as copoEngine
import "../homo/src/slimmc_types" as homoTypes
import "../copo/src/copo_types" as copoTypes
import dispatch_logic
import ../common/run_id
import ../common/build_provenance

const AppName = "slimmc"
const
  CliVersion {.strdefine.} = "5.0.2"
  PyslimmcVersion {.strdefine.} = "5.0.1"
  PyslimmcOptVersion {.strdefine.} = "1.0.0"

proc fail(msg: string) =
  stderr.writeLine("ERROR: " & msg)
  quit(1)

proc printVersions() =
  echo "Slimmc " & CliVersion
  echo "pyslimmc " & PyslimmcVersion
  echo "pyslimmc-opt " & PyslimmcOptVersion

proc buildMode(): string =
  when defined(release): "release"
  elif defined(debug): "debug"
  else: "default"

proc optimizationMode(): string =
  when compileOption("opt", "speed"): "speed"
  elif compileOption("opt", "size"): "size"
  else: "none"

proc compilerBackend(): string =
  when defined(c): "C"
  elif defined(cpp): "C++"
  elif defined(objc): "Objective-C"
  else: "unknown"

proc printBuildInfo() =
  echo ""
  echo "Build"
  echo "  mode:          " & buildMode()
  echo "  optimization:  " & optimizationMode()
  echo "  Nim:           " & NimVersion
  echo "  backend:       " & compilerBackend()
  echo "  target:        " & hostOS & "/" & hostCPU
  echo "  compiled:      " & CompileDate & " " & CompileTime
  if BuildTimestampUtc.len > 0: echo "  build timestamp: " & BuildTimestampUtc
  if BuildGitCommit.len > 0: echo "  git commit:    " & BuildGitCommit
  if BuildGitTag.len > 0: echo "  git tag:       " & BuildGitTag
  if BuildGitDirty in ["true", "false"]: echo "  git dirty:     " & BuildGitDirty

proc verifyUnifiedSlimmcVersion() =
  if homoTypes.AppVersion != CliVersion or copoTypes.AppVersion != CliVersion:
    fail("inconsistent build: CLI, homo and copo must use one Slimmc version")

proc printShortUsage() =
  echo ""
  echo "Usage: " & AppName & " [options] model.model"
  echo "Run `" & AppName & " --help` for details."

proc printHelp() =
  printVersions()
  echo ""
  echo "Usage:"
  echo "  " & AppName & " [options] model.model"
  echo "  " & AppName & " -h | --help"
  echo "  " & AppName & " --version"
  echo ""
  echo "Options:"
  echo "  --check               Validate the model without running it."
  echo "  --debug               Enable engine debug diagnostics."
  echo "  --trace-channels N    Save at most N channel-trace rows."
  echo "  --output-root PATH    Put model result directories under PATH."
  echo ""
  echo "Examples:"
  echo "  " & AppName & " model.model"
  echo "  " & AppName & " --check model.model"
  echo "  " & AppName & " --output-root results --debug model.model"
  echo ""
  echo "Options must precede the model file; the model file is always last."

proc parseRunArgs(args: seq[string]): tuple[modelPath: string, forwardedArgs: seq[string]] =
  var i = 0
  while i < args.len:
    let arg = args[i]
    case arg
    of "--check", "--debug":
      result.forwardedArgs.add(arg)
    of "--trace-channels", "--output-root":
      if i + 1 >= args.len:
        fail(arg & " requires " & (if arg == "--trace-channels": "N" else: "PATH"))
      result.forwardedArgs.add(arg)
      result.forwardedArgs.add(args[i + 1])
      inc i
    else:
      if arg.startsWith("-"):
        fail("unknown option: " & arg)
      if result.modelPath.len > 0:
        fail("only one model file may be specified")
      result.modelPath = arg
      if i != args.high:
        fail("model file must be last; place all options before " & arg)
    inc i

  if result.modelPath.len == 0:
    fail("missing model file")

proc main() =
  verifyUnifiedSlimmcVersion()
  let args = commandLineParams()
  if args.len == 0:
    printVersions()
    printShortUsage()
    quit(0)
  if args.len == 1 and args[0] in ["-v", "-V", "--version"]:
    printVersions()
    printBuildInfo()
    quit(0)
  if args.len == 1 and args[0] in ["-h", "--help"]:
    printHelp()
    quit(0)

  let parsed = parseRunArgs(args)
  let modelPath = parsed.modelPath
  let forwardedArgs = parsed.forwardedArgs

  if not fileExists(modelPath):
    fail("cannot open: " & modelPath)

  try:
    discard validateModelRunId(modelPath)
  except ValueError as e:
    fail(e.msg)

  let n = dispatch_logic.countMonomers(modelPath)
  if n < 0:
    fail("cannot open: " & modelPath)
  if n > 3:
    fail(modelPath & ": " & $n & " `monomer` declarations found; homo needs 0-1, copo needs 2-3")

  let engineDir = dispatch_logic.engineForCount(n)
  let fullArgs = @[modelPath] & forwardedArgs

  try:
    case engineDir
    of "homo": homoEngine.runMain(fullArgs)
    of "copo": copoEngine.runMain(fullArgs)
    else: fail("internal dispatcher error: unhandled engine " & engineDir)
  except CatchableError as e:
    stderr.writeLine("ERROR: " & e.msg)
    quit(1)

when isMainModule:
  main()
