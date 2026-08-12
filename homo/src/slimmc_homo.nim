# slimmc_homo.nim
# Command-line front end for the slimmc homo engine.
#
# Build:
#   nim c --path:src -d:release --opt:speed -o:slimmc src/slimmc_homo.nim
#
# Run:
#   ./slimmc model.model
#   ./slimmc model.model --check
#   ./slimmc model.model --trace
#   ./slimmc model.model --trace-channels
#   ./slimmc model.model 
#   ./slimmc model.model --debug
#
# runMain(args) is also the entry point a monolithic slimmc dispatcher
# (see cli/slimmc_monolith.nim) calls directly in-process, instead of
# execing this file as a separate binary -- it takes an explicit
# seq[string] rather than reading commandLineParams() itself so it has
# no dependency on *this* process's real argv.

import strutils
import os
import slimmc_types
import slimmc_parser
import slimmc_io
import slimmc_kmc

proc printUsage() =
  stderr.writeLine(AppName & " " & AppVersion)
  stderr.writeLine("")
  stderr.writeLine("usage:")
  stderr.writeLine("  slimmc model.model [--check] [--debug] [--trace-channels N] [--output-root PATH]")
  stderr.writeLine("  slimmc --help")
  stderr.writeLine("  slimmc --version")

proc runMain*(args: seq[string]) =
  if args.len == 0:
    printUsage()
    raise newException(ValueError, "missing model file")

  if args.len == 1 and (args[0] == "--help" or args[0] == "-h"):
    printUsage()
    return

  if args.len == 1 and (args[0] == "--version" or args[0] == "-v"):
    echo AppName & " " & AppVersion
    return

  let modelFile = args[0]
  var checkOnly = false
  var opts = RunOptions()
  var outputRoot = ""

  if args.len >= 2:
    var i = 1
    while i < args.len:
      case args[i]
      of "--check": checkOnly = true
      of "--trace-channels":
        # Item 64: --trace-channels requires a positive integer N
        # argument (max rows written to disk); missing/non-positive N
        # is a CLI error, never a silent default.
        if i + 1 >= args.len:
          printUsage()
          raise newException(ValueError, "--trace-channels requires a positive integer argument N")
        let nStr = args[i + 1]
        var n: int64
        try:
          n = parseBiggestInt(nStr)
        except ValueError:
          printUsage()
          raise newException(ValueError, "--trace-channels N: N must be a positive integer, got: " & nStr)
        if n <= 0:
          printUsage()
          raise newException(ValueError, "--trace-channels N: N must be a positive integer, got: " & nStr)
        opts.traceChannelsLimit = n
        inc i
      of "--debug": opts.debug = true
      of "--output-root":
        if i + 1 >= args.len:
          printUsage()
          raise newException(ValueError, "--output-root requires PATH")
        outputRoot = absolutePath(args[i + 1])
        inc i
      else:
        printUsage()
        raise newException(ValueError, "unknown option: " & args[i])
      inc i

  var m = parseModel(modelFile)
  if outputRoot.len > 0:
    let modelDir = parentDir(absolutePath(modelFile))
    let relativeOutput = relativePath(m.outputDir, modelDir)
    if relativeOutput == ".." or relativeOutput.startsWith(".." & DirSep):
      raise newException(ValueError, "internal output_dir resolution escaped the model directory")
    m.outputDir = outputRoot / relativeOutput
  if checkOnly:
    printCheck(m)
  else:
    runSimulation(m, opts)

when isMainModule:
  try:
    runMain(commandLineParams())
  except CatchableError as e:
    stderr.writeLine("ERROR: " & e.msg)
    quit(1)
