# slimmc_copo.nim
# Command-line front end for the slimmc copo engine.
#
# Build:
#   nim c -d:release --opt:speed -o:slimmc-copo src/slimmc_copo.nim
#
# Run:
#   ./slimmc-copo path/to/model.model
#   ./slimmc-copo model.model --check
#
# runMain(args) is also the entry point a monolithic slimmc dispatcher
# (see cli/slimmc_monolith.nim) calls directly in-process, instead of
# execing this file as a separate binary -- it takes an explicit
# seq[string] rather than reading commandLineParams() itself, and raises
# CatchableError on any usage/runtime error instead of calling quit()
# directly (matching slimmc_homo.nim's runMain() exactly), so callers --
# this file's own `when isMainModule:` below, or a merged dispatcher --
# decide independently how to turn that into a process exit code.

import os
import strutils
import copo_types
import copo_parser
import copo_kmc
import copo_io
import copo_storage

proc printUsage() =
  stderr.writeLine(AppName & " " & AppVersion)
  stderr.writeLine("")
  stderr.writeLine("usage:")
  stderr.writeLine("  slimmc-copo model.model [--check] [--debug] [--trace-channels N] [--output-root PATH]")
  stderr.writeLine("  slimmc-copo --help")
  stderr.writeLine("  slimmc-copo --version")

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
  var opts = RunOptions()
  var outputRoot = ""

  if args.len >= 2:
    var i = 1
    while i < args.len:
      case args[i]
      of "--check": opts.checkOnly = true
      of "--debug": opts.debug = true
      of "--output-root":
        if i + 1 >= args.len:
          printUsage()
          raise newException(ValueError, "--output-root requires PATH")
        outputRoot = absolutePath(args[i + 1])
        inc i
      of "--trace-channels":
        # Item 64: --trace-channels requires a positive integer N
        # argument (max rows written to disk); missing/non-positive N
        # is a CLI error, never a silent default. Matches homo exactly.
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
      else:
        printUsage()
        raise newException(ValueError, "unknown option: " & args[i])
      inc i

  var m = parseModel(modelFile)
  if outputRoot.len > 0:
    let modelDir = parentDir(absolutePath(modelFile))
    let relativeOutput = relativePath(m.output_dir, modelDir)
    if relativeOutput == ".." or relativeOutput.startsWith(".." & DirSep):
      raise newException(ValueError, "internal output_dir resolution escaped the model directory")
    m.output_dir = outputRoot / relativeOutput
  if opts.checkOnly:
    printCheck(m)
  else:
    try:
      runSimulation(m, opts)
    except CatchableError as e:
      markFailedStorageV1(m, e.msg, 1)
      raise

when isMainModule:
  try:
    runMain(commandLineParams())
  except CatchableError as e:
    stderr.writeLine("ERROR: " & e.msg)
    quit(1)
