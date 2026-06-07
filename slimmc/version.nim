
import std/strutils

proc capture(cmd: string): string {.compileTime.} =
  let (output, code) = gorgeEx(cmd)
  if code == 0: output.strip() else: "unknown"

const
  prg*   {.strdefine.} = "slimmc"          # nadpisywalne: -d:prg=slimmc-turbo
  extra* {.strdefine.} = ""                # -d:extra="..."
  gtag*  = capture("git describe --abbrev=4 --dirty --always --tags")
  ghash* = capture("git rev-parse HEAD")
  build* = CompileDate & " " & CompileTime
  sys*   = hostOS & "/" & hostCPU
  nimv*  = "Nim " & NimVersion
  gcc*   = capture("gcc --version").splitLines()[0]

