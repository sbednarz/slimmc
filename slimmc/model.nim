# slimmc - a simply and non-general use
# Monte Carlo simulation program of radical polymerization
#
# Copyright (C) 2020-2022 Szczepan Bednarz <sbednarz@pk.edu.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import strutils, tables, sets, algorithm
import strformat

include strings


var data: string

let valid_keyword_parameters = toHashSet(["kd", "f", "ki", "kp", "ktc", "ktd",
    "cI0", "cM0", "cRx0", "cPx0", "cD0", "V_MC", "MwM", "seed"])

let valid_dc_species = toHashSet(["I", "M", "Rx", "Px", "D"])


var parameters_list = initOrderedTable[string, string]()
var breakpoints_list = initOrderedTable[float, seq[string]]()

type
  Breakpoint = tuple
    time: float
    commands: seq[string]

var breakpoints: seq[Breakpoint]


const cmdconc = "1"
const cmdpoly = "2"
const cmddc = "3"
const cmdprint = "4"


#
# parser helpers
#

proc isFloat(s: string): bool =
  ## true if the whole (trimmed) string is a valid float / int literal
  try:
    discard parseFloat(s)
    result = true
  except ValueError:
    result = false

proc lineError(filename: string, lineNo: int, line: string) =
  echo "slimmc: error in model file, check ", filename, " line ", lineNo,
       ": ", line
  quit(-1)

proc cmdError(filename: string, lineNo: int, cmd: string) =
  echo "slimmc: error in model file '", filename, "' line ", lineNo,
       ": syntax error in command '", cmd, "'"
  quit(-1)

proc addBreakpoint(time: float, commands: seq[string]) =
  if breakpoints_list.hasKey(time):
    breakpoints_list[time].add(commands)
  else:
    breakpoints_list[time] = commands


proc parseCommands(cmdStrs: seq[string], filename: string,
                   lineNo: int): seq[string] =
  ## turns a comma-split command list into the flat, encoded command seq.
  ## encoding: conc -> "1"; poly -> "2";
  ##           dc -> "3", species, value;  print -> "4", text
  for raw in cmdStrs:
    let a = raw.strip()

    if a == "conc":
      result.add(cmdconc)

    elif a == "poly":
      result.add(cmdpoly)

    elif a.startsWith("dc"):
      # dc SPECIES VALUE
      let toks = a.splitWhitespace()
      if toks.len == 3 and toks[0] == "dc" and
         toks[1] in valid_dc_species and isFloat(toks[2]):
        result.add(cmddc)
        result.add(toks[1])
        result.add(toks[2])
      else:
        cmdError(filename, lineNo, raw)

    elif a.startsWith("print"):
      # print  |  print progress  |  print 'some text'
      let rest = a[5 .. ^1].strip()
      if rest.len == 0 or rest == "progress":
        result.add(cmdprint)
        result.add("progress")
      elif rest.len >= 2 and rest[0] == '\'' and rest[^1] == '\'':
        result.add(cmdprint)
        result.add(rest[1 .. ^2])          # text between the quotes
      else:
        cmdError(filename, lineNo, raw)

    else:
      cmdError(filename, lineNo, raw)


#
# model file parser
#
proc loadModel(filename: string) =

  data = readFile(filename)
  var lineNo = 0

  for line in splitLines(data):
    inc lineNo
    let s = line.strip()

    # comment ('#' in column 0) or blank line
    if line.len > 0 and line[0] == '#':
      continue
    if s.len == 0:
      continue

    # breakpoint(s): any line containing ':'  (checked before '=' so that a
    # ':' is unambiguous - a '=' may legitimately appear inside a print string)
    if ':' in line:
      # loop:  t0 : dt : N : commands   (maxsplit=3 keeps ':' inside commands)
      let lp = line.split(':', maxsplit = 3)
      if lp.len == 4 and isFloat(lp[0].strip) and
         isFloat(lp[1].strip) and isFloat(lp[2].strip):
        let t0 = parseFloat(lp[0].strip)
        let dt = parseFloat(lp[1].strip)
        let n  = int(parseFloat(lp[2].strip))
        let commands = parseCommands(split_cmds(lp[3]), filename, lineNo)
        for k in 0 .. n:
          addBreakpoint(t0 + dt * k.float, commands)
        continue

      # single:  time : commands
      let sp = line.split(':', maxsplit = 1)
      if not isFloat(sp[0].strip):
        lineError(filename, lineNo, line)
      let time = parseFloat(sp[0].strip)
      let commands = parseCommands(split_cmds(sp[1]), filename, lineNo)
      addBreakpoint(time, commands)
      continue

    # parameter:  name = value
    if '=' in line:
      let eq = line.find('=')
      let key = line[0 ..< eq].strip()
      let val = line[eq + 1 .. ^1].strip()
      if key notin valid_keyword_parameters:
        echo "slimmc: error in model file, check ", filename, " line ", lineNo,
             ", parameter ", key
        quit(-1)
      if not isFloat(val):
        lineError(filename, lineNo, line)
      parameters_list[key] = val
      continue

    # pre-simulation command:  list parameters | breakpoints | initialstate
    let toks = s.splitWhitespace()
    if toks.len == 2 and toks[0] == "list":
      case toks[1]
      of "parameters":   flags["listparameters"] = 1
      of "breakpoints":  flags["listbreakpoints"] = 1
      of "initialstate": flags["listinitialstate"] = 1
      else: lineError(filename, lineNo, line)
      continue

    lineError(filename, lineNo, line)


proc setVariable(variable: string, value: float) =
  case variable:
    of "kd":
      kd = value
    of "f":
      f = value
    of "ki":
      ki = value
    of "kp":
      kp = value
    of "ktc":
      ktc = value
    of "ktd":
      ktd = value
    of "cI0":
      cI0 = value
    of "cM0":
      cM0 = value
    of "cRx0":
      cRx0 = value
    of "cPx0":
      cPx0 = value
    of "cD0":
      cD0 = value
    of "V_MC":
      V_MC = value
    of "MwM":
      MwM = value
    of "seed":
      seed = (uint32)value
    else:
      echo "slimc: unrecognized parameter ", variable
      quit(-1)

proc initVariables() =
  for variable, value in parameters_list.pairs:
    setVariable(variable, parseFloat(value))

proc printVariables() =
  echo "simulation parameters:"
  echo "kd=", kd, " 1/s"
  echo "f=", f
  echo "ki=", ki, " L/(mol*s)"
  echo "kp=", kp, " L/(mol*s)"
  echo "ktc=", ktc, " L/(mol*s)"
  echo "ktd=", ktd, " L/(mol*s)"
  echo "cI0=", cI0, " mol/L"
  echo "cM0=", cM0, " mol/L"
  echo "cRx0=", cRx0, " mol/L"
  echo "cPx0=", cPx0, " mol/L"
  echo "cD0=", cD0, " mol/L"
  echo "MwM=", MwM, " g/mol"
  echo "V_MC=", V_MC, " L"
  if seed == 0:
    echo "seed=random"
  else:
    echo "seed=", seed


proc initBreakpoints() =
  if breakpoints_list.len == 0:
    echo "error: breakpoint(s) must be defined, check the model file"
    quit(-1)
  breakpoints_list.sort do (a, b: (float, seq[string])) -> int: cmp(a[0], b[0])
  for t, a in breakpoints_list.pairs:
    breakpoints.add((time: t, commands: a))


proc printBreakpoints() =
  echo "breakpoints list:"
  for i, b in breakpoints:
    var t = b.time
    var a = b.commands
    var j = 0
    var str = &"{i+1} t={t:.12e}s => "
    while j < a.len:
      if a[j] == cmdconc:
        str = str & "conc"
      elif a[j] == cmdpoly:
        str = str & "poly"
      elif a[j] == cmddc:
        str = str & &"dc {a[j+1]} {a[j+2]}"
        j = j+2
      elif a[j] == cmdprint:
        if a[j+1] == "progress":
          str = str & &"print progress"
        else:
          str = str & &"print \'{a[j+1]}\'"
        j = j+1
      if j < a.len-1:
        str = str &  ", "
      inc(j)
    echo str

