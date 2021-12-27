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

import regex, strutils, tables, sets, algorithm
import strformat

include strings


var data: string
var m: RegexMatch

let valid_keyword_parameters = toHashSet(["kd", "f", "ki", "kp", "ktc", "ktd",
    "cI0", "cM0", "cRx0", "cPx0", "cD0", "V_MC", "MwM", "seed"])


var parameters_list = initOrderedTable[string, string]()
#var breakpoints_list = initOrderedTable[string, seq[string]]()
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
# model file parser
#
proc loadModel(filename: string) =

  data = readFile(filename)
  var i = 0

  for line in splitLines(data):
    inc i
    
    # skip a comment => "#" starts a comment line
    if match(line, re"^(\#).*", m):
      continue

    # skip an empty line
    if match(line, re"^\s*$", m):
      continue

    # parameter = val
    if match(line, re"\s*([a-zA-Z_0-9]+)\s*\=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*", m):
      if valid_keyword_parameters.contains(line[m.group(0)[0]]):
        parameters_list[line[m.group(0)[0]]] = line[m.group(1)[0]]
        continue
      else:
        echo "slimmc: error in model file, check ", filename, " line ", i,
            ", parameter ", line[m.group(0)[0]]
        quit(-1)

    # commands to be run before simulation starts => command arg1
    if match(line, re"\s*([a-zA-Z_0-9]+)\s+([a-zA-Z_0-9]+)", m):
      var cmd = line[m.group(0)[0]]
      var arg = line[m.group(1)[0]]

      if cmd=="list":
        if arg=="parameters":
          flags["listparameters"]=1
        elif arg=="breakpoints":
          flags["listbreakpoints"]=1
        elif arg=="initialstate":
          flags["listinitialstate"]=1
        else:
          echo "slimmc: error in model file, check ", filename, " line ", i, " command ",cmd, " arg ", arg
          quit(-1)
        continue



    # a series of breakpoints (loop) => start time (s) : time step (s) : number of steps : comma seperated commands
    if match(line, re"\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*(.*)", m):

      var commands: seq[string]
      var t0 = parseFloat(line[m.group(0)[0]])
      var dt = parseFloat(line[m.group(2)[0]])
      var N = ((int)parseFloat(line[m.group(4)[0]]))
      var commands_str = split_cmds(line[m.group(6)[0]])
      for a in commands_str:
        if match(a, re"^\s*(conc)\s*$"):
          commands.add(cmdconc)
          continue
        if match(a, re"^\s*(poly)\s*$"):
          commands.add(cmdpoly)
          continue
        if match(a, re"^\s*(dc)\s+(I|M|Rx|Px|D)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*$", m):
          commands.add(cmddc)
          commands.add(a[m.group(1)[0]])
          commands.add(a[m.group(2)[0]])
          continue
        if match(a, re"^\s*(print)\s+(progress)\s*$", m):
          commands.add(cmdprint)
          commands.add(a[m.group(1)[0]])
          continue
        if match(a, re"^\s*(print)\s*$", m):
          commands.add(cmdprint)
          commands.add("progress")
          continue
        if match(a, re"^\s*(print)\s+\'([a-zA-Z0-9_\!\?\-\+\*\,\.\%\#\;\:\=\(\)\[\]\<\>\/\\\ ]*)\'\s*$", m):
          commands.add(cmdprint)
          commands.add(a[m.group(1)[0]])
          continue

        echo "slimmc: error in model file \'", filename, "\' line ", i,
            ": syntax error in command \'", a, "\'"
        quit(-1)
  
      var i = 0
      while i < N+1:
        var time = t0 + dt*((float)i)
        if breakpoints_list.hasKey(time):
          breakpoints_list[time].add(commands)
        else:
          breakpoints_list[time]=commands
        inc i
      continue


    # a breakpoint => time : comma seperated commands
    if match(line, re"\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\:(.*)", m):
      var commands: seq[string]
      var time = parseFloat(line[m.group(0)[0]])
      var commands_str = split_cmds(line[m.group(2)[0]])
      for a in commands_str:
        if match(a, re"^\s*(conc)\s*$"):
          commands.add(cmdconc)
          continue
        if match(a, re"^\s*(poly)\s*$"):
          commands.add(cmdpoly)
          continue
        if match(a, re"^\s*(dc)\s+(I|M|Rx|Px|D)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*$", m):
          commands.add(cmddc)
          commands.add(a[m.group(1)[0]])
          commands.add(a[m.group(2)[0]])
          continue
        if match(a, re"^\s*(print)\s+(progress)\s*$", m):
          commands.add(cmdprint)
          commands.add(a[m.group(1)[0]])
          continue
        if match(a, re"^\s*(print)\s*$", m):
          commands.add(cmdprint)
          commands.add("progress")
          continue
        if match(a, re"^\s*(print)\s+\'([a-zA-Z0-9_\!\?\-\+\*\,\.\%\#\;\:\=\(\)\[\]\<\>\/\\\ ]*)\'\s*$", m):
          commands.add(cmdprint)
          commands.add(a[m.group(1)[0]])
          continue


        echo "noplp"
        echo "slimmc: error in model file \'", filename, "\' line ", i,
            ": syntax error in command \'", a, "\'"
        quit(-1)

      if breakpoints_list.hasKey(time):
        breakpoints_list[time].add(commands)
      else:
        breakpoints_list[time]=commands
      continue


    echo "slimmc: error in model file, check ", filename, " line ", i,": ", line
    quit(-1)



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


