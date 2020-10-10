# slimmc - a simply and non-general use
# Monte Carlo simulation program of radical polymerization
#
# Copyright (C) 2020 Szczepan Bednarz <sbednarz@pk.edu.p>
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

import nregex, strutils, tables, sets, algorithm
import strformat


var data: string
var m: RegexMatch

let valid_keyword_parameters = toHashSet(["kd", "f", "ki", "kp", "ktc", "ktd",
    "cI0", "cM0", "cRx0", "cPx0", "cD0", "V_MC", "MwM"])


var parameters_list = initOrderedTable[string, string]()
var breakpoints_list = initOrderedTable[string, seq[string]]()

type
  Breakpoint = tuple
    time: float
    actions: seq[string]

var breakpoints: seq[Breakpoint]


const cmdconc = "1"
const cmdpoly = "2"
const cmddc = "3"




proc loadModel(filename: string) =

  data = readFile(filename)
  var i = 0

  for line in splitLines(data):
    inc i
    
    # a comment
    if match(line, re"^(\#).*", m):
      continue

    # an empty line
    if match(line, re"^\s*$", m):
      continue

    # parameter = val
    if match(line, re"\s*([a-zA-Z_0-9]+)\s*\=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*", m):
      if valid_keyword_parameters.contains(line[m.group(0)[0]]):
        parameters_list[line[m.group(0)[0]]] = line[m.group(1)[0]]
        continue
      else:
        echo "slimmc: error in model file, check ", filename, " line ", i,
            " : parameter ", m.group(0)
        quit(-1)

    # command arg1
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
          echo "slimmc: error in model file, check ", filename, " line ", i
          quit(-1)
        continue

      elif cmd=="output":
        if arg=="everystep":
          flags["outputeverystep"]=1
        elif arg=="percent":
          flags["outputpercent"]=1
        else:
          echo "slimmc: error in model file, check ", filename, " line ", i,
              " : command ", m.group(0)
          quit(-1)
        continue
      else:
        echo "slimmc: error in model file, check ", filename, " line ", i,
            " : parameter ", m.group(0)
        quit(-1)



    # a series of breakpoints (loop) 
    if match(line, re"\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\:\s*(.*)", m):

      var actions: seq[string]
      var t0 = parseFloat(line[m.group(0)[0]])
      var dt = parseFloat(line[m.group(2)[0]])
      var N = ((int)parseFloat(line[m.group(4)[0]]))
      var actions_str = split(line[m.group(6)[0]], ",")
      for a in actions_str:
        if match(a, re"^\s*(conc)\s*$"):
          actions.add(cmdconc)
          continue
        if match(a, re"^\s*(poly)\s*$"):
          actions.add(cmdpoly)
          continue
        if match(a, re"^\s*(dc)\s+(I|M|Rx|Px|D)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*$", m):
          actions.add(cmddc)
          actions.add(a[m.group(1)[0]])
          actions.add(a[m.group(2)[0]])
          continue

        echo "slimmc: error in model file \'", filename, "\' line ", i,
            ": syntax error in command \'", a, "\'"
        quit(-1)
  
      var i = 0
      while i < N:
        var time = t0 + dt*((float)i)
        var time_str = &"{time:.12f}"
        if breakpoints_list.hasKey(time_str):
          breakpoints_list[time_str].add(actions)
        else:
          breakpoints_list[time_str]=actions
        inc i
      continue


    # a breakpoint
    if match(line, re"\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\:(.*)", m):
      var actions: seq[string]
      var time = parseFloat(line[m.group(0)[0]])
      var actions_str = split(line[m.group(2)[0]], ",")
      for a in actions_str:
        if match(a, re"^\s*(conc)\s*$"):
          actions.add(cmdconc)
          continue
        if match(a, re"^\s*(poly)\s*$"):
          actions.add(cmdpoly)
          continue
        if match(a, re"^\s*(dc)\s+(I|M|Rx|Px|D)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*$", m):
          actions.add(cmddc)
          actions.add(a[m.group(1)[0]])
          actions.add(a[m.group(2)[0]])
          continue

        echo "noplp"
        echo "slimmc: error in model file \'", filename, "\' line ", i,
            ": syntax error in command \'", a, "\'"
        quit(-1)

      var time_str = &"{time:.12f}"
      if breakpoints_list.hasKey(time_str):
        breakpoints_list[time_str].add(actions)
      else:
        breakpoints_list[time_str]=actions
      continue


    echo "slimmc: error in model file \'", filename, "\' line ", i
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


proc initBreakpoints() =
  if breakpoints_list.len == 0:
    echo "error: breakpoint(s) must be defined, check the model file"
    quit(-1)
  breakpoints_list.sort do (a, b: (string, seq[string])) -> int: cmp(a[0], b[0])
  for t, a in breakpoints_list.pairs:
    breakpoints.add((time: parseFloat(t), actions: a))


proc printBreakpoints() =
  echo "breakpoints list:"
  for i, b in breakpoints:
    var t = b.time
    var a = b.actions
    var j = 0
    var str = &"{i+1} t={t:.12f}s => "
    while j < a.len:
      if a[j] == cmdconc:
        str = str & "conc"
      elif a[j] == cmdpoly:
        str = str & "poly"
      elif a[j] == cmddc:
        str = str & &"dc {a[j+1]} {a[j+2]}"
        j = j+2
      if j < a.len-1:
        str = str &  ", "
      inc(j)
    echo str


