# slimmc - a simply and non-general use
# Monte Carlo simulation program of radical polymerization kinetics
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


import random/mersenne
import times
import math
import strformat
import streams
import parseopt
import os

include variables
include model
include kmc
include version

proc printHelpMsg() =
  echo &"""
Usage: {prg} [option] modelfile 
Run kinetics simulation of radical polymerizaton.
More info at: https://github.com/sbednarz/slimmc

Options:

         -t or --test  check model file syntax, do not run simulation
         -h            display this help and exit
         --version     output version information and exit"""

  if extra!="":
    echo extra



proc printVersionMsg() =
  echo &"""
{prg} {gtag} (built {build})
git sha: {ghash}
compiled on {sys} by {nimv} and {gcc}"""

  if extra!="":
    echo extra



proc printWelcomeMsg(modelfile: string) =
  echo "[slimmc]"
  echo "processing model file: ",modelfile

var modelfile = ""


if paramCount() == 0:
  printHelpMsg()
  quit()

var p = initOptParser(cmdLine="")
while true:
  p.next()
  case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      if p.key == "h":
        printHelpMsg()
        quit()
      elif p.key == "version":
        printVersionMsg()
        quit()
      elif p.key == "t":
        flags["modelsyntaxtest"]=1
      elif p.key == "test":
        flags["modelsyntaxtest"]=1
      else:
        printHelpMsg()
        quit()
    of cmdArgument:
      modelfile = p.key


printWelcomeMsg(modelfile)

loadModel(modelfile)

initVariables()
if flags["listparameters"]==1:
  printVariables()

initBreakpoints()
if flags["listbreakpoints"]==1:
  printBreakpoints()

runSimulation()
