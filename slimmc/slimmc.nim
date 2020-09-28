# slimmc - a simply and non-general use
# Monte Carlo simulation program of radical polymerization kinetics
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

# To run: nim c -r slimmc.nim

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
  echo """
Usage: slimmc [option] modelfile 
Run kinetics simulation of radical polymerizaton.

Options:
         -h           display this help and exit
         --version    output version information and exit
  """

proc printVersionMsg() =
  echo &"""
slimmc version {gtag} {ghash} (built {build}) written by Szczepan Bednarz
  """

proc printWelcomeMsg(modelfile: string) =
  echo "[slimmc]"
  echo "processing model file: ",modelfile

var modelfile = ""

#errors.txt

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
      else:
        printHelpMsg()
        quit()
    of cmdArgument:
      modelfile = p.key


printWelcomeMsg(modelfile)

loadModel(modelfile)

initVariables()
printVariables()

initBreakpoints()
printBreakpoints()

runSimulation()
