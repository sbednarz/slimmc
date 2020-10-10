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

import random/mersenne
import tables

# Avogadro const
const N_A = 6.023e23

# the simulation control volume (L)
var V_MC = 1.0e-17

# V_MC ~ N / (c_radicals * NA)
# c_radicals = 1e-7
# 10/(1e-7*6.023e23) = 1.66e-16 L

var
  cI0 = 0.0  #mol/L
  cRx0 = 0.0 #mol/L
  cM0 = 0.0  #mol/L
  cPx0 = 0.0 #mol/L
  cD0 = 0.0  #mol/L

  MwM = 0.0 #g/mol - molecular weight of the monomer
  
  kd = 0.0
  f = 0.0
  ki = 0.0
  kp = 0.0
  ktc = 0.0
  ktd = 0.0

var
  nI, nRx, nM, nPx, nD: int
  n_mc: int
  kd_MC, ki_MC, kp_MC, ktc_MC, ktd_MC: float
  rng: MersenneTwister

# memory usage: uint8 vs. int
var
  macroP: seq[seq[int]]
  macroD: seq[seq[int]]
  chain: seq[int]

macroP = newSeq[seq[int]]()
macroD = newSeq[seq[int]]()


var
  cldP: seq[int]
  cldD: seq[int]
  nChainsP: int
  nChainsD: int
  PnP: float
  PnD: float
  PwP: float
  PwD: float

# homopolymerization
const monomerunit = 1

var flags = initTable[string, int]()
flags["listparameters"]=0
flags["listbreakpoints"]=0
flags["listinitialstate"]=0
flags["outputeverystep"]=0
flags["outputpercent"]=0

