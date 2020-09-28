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
import times
import math
import strformat
import streams
import times

include commands

# Implemented scheme of radical polymerization:
#1  I => 2Rx f*kd
#2  M + Rx => Px ki
#3  Px + M => Px kp
#4  Px + Px => D ktc
#5  Px + Px => D + D ktd


proc initSimulation() =

  # pseudorandom number generator initialization
  var seed = uint32(toUnix(getTime()))
  rng = initMersenneTwister(seed)

  kd_MC = kd
  ki_MC = ki/(V_MC*N_A)
  kp_MC = kp/(V_MC*N_A)
  ktc_MC = ktc/(V_MC*N_A)
  ktd_MC = ktd/(V_MC*N_A)

  nI = int(cI0*V_MC*N_A)
  nRx = int(cRx0*V_MC*N_A)
  nM = int(cM0*V_MC*N_A)
  nPx = int(cPx0*V_MC*N_A)
  nD = int(cD0*V_MC*N_A)


proc printPostInitMsg() =
  echo "simulation initialized"
  echo "initial state:"
  echo "nI=", nI
  echo "nRx=", nRx
  echo "nM=", nM
  echo "nPx=", nPx
  echo "nD=", nD
  echo "N=", nI+nRx+nM+nPx+nD
  echo "simulation is running..."


proc gotoTime(t0: float, t_step: float) =

  var
    tt: float
    Rd_MC, Ri_MC, Rp_MC, Rtc_MC, Rtd_MC: float
    sumR, pd, pi, pp, ptc, ptd, r, tau: float

  tt = t0
  while (tt < t_step):

    Rd_MC = kd_MC*((float)nI) #*2?
    Ri_MC = ki_MC*((float)nRx)*((float)nM)
    Rp_MC = kp_MC*((float)nPx)*((float)nM)
    Rtc_MC = ktc_MC*((float)nPx)*((float)nPx)
    Rtd_MC = ktd_MC*((float)nPx)*((float)nPx)

    sumR = Rd_MC + Ri_MC + Rp_MC + Rtc_MC + Rtd_MC

    pd = Rd_MC/sumR
    pi = Ri_MC/sumR
    pp = Rp_MC/sumR
    ptc = Rtc_MC/sumR
    ptd = Rtd_MC/sumR

    r = rng.random()

    # decomposition
    if r <= pd:
      if nI > 0:
        nI = nI - 1
        r = rng.random()
        if r <= f:
          nRx = nRx + 2

    # initiation
    elif r <= pd+pi:
      if nRx > 0 and nM > 0:
        nRx = nRx - 1
        nM = nM - 1
        nPx = nPx + 1
        chain = newSeq[int](1)
        chain[0] = monomerunit
        macroP.add(chain)

    # propagation
    elif r <= pd+pi+pp:
      if nM > 0:
        nM = nM - 1
        var maxN = macroP.len
        var i = (int)(rng.random()*(float)maxN)
        chain = newSeq[int](1)
        chain[0] = monomerunit
        macroP[i].add(chain)

    # termination by comb
    elif r <= pd+pi+pp+ptc:
      if nPx > 1:
        nPx = nPx - 2
        nD = nD + 1

        var maxN = macroP.len
        var i = (int)(rng.random()*(float)maxN)
        var j = i
        while (j == i):
          j = (int)(rng.random()*(float)maxN)
        chain = macroP[i] & macroP[j]
        macroD.add(chain)
        macroP.delete(i)
        macroP.delete(j)


    # termination by disp
    elif r <= pd+pi+pp+ptc+ptd:
      if nPx > 1:
        nPx = nPx - 2
        nD = nD + 2

        var maxN = macroP.len
        var i = (int)(rng.random()*(float)maxN)
        var j = i
        #select 2 different macroradicals
        while (j == i):
          j = (int)(rng.random()*(float)maxN)
        macroD.add(macroP[i])
        macroD.add(macroP[j])
        macroP.delete(i)
        macroP.delete(j)


    r = rng.random()
    tau = -ln(r)/sumR
    tt = tt + tau
    n_mc = n_mc + 1

#minimal val of tau?
#sumR =>0 tau=>?
#simulation freezing?
#t_step vs. tau


proc runSimulation() =
  initSimulation()
  printPostInitMsg()

  var time0 = 0.0
  var time1 = 0.0
  var step = 0
  var nSteps = breakpoints.len

  var realtime0 = now()
  var realtime1 = now()
  var duration = realtime1 - realtime0

  if breakpoints[step].time == 0:
    runActions(step)
    echo "t=", time1, "s (", step+1,"/",nSteps,")"

  while true:
    inc (step)
    if step == nSteps:
      break

    time1 = breakpoints[step].time
    gotoTime(time0, time1)
    runActions(step)
    time0 = time1
    realtime1 = now()
    duration = realtime1 - realtime0
    echo "t=", time0, "s (", step+1,$"/",nSteps,") ", duration.inSeconds(),"s (~",duration.inHours(),"min)" 
