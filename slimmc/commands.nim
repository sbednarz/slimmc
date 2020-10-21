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

import math
import streams

# TODO
# stats commands (macromol memory, mcsteps ...)

var 
  hlogmP: seq[float]
  dwP: seq[float]
  hlogmD: seq[float]
  dwD: seq[float]

#----- concentration of species

proc mc2conc(): (float, float, float, float, float) =
  var cI = ((float)nI)/(N_A*V_MC)
  var cRx = ((float)nRx)/(N_A*V_MC)
  var cM = ((float)nM)/(N_A*V_MC)
  var cPx = ((float)nPx)/(N_A*V_MC)
  var cD = ((float)nD)/(N_A*V_MC)
  return (cI, cRx, cM, cPx, cD)


proc conc(step: int, time: float) =
  
  var cI, cRx, cM, cPx, cD: float
  (cI, cRx, cM, cPx, cD) = mc2conc()

  var f = newFileStream("cI.txt", fmAppend)
  f.writeLine(&"{time:.12e}\t{cI:.12e}\t{nI}")
  f.close()
  f = newFileStream("cRx.txt", fmAppend)
  f.writeLine(&"{time:.12e}\t{cRx:.12e}\t{nRx}")
  f.close()
  f = newFileStream("cM.txt", fmAppend)
  f.writeLine(&"{time:.12e}\t{cM:.12e}\t{nM}")
  f.close()
  f = newFileStream("cPx.txt", fmAppend)
  f.writeLine(&"{time:.12e}\t{cPx:.12e}\t{nPx}")
  f.close()
  f = newFileStream("cD.txt", fmAppend)
  f.writeLine(&"{time:.12e}\t{cD:.12e}\t{nD}")
  f.close()



#----- macromolecules properties

proc estimate_macro_memory(): int =
  var size = 0
  for e in macroP:
    size = size + e.len
  for e in macroD:
    size = size + e.len
  return size


proc calcMacroProperties() =
  #TODO Pn Pw
  nChainsP = macroP.len
  for unit in macroP:
    var DP = unit.len
    var c = cldP.len
    if c <= DP:
      cldP.add(newSeq[int](DP-c+1))
    inc cldP[DP]
  nChainsD = macroD.len
  for unit in macroD:
    var DP = unit.len
    var c = cldD.len
    if c <= DP:
      cldD.add(newSeq[int](DP-c+1))
    inc cldD[DP]
 
  hlogmP = newSeq[float]()
  dwP = newSeq[float]()
  var i = 0
  for dp,n in cldP:
    hlogmP.add(log10(((float)dp)*MwM))
    dwP.add(((float)n)*((float)dp)*((float)dp)*MwM*MwM)
  
  hlogmD = newSeq[float]()
  dwD = newSeq[float]()
  i = 0
  for dp,n in cldD:
    hlogmD.add(log10(((float)dp)*MwM))
    dwD.add(((float)n)*((float)dp)*((float)dp)*MwM*MwM)
  

proc writeMacroProperties(time: float) =
  var timelabel = &"{time:.12e}"
  var f = newFileStream("P."&timelabel&".cld.txt", fmWrite)
  f.writeLine(&"#nChains: {nChainsP}")
  for dp, n in cldP:
    if dp > 0:
      f.writeLine(&"{dp}\t{n}")
  f.close()
  f = newFileStream("D."&timelabel&".cld.txt", fmWrite)
  f.writeLine(&"#nChains: {nChainsD}")
  for dp, n in cldD:
    if dp > 0:
      f.writeLine(&"{dp}\t{n}")
  f.close()

  var i = 0
  f = newFileStream("P."&timelabel&".HlogM.txt", fmWrite)
  while i < hlogmP.len:
    f.writeLine(&"{hlogmP[i]:.12e}\t{dwP[i]:.12e}")
    inc(i)
  f.close()
  i=0
  f = newFileStream("D."&timelabel&".HlogM.txt", fmWrite)
  while i < hlogmD.len:
    f.writeLine(&"{hlogmD[i]:.12e}\t{dwD[i]:.12e}")
    inc(i)
  f.close()



proc poly(step: int, time: float) =
  calcMacroProperties()
  writeMacroProperties(time)



#--- dc = change concentration of species

proc dc(step: int, time: float, species: string, dcval: float) =
    
  var cI, cRx, cM, cPx, cD: float
  var new_cI, new_cRx, new_cM, new_cPx, new_cD: float

  # instantaneous concentrations
  (cI, cRx, cM, cPx, cD) = mc2conc()

  # the change of the concentration (dcval could be positive or negative val)
  case species:
    of "I":
      new_cI = cI + dcval
      nI = int(new_cI*V_MC*N_A)
    of "Rx":
      new_cRx = cRx + dcval
      nRx = int(new_cRx*V_MC*N_A)
    of "M":
      new_cM = cM + dcval
      nM = int(new_cM*V_MC*N_A)
    of "Px":
      new_cPx = cPx + dcval
      nPx = int(new_cPx*V_MC*N_A)
    of "D":
      new_cD = cD + dcval
      nD = int(new_cD*V_MC*N_A)

  #echo "nI=", nI
  #echo "nRx=", nRx
  #echo "nM=", nM
  #echo "nPx=", nPx
  #echo "nD=", nD
  #echo "N=", nI+nRx+nM+nPx+nD



#----- run actions at a brakpoint

proc runActions(step: int) =
  var b = breakpoints[step]
  var t = b.time
  var a = b.actions
  var i = 0
  while i < a.len:
    if a[i] == cmdconc:
      conc(step, t)
    elif a[i] == cmdpoly:
      poly(step, t)
    elif a[i] == cmddc:
      var species = a[i+1]
      var dcval = parseFloat(a[i+2])
      dc(step, t, species, dcval)
      i = i+2
    inc(i)
