#!/usr/bin/env python3

import sys, os

if len(sys.argv)<3 or not sys.argv[1].endswith('.fitoutput') or not os.path.isfile(sys.argv[1]):
  print("supply valid configurational file")
  sys.exit(111)
  
os.system('./switch.sh '+sys.argv[1] + '&& make -j4')

import lhepgen


outs = []
for fname in sys.argv[2:]:
  print(f"calculating {fname}...")
  qxs = set()
  with open(fname) as ff:
    for line in ff:
      lls = line.split()
      q2,xb = float(lls[0]), float(lls[1])
      qxs.add((q2,xb))


  prepath="pi0P"
  if 'neutron' in fname:
    lhepgen.SetReactionPar(3)
    prepath="pi0N"
  elif 'eta' in fname:
    lhepgen.SetReactionPar(2)
    prepath="etaP"
  else:
    lhepgen.SetReactionPar(1)

  outs = []
  for qsq,xbj in qxs:
    w2 = lhepgen.getWsq(qsq,xbj)
    xi = lhepgen.compassxi(qsq,xbj)

    pname = "preparation/{}/preparation_{:.4f}_{:.4f}.dat".format(prepath,qsq,w2)
    if not os.path.isfile(pname):
      print(f"{pname} is missing...")
      continue
    lhepgen.loadPreparationFromFile("preparation/{}/preparation_{:.4f}_{:.4f}.dat".format(prepath,qsq,w2),qsq,xi)

    tt = lhepgen.getTmin(qsq, xbj)
    while tt>-2:
      myAmps = lhepgen.getAmplitude(qsq,xi,xbj,tt);

      sigT = lhepgen.getCX(myAmps, w2)
      sigL = lhepgen.getCXL(myAmps, w2)
      sigTT = lhepgen.getCXTT(myAmps, w2, 0)
      sigLT = lhepgen.getCXLT(myAmps, w2, 0)

      tt -= 0.05

      outs.append((fname,qsq,xbj,tt,sigT,sigL,sigTT,sigLT))


with open(sys.argv[1]+'.dat', 'w') as ff:
  for out in outs:
    fmt = " {:.5f}"*(len(out)-1)
    fmt = "{} "+fmt+"\n"
    ff.write(fmt.format(*out))
