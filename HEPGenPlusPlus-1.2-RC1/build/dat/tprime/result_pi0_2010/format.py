#!/usr/bin/env python3

import sys,math,lhepgen

data = []
m_targ,m_meson = 0.93827208816,0.1349767

xb = 0.36
for line in sys.stdin:
  lls = line.split()
  Q2 = float(lls[0])
  mtp = float(lls[1])

  st,slt,stt = [float(ll) for ll in lls[2::5]]
  dst,dslt,dstt = [float(ll) for ll in lls[3::5]]

  tmin = abs(lhepgen.getTmin(Q2,xb))
  mt = mtp + tmin

  out = ["{:.2f}".format(vv) for vv in [Q2,xb,Q2,xb]]
  out += ["{:.5f}".format(mt)]
  out += ["{:.4f}".format(vv*1000) for vv in [st,dst,0,slt,dslt,0,stt,dstt,0]]
  out += ["0"]
  print("\t".join(out))

