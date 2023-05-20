#!/usr/bin/env python3

import sys,math,lhepgen

data = []
m_targ,m_meson = 0.939565420,0.1349767

Q2,xb = 1.75, 0.36
with open(sys.argv[1]) as ff:
  for line in ff:
    lls = line.split()
    E0,mtp,s0,ds0 = [float(ll) for ll in lls[:4]]
    tmin = abs(lhepgen.getTmin(Q2,xb))
    mt = mtp + tmin

    out = ["{:.2f}".format(vv) for vv in [Q2,xb,Q2,xb]]
    out += ["{:.4f}".format(vv*1000) for vv in [mt,s0,ds0]]
    out += ["0\t0 0 0\t0 0 0", "{:.4f}".format(E0)]
    print("\t".join(out))

