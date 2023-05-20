#!/usr/bin/env python3

import sys,math,lhepgen

data = []
m_targ,m_meson = 0.939565420,0.1349767

Q2,xb = 1.75, 0.36
with open(sys.argv[1]) as ff:
  for line in ff:
    lls = line.split()
    mtp = float(lls[0])
    st,sl,stt,slt = [float(ll) for ll in lls[1::5]]
    dst,dsl,dstt,dslt = [float(ll) for ll in lls[2::5]]

    tmin = abs(lhepgen.getTmin(Q2,xb))

    mt = mtp + tmin

    out = ["{:.2f}".format(vv) for vv in [Q2,xb,Q2,xb]]
    out += ["{:.5f}".format(mt)]
    out += ["{:.4f}".format(vv*1000) for vv in [st,dst,0,slt,dslt,0,stt,dstt,0]]
    out += ["0"]
    print("\t".join(out))

