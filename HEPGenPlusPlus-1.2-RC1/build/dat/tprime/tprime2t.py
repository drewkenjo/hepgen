#!/usr/bin/env python3

import sys,math

m_targ,m_meson = 0.93827208816,0.1349767

with open(sys.argv[1]) as ff:
  for line in ff:
    lls = line.split()
    qnom,xnom,Q2,xb,mt,s0,ds0,dds0,slt,dslt,ddslt,stt,dstt,ddstt = [float(ll) for ll in lls]

    W2 = Q2*(1/xb-1)+m_targ*m_targ;
    m12 =-Q2;
    m22 = m_targ*m_targ;
    m32 = m_meson*m_meson;
    m42 = m_targ*m_targ;
    
    e1cm = (W2 + m12 - m22) / (2*math.sqrt(W2));
    e3cm = (W2 + m32 - m42) / (2*math.sqrt(W2));
    
    p1cm = math.pow(e1cm,2) - m12;
    p3cm = math.pow(e3cm,2) - m32;
    
    p1cm = math.sqrt(p1cm);
    p3cm = math.sqrt(p3cm);

    tmin = math.pow(m12-m32-m22+m42,2)/4.0/W2 - math.pow(p1cm-p3cm,2.0);

    mt = mt + abs(tmin)

    lls[4] = "{:.4f}".format(mt)

    print(tmin)

    print(" ".join(lls))


