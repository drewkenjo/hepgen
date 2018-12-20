#!/usr/bin/python

import math, sys
import lhepgen

# Get the xsec for some specific kinematic point.
# While Qsq is mandatory, you can choose between nu/W/xbj and t/t'.
def get_xsec(Qsq, nu=None, W=None, xbj=None, t=None, tprim=None, phi=0.0, Ein=5.75, table=None):
    print(sys.argv)
    m_P = 0.93827203

    # Create an instance of 'hWeightInterface'. It is used to define the kinematics  used in the calculation.
    w = lhepgen.hWeightInterface()
    w.MuIn.setEnergy(Ein)
    w.beamE = Ein
    w.qsq = Qsq

    # Depending on the kinematic variables input to the function, calculate the other variables/
    if W is None and xbj is None:
        w.nu=nu
        w.xbj=w.qsq/(2.*m_P*w.nu)
        W=lhepgen.getWsq(Qsq,w.xbj) # actually returns sqrt(wsq) !!
        w.w2 = W**2
        if w.xbj<0 or w.xbj>1. or W<m_P or w.nu<0. or w.nu>Ein:
          return 0.0
    elif W is not None:
        w.w2 = W**2
        w.xbj = lhepgen.getXbj(Qsq,W)
        w.nu = w.qsq/(2.*m_P*w.xbj)
        if w.xbj<0 or w.xbj>1. or W<m_P or w.nu<0. or w.nu>Ein:
          return 0.0 
    else:
        w.xbj = xbj
        W = lhepgen.getWsq(Qsq,w.xbj) # actually returns sqrt(wsq) !!
        w.w2 = W**2
        w.nu = w.qsq/(2.*m_P*w.xbj)
        if w.xbj<0 or w.xbj>1. or W<m_P or w.nu<0. or w.nu>Ein:
          return 0.0 

    w.y = w.nu/Ein
    if w.y<0.0 or w.y>1.0:
       return 0.0    

    # Calculate xi and t0
    xi = lhepgen.compassxi(w.xbj)
    t0 = -4.*(m_P**2)*(xi**2)/(1.-(xi**2))
    if t>t0 :
         return 0.0
    # The interface for the amplitudes uses t
    if tprim is None:
        w.t = t
        w.tprim = t0-w.t
    else:
        w.tprim = tprim
        w.t = -w.tprim+t0

    # Calculate the subprocess amplitude of the process. Use pre-calculated subprocess tables if 'table' is set.
    if table is not None:
        myAmps = table.getAmpsForKine(w.qsq, w.xbj, w.t)
    else:
        lhepgen.prepareConvolution(w.qsq, xi)
        myAmps = lhepgen.getAmplitude(w.qsq, xi, w.xbj, w.t)

    # Calculate photon polarization factor epsilon
    epsilon = lhepgen.hreweightKine.getEpsilon(w)
    if epsilon >1.0 or epsilon<0.0:
       return 0.0
    # Calculate the value of each structure function.
    valsig = lhepgen.getCX(myAmps, W)
    valsigL = lhepgen.getCXL(myAmps, W)
    valsigTT, valsigLT = 0, 0
    # Only calculate sigma_LT and sigma_TT if phi is given.
    if phi is not None:
        valsigTT = lhepgen.getCXTT(myAmps, W, phi)
        valsigLT = lhepgen.getCXLT(myAmps, W, phi)
#    sig=valsig + valsigL*epsilon + epsilon*valsigTT + math.sqrt(2*epsilon*(1+epsilon))*valsigLT
    sig=valsig + valsigL*epsilon 

    print ("{0:3.2f} {1:5.3f} {2:5.3f} {3:5.3f} {4:5.3f} {5:5.3f} {6:6.3f} {7:6.3f} {8:8.3f} {9:8.3f} {10:8.3f} {11:8.3f} {12:8.3f} {13:8.3f}".format(Ein, W, w.xbj, w.qsq, -t, -t0, -w.tprim, xi, phi, valsig, valsigL, valsigTT, valsigLT, sig))  

    return sig

#================================================================================================#

if __name__ == "__main__":
    # hepgen++ provides a grid of pre-calculated subprocess amplitudes. Using them
    # speeds up the calculations alot, however at the expense of some accuracy.

    use_subpoc_table = False
#    use_subpoc_table = True

    # If requested, load subprocess tables from disk (side note: location of
    # the tables is '$HEPGEN/share/pi0_cache/output/').
    subpoc_table = None
    if use_subpoc_table:
        print('---------------------------------------------------------------------------------')
        print('===')
        print('===                      WITH SUBPROCESS TABLES')
        print('===')
        print('---------------------------------------------------------------------------------')
        subpoc_table = lhepgen.gkSubProcessTableCache()
        subpoc_table.loadCache()
    else:
        print('---------------------------------------------------------------------------------')
        print('===')
        print('===                      WITHOUT SUBPROCESS TABLES')
        print('===')
        print('---------------------------------------------------------------------------------')

# get the xsec  for a specific kinematic point and print it

#    Q2list=(1.0, 2.,2.2, 3., 4., 5.)
#    xBlist=(0.1,0.2,0.28,0.3,0.4,0.5)
#    tlist= (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.5,2.)


    print ("{0:5} {1:3} {2:6} {3:6} {4:6} {5:4} {6:6} {7:4} {8:5} {9:5} {10:7} {11:7} {12:7} {13:7}".format(    " Ein", " W ","   xB","  Q2", "  -t", "-t0", "-tprim","  xi ","     phi", "  sigT", "     sigL", "    sigTT  ", " sigLT  ", "sig"))  

#    for Q2 in Q2list:
#       for xB in xBlist:

#        pi0/eta out of proton
    tlist= (0.12,0.17,0.25,0.35,0.49,0.78,1.22)
    for tl in tlist:
       get_xsec(Qsq=2.21, xbj=0.275, t=-tl, Ein=5.75, table=subpoc_table) 




#-------  pi0 out of neutron 
#    tlist= (0.02,0.07,0.12,0.17) # this is t0-t
#    m_P = 0.93827203
#    xbj=0.360
#    xi = xbj/(2-xbj)
#    t0 = -4.*(m_P**2)*(xi**2)/(1.-(xi**2))
#--------
#    for tl in tlist:
#-#       get_xsec(Qsq=2.21, xbj=0.275, t=-tl, Ein=5.75, table=subpoc_table) 
#       get_xsec(Qsq=1.75, xbj=0.360, t=t0-tl, Ein=5.75, table=subpoc_table) 
         
