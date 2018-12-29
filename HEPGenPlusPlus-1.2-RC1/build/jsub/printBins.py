#!/usr/bin/python3

import os, numpy

qsq_start = 1.0
qsq_end = 12.0
qsq_bins = 110

w_start = 1.0
w_end = 4.0
w_bins = 150

QSQ_LIST="QSQ "
W_LIST="W "

#execute the commands
for qsq_val in numpy.linspace(qsq_start, qsq_end, qsq_bins, endpoint=False):
  for w_val in numpy.linspace(w_start, w_end, w_bins, endpoint=False):
    #os.system("~/batchJobGentoo.sh /sc/userdata/regali/slcBatchHepGen/hepgen_EXP/install/bin/gkCrossSection "+str(qsq_val)+" "+str(w_val)+" -0.1 0 -a")
    #print("gkCrossSection {:.3f} {:.3f} -0.1 0 -a".format(qsq_val, w_val))
    print("{:.4f} {:.4f}".format(qsq_val, w_val))

#create the cache folder index

for qsq_val in numpy.linspace(qsq_start, qsq_end, qsq_bins, endpoint=False):
  QSQ_LIST = QSQ_LIST+ " {:.4f}".format(qsq_val)
for w_val in numpy.linspace(w_start, w_end, w_bins, endpoint=False):
  W_LIST = W_LIST+" {:.4f}".format(w_val)


DATAFOLDER="DATAFOLDER ./output"
with open('tableIndex', 'w') as ff:
  ff.write(QSQ_LIST)
  ff.write('\n')
  ff.write(W_LIST)
  ff.write('\n')
  ff.write(DATAFOLDER)
  ff.write('\n')
  ff.write('\n')

#print("--- resulting datacache indexFile ---")
#print("--- Mind: export HEPGEN_CREATESPA to the folder where the amplitudes should be saved beforehand!!!")
