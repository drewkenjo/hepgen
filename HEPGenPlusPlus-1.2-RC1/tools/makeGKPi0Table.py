#!/usr/bin/python3

import os

qsq_start = 1.0
qsq_end = 26.0
qsq_bins = 100

w_start = 1.0
w_end = 21.0
w_bins = 100

QSQ_LIST="QSQ "
W_LIST="W "
DATAFOLDER="DATAFOLDER ./output"

#execute the commands
for qsq in range(0,qsq_bins):
  qsq_val = qsq*(qsq_end-qsq_start)/qsq_bins+qsq_start
  for w in range(0,w_bins):
    w_val = w*(w_end-w_start)/w_bins+w_start
    #os.system("~/batchJobGentoo.sh /sc/userdata/regali/slcBatchHepGen/hepgen_EXP/install/bin/gkCrossSection "+str(qsq_val)+" "+str(w_val)+" -0.1 0 -a")
    print("gkCrossSection "+str(qsq_val)+" "+str(w_val)+" -0.1 0 -a")

#create the cache folder index

for qsq in range(0,qsq_bins):
  qsq_val = qsq*(qsq_end-qsq_start)/qsq_bins+qsq_start
  QSQ_LIST = QSQ_LIST+ " "+str(qsq_val)
for w in range(0,w_bins):
    w_val = w*(w_end-w_start)/w_bins+w_start
    W_LIST = W_LIST+" " + str(w_val)


print("--- resulting datacache indexFile ---")
print(QSQ_LIST)
print(W_LIST)
print(DATAFOLDER)

print("--- Mind: export HEPGEN_CREATESPA to the folder where the amplitudes should be saved beforehand!!!")
