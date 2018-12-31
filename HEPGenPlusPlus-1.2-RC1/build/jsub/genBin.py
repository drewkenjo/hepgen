import random

m = 0.93827203
xb = -1
while xb<0.1 or xb>1 or t0>0.1:
  q2 = random.uniform(1,10)
  w = random.uniform(1,4)

  xb = q2/(w*w - m*m + q2)
  t0 = xb*xb*m*m/(1-xb)

print('{:.4f} {:.4f}'.format(q2, w))
