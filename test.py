
import math
import itertools 
from fractions import Fraction as frac

p = [frac(1,2), frac(1,2), frac(1,4)]

nullevt = frac(1,1)
for evt in p:
    nullevt *= 1-evt

event = {}
for k in range(1,len(p)+1):
    for X in itertools.combinations(range(len(p)), k):
        ans = frac(1,1)
        for i in range(len(p)):
            ans*=(p[i] if i in X else 1-p[i])
        event[X] = ans / (1-nullevt)

for e in event:
    print(e, event[e])

scale = [2,2,4]
bi = [frac(0,1)]*len(p)
for k in range(1,len(p)+1):
    for X in itertools.combinations(range(len(p)), k):
        for i in range(len(p)):
            if i in X:
                bi[i] += event[X] * scale[i] / sum(scale[j] for j in X)
print(bi)


