import re
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt

random_instance_regex = re.compile(r'instances/p\d+_\d+\.dat')

instances = [f for f in glob.glob('instances/*.dat')]
randoms = [f for f in instances if re.fullmatch(random_instance_regex, f)]
domains = [f for f in instances if not re.fullmatch(random_instance_regex, f)]

data = []

for name in instances:
    bname = os.path.basename(name)
    with open(name) as f:
        n,m = map(int,f.readline().split())
        degs = [int(line.split(maxsplit=1)[0]) for line in f]
    if sum(degs)<100:
        continue

    stype = not bool(re.fullmatch(random_instance_regex, name))
    stype += bname.startswith('cost_matrix')
    data.append([m/n if m>n else -n/m, sum(degs)/m, stype])

    info = f'{n=:>7} {m=:>7} size={sum(degs):>9} {bname}'
    print(info)

    if stype!=1:
        continue
    plt.step(sorted(degs)[::-1], range(m))
    plt.xlabel('deg')
    plt.ylabel('#edges with >= degree')
    plt.title(info[:50])
    plt.show()
    plt.clf()

df = pd.DataFrame(data, columns=['ratio', 'avg_deg', 'stype'])
print(df)

colors = {True: 'blue', False: 'red'}
df[df.stype==0].plot.scatter('avg_deg', 'ratio', c='blue', label='random')
df[df.stype==1].plot.scatter('avg_deg', 'ratio', c='red', label='real', ax=plt.gca())
df[df.stype==2].plot.scatter('avg_deg', 'ratio', c='green', label='cost_matrix', ax=plt.gca())
plt.semilogx()
plt.yscale('symlog')
plt.savefig('summary.pdf')
plt.show()


