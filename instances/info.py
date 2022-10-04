import re
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle

def coefficientOfVariation(lst):
    return np.std(lst) / np.mean(lst)

def collectInstanceStats():
    random_instance_regex = re.compile(r'instances/p\d+_\d+\.dat')
    instances = [f for f in glob.glob('instances/*.dat')]
    randoms = [f for f in instances if re.fullmatch(random_instance_regex, f)]
    domains = [f for f in instances if not re.fullmatch(random_instance_regex, f)]

    data = []
    for name in instances:
        bname = os.path.basename(name)
        with open(name) as f:
            n,m = map(int,f.readline().split())
            edges = [list(map(int,line.split()[1:])) for line in f]
            edgedegs = [len(edge) for edge in edges]
            nodedegs = [0]*n
            for edge in edges:
                for node in edge:
                    nodedegs[node] += 1

        stype = not bool(re.fullmatch(random_instance_regex, name))
        stype += bname.startswith('cost_matrix')

        assert(sum(nodedegs) == sum(edgedegs))

        data.append([
            m/n if m>n else -n/m, 
            stype, 
            n,
            m,
            sum(nodedegs),
            sum(nodedegs)/n, 
            sum(edgedegs)/m, 
            np.var(nodedegs), 
            np.var(edgedegs),
            coefficientOfVariation(nodedegs),
            coefficientOfVariation(edgedegs),
        ])

        info = f'{n=:>7} {m=:>7} size={sum(edgedegs):>9} {bname}'
        print(info)

        continue # after this is some plotting stuff
        plt.step(sorted(degs)[::-1], range(m))
        plt.xlabel('deg')
        plt.ylabel('#edges with >= degree')
        plt.title(info[:50])
        plt.show()
        plt.clf()
    return data

summary_file = 'summary.pkl'
if os.path.isfile(summary_file):
    with open(summary_file, 'rb') as f:
        data = pickle.load(f)
else:
    data = collectInstanceStats()
    with open(summary_file, 'wb') as f:
        pickle.dump(data, f)

df = pd.DataFrame(data, columns=['ratio', 'stype', 'n', 'm', 'size', 'node_avg_deg', 'edge_avg_deg', 'node_var', 'edge_var', 'node_variation', 'edge_variation'])
#print(df)

relevant = df[(df['size']>100) & (df.stype==1)]
def plotScatter(data, x,y,name):
    data[data.stype==0].plot.scatter(x, y, c='blue', label='random')
    data[data.stype==1].plot.scatter(x, y, c='red', label='real', ax=plt.gca())
    #data[data.stype==2].plot.scatter(x, y, c='green', label='cost_matrix', ax=plt.gca())
    plt.semilogx()
    #plt.yscale('symlog')
    plt.savefig(f'{name}.pdf')
    plt.clf()
    #plt.show()

print(relevant.node_variation.describe())
print(relevant.edge_variation.describe())

#plotScatter(df,'node_avg_deg', 'node_var', 'node_dist')
#plotScatter(df, 'edge_avg_deg', 'edge_var', 'edge_dist')
plotScatter(df, 'size', 'node_variation', 'node_var')
plotScatter(df, 'size', 'edge_variation', 'edge_var')

plt.clf()
relevant.node_variation.hist(bins=1000, density=1, cumulative=True, histtype='step')
relevant.edge_variation.hist(bins=1000, density=1, cumulative=True, histtype='step', ax=plt.gca())
plt.show()


