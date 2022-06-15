
import subprocess
import json
import sys
import random
import pandas as pd
import matplotlib.pyplot as plt

random.seed(378)

GENERATOR_CLI = '../build/genbip'
SOLVER_CLI = '/home/chris/repos/findminhs/target/release/findminhs'
header = ['n','m','deg','T','univar','unicls','seed','opt','steps','time']

def measure(n=100, m=100, deg=10, T=0.0, univar=0, unicls=1):
    seed = random.randint(0,10**9)
    genCall = [GENERATOR_CLI, *map(str, [n, m, deg, seed, T, univar, unicls])]
    solCall = [SOLVER_CLI, 'solve', '-r', 'report.json', 'bigirg.graph', 'settings.json']
    with open('bigirg.graph', 'w') as f:
        genlog = subprocess.run(genCall, stdout=f, stderr=subprocess.PIPE)
    subprocess.run(solCall, stderr=subprocess.DEVNULL)
    with open('report.json', 'r') as f:
        report = json.load(f)
    opt = report['opt']
    steps = report['branching_steps']
    rtime = report['runtimes']['total']
    print(f"T {T: >7}   OPT {opt: >7}   STEPS {steps: >7}   TIME {rtime: >7.4f}", genlog.stderr)
    return [n,m,deg,T,univar,unicls,seed,opt,steps,rtime]

def singleSpread(reps=5, **kwargs):
    return [measure(**kwargs) for _ in range(reps)]


data = []
for deg in range(5, 31,2):
    print(deg)
    data += singleSpread(reps=25, deg=deg, univar=0, unicls=1, n=3000, m=3000)
df = pd.DataFrame(data, columns=header)
df.boxplot(column='steps', by='deg')
plt.show()
