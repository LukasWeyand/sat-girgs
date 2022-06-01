
import subprocess
import json
import sys
import random
import pandas as pd
import matplotlib.pyplot as plt

random.seed(378)

GENERATOR_CLI = '../build/genbip'
SOLVER_CLI = '/home/chris/repos/findminhs/target/release/findminhs'
data = []
header = ['n','m','deg','T','univar','unicls','seed','opt','steps','time']

def measure(n=100, m=100, deg=10, T=0.0, univar=1, unicls=0):
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
    data.append([n,m,deg,T,univar,unicls,seed,opt,steps,rtime])
    return rtime

def singleSpread(reps=5, **kwargs):
    for i in range(reps):
        measure(**kwargs)

def expdist():
    for univar in [0,1]:
        for unicls in [0,1]:
            print(f'{univar=} {unicls=}')
            singleSpread(reps=20, univar=univar, unicls=unicls)
    df = pd.DataFrame(data, columns=header)
    for col in ['opt','steps','time']:
        df.boxplot(column=col, by=['univar', 'unicls'])
        plt.savefig(f'plots/dist-{col}.pdf')

def expratio():
    for n in range(50,201,10):
        singleSpread(reps=20, n=n)
    df = pd.DataFrame(data, columns=header)
    for col in ['opt','steps','time']:
        df.boxplot(column=col, by='n')
        plt.savefig(f'plots/ratio-{col}.pdf')

def exptemp():
    for T in range(0, 10):
        singleSpread(reps=20, T=T/10)
    df = pd.DataFrame(data, columns=header)
    for col in ['opt','steps','time']:
        df.boxplot(column=col, by='T')
        plt.savefig(f'plots/temp-{col}.pdf')

def expdeg():
    for deg in range(5, 21):
        print(deg)
        singleSpread(reps=50, deg=deg)
    df = pd.DataFrame(data, columns=header)
    for col in ['opt','steps','time']:
        df.boxplot(column=col, by='deg')
        plt.savefig(f'plots/deg-{col}.pdf')

def doall():
    for univar in [0,1]:
        for unicls in [0,1]:
            if univar and unicls:
                continue

            # ratio
            data = []
            for n in range(50,201,10):
                singleSpread(reps=20, n=n, univar=univar, unicls=unicls)
            df = pd.DataFrame(data, columns=header)
            for col in ['opt','steps','time']:
                df.boxplot(column=col, by='n')
                plt.savefig(f'plots/{univar}{unicls}-ratio-{col}.pdf')
            
            # temp
            data = []
            for T in range(0, 10):
                singleSpread(reps=20, T=T/10, univar=univar, unicls=unicls)
            df = pd.DataFrame(data, columns=header)
            for col in ['opt','steps','time']:
                df.boxplot(column=col, by='T')
                plt.savefig(f'plots/{univar}{unicls}-temp-{col}.pdf')

            # deg
            data = []
            for deg in range(5, 21):
                print(deg)
                singleSpread(reps=50, deg=deg, univar=univar, unicls=unicls)
            df = pd.DataFrame(data, columns=header)
            for col in ['opt','steps','time']:
                df.boxplot(column=col, by='deg')
                plt.savefig(f'plots/{univar}{unicls}-deg-{col}.pdf')

