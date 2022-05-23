
import subprocess
import json
import sys

def measure(T, seed=378):
    with open('bigirg.graph', 'w') as f:
        subprocess.run(['./genbip', '100', '1000', '10', str(seed), str(T)], stdout=f)
    subprocess.run(['/home/chris/repos/findminhs/target/release/findminhs', 
        'solve', '-r', 'report.json', 'bigirg.graph', 'settings.json'], stderr=subprocess.DEVNULL)
    with open('report.json', 'r') as f:
        report = json.load(f)
    
    runtime = float(report['runtimes']['total'])
    print(f"T {T: >7}   OPT {report['opt']: >7}   STEPS {report['branching_steps']: >7}   TIME {runtime: >7.4f}")
    return runtime


times = []
for T in range(0, 10):
    reps = 5
    avg = 0.0
    for i in range(reps):
        avg += measure(T/10, 1000*T+i+378)
    times.append(avg / reps)

print(*times, sep='\n')

