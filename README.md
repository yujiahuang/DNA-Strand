## Quickstart

Install Simulated Annealing first
```
pip install -e git+https://github.com/perrygeo/simanneal.git  # latest from github
```

Run
```
python dna-strand.py [-r/-n]
```
`-r` is to randomly generate initial strands, `-n` is to use NuPack's strands initial strands. The default is random

Code Description

We created a class `GenerateStrand` that inherits from `simanneal.Annealer`. State is a list of 3 domains (a, b, c for Amplifier case). Definitions of moves and costs(energy) are in the report.
The annealing parameters are:
```
steps = 60000
Tmin = 1.0
```

Code is on github: `https://github.com/yujiahuang/DNA-Strand`





