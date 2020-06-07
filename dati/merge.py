import pandas as pd

a = 2
b = 3

hA = pd.read_csv(f"histo/h{a}.dat", sep=" ", header=None)
hB = pd.read_csv(f"histoB/hB{b}.dat", sep = " ", header=None)

fA = pd.read_csv(f"fract/f{a}.dat", sep=" ", header=None)
fB = pd.read_csv(f"fractB/fB{b}.dat", sep = " ", header=None)

h = hA + hB
f = fA + fB

run_nuovo = 3

h.to_csv(f"histo/h{run_nuovo}.dat", sep=" ", index=False, header=False)
f.to_csv(f"fract/f{run_nuovo}.dat", sep=" ", index=False, header=False)


