import pandas as pd
import matplotlib.pyplot as plt

free = pd.read_csv("en/en2.dat", sep = " ", header = None)
p = pd.read_csv("prob/pdn2.dat", sep = " ", header = None)
p = p.iloc[:, 1]


plt.plot(p)
plt.show()

en = np.log(0.0473)*np.arange(p.shape[0])-np.log(p/p[0])
plt.plot(en)

