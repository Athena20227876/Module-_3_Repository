# GLP1_GLP1R.py

import numpy as np
import matplotlib.pyplot as plt

def Fact(G, S, Kd1, Kd2):
    """Fractional receptor activation: (Kd2*G + Kd1*S) / (Kd1*Kd2 + Kd2*G + Kd1*S)"""
    return (Kd2*G + Kd1*S) / (Kd1*Kd2 + Kd2*G + Kd1*S)

Kd1 = 0.17  # GLP1 Kd [nM]
Kd2 = 0.38  # Semaglutide Kd [nM]
G_fasting = 0.015 # Fasting [GLP1](nM)
G_fed = 0.04  # After meal [GLP1](nM)

S = np.logspace(-3, 3, 500)     # ligand concentration, 0.001 to 1000 [nM]

plt.figure()
plt.semilogx(S, Fact(G_fasting, S, Kd1, Kd2), lw=2, label="Fasting GLP1")
plt.semilogx(S, Fact(G_fed, S, Kd1, Kd2), lw=2, label="Fed GLP1")
plt.axvline(1.0, ls="--", c="gray", lw=1)
plt.ylim(0, 1)
plt.xlabel("[Semaglutide] (nM, log scale)")
plt.ylabel("Fraction active (Fact)")
plt.legend()
plt.savefig("GLP1_GLP1R.png", dpi=200)
plt.show()


# Interpretation: With no drug, eating roughly doubles receptor activation, but with the drug at a typical 1nM semaglutide level, activation is 
# at a similat level of around 73%, showing that semaglutide drives most of the receptor activation, regardless of whether the person has eaten or not.