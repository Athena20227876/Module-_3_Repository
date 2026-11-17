# Find fit for Cavalitide binding data

# Step 1: Import packages needed for Task 1
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Step 2: Create model
def Fact(L1, L2, Kd1, Kd2):
    """Fraction of active receptors with two agonists competing for one receptor.
    L1, Kd1 = GLP1 (natual ligand); L2, Kd2 = Cavalitide (drug being tested)"""
    return (Kd2*L1 + Kd1*L2)/(Kd1*Kd2 + Kd2*L1 + Kd1*L2) 

# Step 3: State known parameters
Kd1 = .17                                                               # GLP1 dissociation constant [nM]
L1 = .015                                                               # Concentration of the first ligand, GLP1 after fasting [nM]
L2 = np.logspace(-3, 2 ,500)                                            # Array of 500 Cavalitide concentrations between .001 and 100 nM, smoothes model curve
# Kd2 = variable we're solving for

# Step 4: Load the data into the file

data = np.loadtxt("cavalitide_data.csv", delimiter=",", skiprows=1)     # Load data set, comma separations, skip first header row
L2_loaded, Fact_loaded = data[:,0], data[:,1]                           # Separates L2 (Cavalitide concentrations) and active fraction of drug

# Step 5: Approach 1 - brute force optimization of Kd2

Kd2_candidates = np.logspace(-3, 2, 1000)                               # 1000 candidate Kd values between .001 and 100 nM
SSE_brute = np.inf                                                      # starts best error at infinity
for Kd2_test in Kd2_candidates:
    SSE = np.sum((Fact_loaded - Fact(L1, L2_loaded, Kd1, Kd2_test))**2) # Sum of squared errors formula
    if SSE < SSE_brute:                                                 # If new SSE is lower than initial SSE, that is saved with its Kd
        SSE_brute = SSE
        Kd2_brute = Kd2_test                                            # Lowest SSE is saved into Kd2_brute
print(f"Brute-force Kd2 = {Kd2_brute:.4f} nM (SSE = {SSE_brute:.4f})")  # Prints lowest SSE Kd

# Step 6: Approach 2 - curve fit of Kd2
popt, pcov = curve_fit(lambda L2, Kd2: Fact(L1, L2, Kd1, Kd2),          # L2_loaded is x, Fact_loaded is y
                       L2_loaded, Fact_loaded, p0=[.1], bounds=(0, np.inf))
Kd2_fit = popt[0]
Kd2_error = np.sqrt(pcov[0,0])                                          # Find error of measurement
SSE_fit = np.sum(Fact_loaded - Fact(L1, L2_loaded, Kd1, Kd2_fit)**2)
print(f"Curve_fit Kd2 = {Kd2_fit:.4f} +/- {Kd2_error:.4f} nM (SSE = {SSE_fit:.4f})")

# Step 7: Compare data from Cavalitide Kd with Semaglutide Kd

Kd_semabrute = .335                                                     # Values obtained from semaglutidePK.py file Approach 1
Kd_semabrute_SSE = .0047                                                # Values obtained from semaglutidePK.py file Approach 1
Cav_vs_sema_brute = (Kd_semabrute / Kd2_fit)

Kd_semacurvefit = .350                                                  # Values obtained from semaglutidePK.py file Approach 2
Kd_semacurvefit_SSE = .0044                                             # Values obtained from semaglutidePK.py file Approach 2
Cav_vs_sema_curvefit = (Kd_semacurvefit / Kd2_fit)

Cav_vs_sema_avg = (Cav_vs_sema_brute + Cav_vs_sema_curvefit) / 2        # Take and print average of two approaches
print(f"Cavalitide binds ~{Cav_vs_sema_avg:.2f} times more tightly than semaglutide")

# Step 8: Plot data and model fits to create figure for data

plt.figure()                                                            # Create empty figure
plt.semilogx(L2_loaded, Fact_loaded, "o", c="tab:purple",
            label="Experimental data")                                  # Plots measured points with circle markers
plt.semilogx(L2, Fact(L1, L2, Kd1, Kd2_brute), ":", lw=4,
             label=f"Brute-Force Fit (Kd2 = {Kd2_brute:.4f} nM)")       # Plots model curve with brute-force Kd
plt.semilogx(L2, Fact(L1, L2, Kd1, Kd2_fit), "--", lw=2,
             label=f"Curve-Fit Fit (Kd2 = {Kd2_fit:.4f} nM)")           # Plots model curve with curve-fit Kd
plt.ylim(-.1, 1.1)
plt.xlabel("[Cavalitide] (nM)")
plt.ylabel("Active Fraction")
plt.title("GLP1R Activation by Cavalitide with GLP1 Present")
plt.legend(loc="upper left")
plt.show()                                                              # Create plot with above criteria