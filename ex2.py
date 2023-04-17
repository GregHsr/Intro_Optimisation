import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimize

# Question 1
# Lecture du fichier

dt = 0.5

Q_obs = np.loadtxt("ALIOU93_Q.txt", skiprows=1)
P_obs = np.loadtxt("ALIOU93_P.txt", skiprows=1)

temps = [dt*k for k in range(0, len(Q_obs))]

# Tracer les données

# plt.figure(1)
# # Echelle 1
# plt.plot(temps, Q_obs, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red')
# plt.ylabel("Débit (m3/s)")

# # Echelle 2
# plt.twinx()
# plt.plot(temps, P_obs, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue')
# plt.ylabel("Pluie (mm/0.5h)")

# plt.grid()
# plt.xlabel("Temps (heure)")

# plt.show()

# Question 2

k = 1
a = 1

temps_mod = temps[3551:3850]
P_obs_mod = P_obs[3551:3850]
Q_obs_mod = Q_obs[3551:3850]

def model(temps, k, a):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q_obs[0] * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt, pcov = optimize.curve_fit(model, temps_mod, Q_obs_mod, (k,a))
print(param_opt)

plt.figure(2)
valeurs_calculees = model(temps_mod, param_opt[0], param_opt[1])

plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.show()

# Calcul des résidus

def residus(param, valeurs_calculees, Q_obs):
    resi = 0
    for k in range(len(Q_obs)):
        resi += (Q_obs[k] - valeurs_calculees[k])**2
    return resi/len(valeurs_calculees)

print(residus(param_opt, valeurs_calculees, Q_obs_mod))

# plt.figure(3)

# valeurs_tot = model(temps, param_opt[0], param_opt[1])

# plt.plot(temps, Q_obs, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
# plt.plot(temps, valeurs_tot, linestyle="None", linewidth=0.3, marker='.', ms=3, color='blue',label="Données calculées")
# plt.ylabel("Débit (m3/s)")
# plt.xlabel("Temps (heure)")
# plt.legend()
# plt.show()

