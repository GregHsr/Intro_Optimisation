import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy.special import gamma

"""
########## EXERCICE 1 ##########

#-------- data ----------------------------------------

t_obs=np.array([0., 100., 198., 306.])
y_obs=np.array([1. , 5., 48., 505. ])

#------------------------------------------------------

param_init=np.array([0])
def model(t, p0):
    return  p0*np.exp(+t/50) 

param_opt, pcov = optimize.curve_fit(model, t_obs, y_obs, param_init) 

print(param_opt)

plt.figure(1)
plt.plot(t_obs, y_obs, linestyle="None", linewidth=4.0, marker='.', ms=10, color='red')
plt.xlabel("temps [min]")
plt.ylabel("Population (en milliers)")
plt.grid()


plt.figure(2)
plt.plot(t_obs, model(t_obs, param_opt[0]), linestyle='-', linewidth=1.0, marker='o', ms=6, color='black',label="Données calculées")
plt.plot(t_obs, y_obs, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red', label="Données expérimentales")
plt.xlabel("temps [min]")
plt.ylabel("Population (en milliers)")
plt.grid()
plt.legend()
"""

########## EXERCICE 2 ##########

#-------- data ----------------------------------------

dt = 0.5

Q_obs = np.loadtxt("ALIOU93_Q.txt", skiprows=1)
P_obs = np.loadtxt("ALIOU93_P.txt", skiprows=1)

temps = [dt*k for k in range(0, len(Q_obs))]

temps_mod = temps[3551:3850]
P_obs_mod = P_obs[3551:3850]
Q_obs_mod = Q_obs[3551:3850]

#-------- Fonctions ----------------------------------

def residus(valeurs_calculees, Q_obs):
    resi = 0
    for k in range(len(Q_obs)):
        resi += (Q_obs[k] - valeurs_calculees[k])**2
    return (resi/len(valeurs_calculees))

#-------- Models -------------------------------------

## Modèle 1

Q0_1 = 1
k_1 = 1
a_1 = 1

def model_1(temps, k, a):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q0_1 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt, pcov = optimize.curve_fit(model_1, temps_mod, Q_obs_mod,(1,1))

print(param_opt)
print(residus(model_1(temps_mod, param_opt[0],param_opt[1]), Q_obs_mod))

plt.figure(3)
valeurs_calculees_1 = model_1(temps_mod, param_opt[0],param_opt[1])

plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_1, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()



## Modèle 2

# Métode 1

Q0_21 = Q_obs_mod[0]
k_21 = 1
a_21 = 1

def model_2(temps, k, a, Q0):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q0 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt21, pcov = optimize.curve_fit(model_2, temps_mod, Q_obs_mod,(k_21,a_21,Q0_21))
print("param_opt 1", param_opt21)
print("residus 1", residus(model_2(temps_mod, param_opt21[0],param_opt21[1],param_opt21[2]), Q_obs_mod))

# Métode 2

Q0_22 = 1
k_22 = 1
a_22 = 1

def model_2(temps, k, a, Q0):
    if Q0 < 0:
        pass
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q0 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt22, pcov = optimize.curve_fit(model_2, temps_mod, Q_obs_mod,(k_22,a_22,Q0_22))
print("param_opt 2", param_opt22)
print("residus 2", residus(model_2(temps_mod, param_opt22[0],param_opt22[1],param_opt22[2]), Q_obs_mod))

# Comparaison graphique

plt.figure(4)
plt.subplot(211)
valeurs_calculees_21 = model_2(temps_mod, param_opt21[0],param_opt21[1],param_opt21[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_21, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("Q0 = Q_obs_mod[0]")

plt.subplot(212)
valeurs_calculees_22 = model_2(temps_mod, param_opt22[0],param_opt22[1],param_opt22[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_22, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("Q0 > 0")



# Changement de la fenêtre d'observation

temps_mod2 = temps[1300:1500]
P_obs_mod2 = P_obs[1300:1500]
Q_obs_mod2 = Q_obs[1300:1500]


Q0_31 = Q_obs_mod2[0]
k_31 = 1
a_31 = 1

def model_31(temps, k, a):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod2[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q0_31 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt31, pcov = optimize.curve_fit(model_31, temps_mod2, Q_obs_mod2,(k_31,a_31))
print("param_opt 31", param_opt31)
print("residus 31", residus(model_31(temps_mod2, param_opt31[0],param_opt31[1]), Q_obs_mod2))

temps_mod3 = temps[1300:1800]
P_obs_mod3 = P_obs[1300:1800]
Q_obs_mod3 = Q_obs[1300:1800]

def model_32(temps, k, a):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod3[s]*np.exp(-k*(temps[t]-temps[s]))
        Retour.append(Q0_31 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt32, pcov = optimize.curve_fit(model_32, temps_mod3, Q_obs_mod3,(k_31,a_31))
print("param_opt 32", param_opt32)
print("residus 32", residus(model_32(temps_mod3, param_opt32[0],param_opt32[1]), Q_obs_mod3))


plt.figure(5)
plt.subplot(211)
valeurs_calculees_31 = model_31(temps_mod2, param_opt31[0],param_opt31[1])
plt.plot(temps_mod2, Q_obs_mod2, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod2, valeurs_calculees_31, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("t = [1300:1500]")

plt.subplot(212)
valeurs_calculees_32 = model_32(temps_mod3, param_opt32[0],param_opt32[1])
plt.plot(temps_mod3, Q_obs_mod3, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod3, valeurs_calculees_32, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("t = [1300:1800]")


# Modèle 4

Q0_4 = 1
k_4 = 0.104
a_4 = 0.114
tau_4 = 0

def model_4(temps, k, a, Q0):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]- tau_4))
        Retour.append(Q0 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt4, pcov = optimize.curve_fit(model_4, temps_mod, Q_obs_mod,(k_4,a_4,Q0_4))
print("param_opt 4", param_opt4)
print("residus 4", residus(model_4(temps_mod, param_opt4[0],param_opt4[1],param_opt4[2]), Q_obs_mod))

plt.figure(6)
plt.subplot(221)
valeurs_calculees_4 = model_4(temps_mod, param_opt4[0],param_opt4[1],param_opt4[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_4, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("tau = 0")

plt.subplot(222)
tau_4 = 5
valeurs_calculees_4 = model_4(temps_mod, param_opt4[0],param_opt4[1],param_opt4[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_4, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("tau = 5")

plt.subplot(223)
tau_4 = 10
valeurs_calculees_4 = model_4(temps_mod, param_opt4[0],param_opt4[1],param_opt4[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_4, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("tau = 10")

plt.subplot(224)
tau_4 = 20
valeurs_calculees_4 = model_4(temps_mod, param_opt4[0],param_opt4[1],param_opt4[2])
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_4, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("tau = 20")


# Comparaison avec le modèle 1

tau_4 = 3

def modele_41(temps, k, a, Q0, tau):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s]*np.exp(-k*(temps[t]-temps[s]+tau))
        Retour.append(Q0 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt41, pcov = optimize.curve_fit(modele_41, temps_mod, Q_obs_mod,(k_4,a_4,Q0_4,tau_4))
print("param_opt 41", param_opt41)

valeurs_calculees_41 = modele_41(temps_mod, param_opt41[0],param_opt41[1],param_opt41[2],param_opt41[3])
print("residus 41", residus(valeurs_calculees_41, Q_obs_mod))

plt.figure(7)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_1, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées (sans tau)")
plt.plot(temps_mod, valeurs_calculees_41, linestyle="None", linewidth=1.0, marker='+', ms=7, color='green',label="Données calculées (avec tau)")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()


# Modèle à plusieurs réservoirs

Q0_5 = 1
k_5 = 0.104
a_5 = 0.114
tau_5 = 0.441
M_5 = 1

def model_5(temps, k, a, tau, M):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs_mod[s] * np.exp(-k*(temps[t]-temps[s])) * (k*(temps[t]-temps[s]))**(M-1) / gamma(M)
        Retour.append(Q0_5 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

param_opt5, pcov = optimize.curve_fit(model_5, temps_mod, Q_obs_mod,(k_5,a_5,tau_5,M_5))
print("param_opt 5", param_opt5)

valeurs_calculees_5 = model_5(temps_mod, param_opt5[0],param_opt5[1],param_opt5[2],param_opt5[3])
print("residus 5", residus(valeurs_calculees_5, Q_obs_mod))

plt.figure(8)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_5, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()

# Tests avec plusieurs réservoirs

plt.figure(9)
plt.subplot(221)
M_51 = 1
valeurs_calculees_51 = model_5(temps_mod, param_opt5[0],param_opt5[1],param_opt5[2],M_51)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_51, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("M = 1")

plt.subplot(222)
M_52 = 5
valeurs_calculees_52 = model_5(temps_mod, param_opt5[0],param_opt5[1],param_opt5[2],M_52)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_52, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("M = 5")

plt.subplot(223)
M_53 = 10
valeurs_calculees_53 = model_5(temps_mod, param_opt5[0],param_opt5[1],param_opt5[2],M_53)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_53, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("M = 10")

plt.subplot(224)
M_54 = 20
valeurs_calculees_54 = model_5(temps_mod, param_opt5[0],param_opt5[1],param_opt5[2],M_54)
plt.plot(temps_mod, Q_obs_mod, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
plt.plot(temps_mod, valeurs_calculees_54, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()
plt.title("M = 20")

# Résultat final

def model_fin(temps, k, a, tau, M):
    Retour = []
    for t in range(len(temps)):
        integ = 0
        for s in range(t):
            integ += dt*P_obs[s] * np.exp(-k*(temps[t]-temps[s])) * (k*(temps[t]-temps[s]))**(M-1) / gamma(M)
        Retour.append(Q0_5 * np.exp(-k*temps[t]) + a*integ)
    return  Retour

plt.figure(10)
plt.plot(temps, Q_obs, linestyle="None", linewidth=1.0, marker='o', ms=4, color='red',label="Données expérimentales")
valeurs_calculees_10 = model_fin(temps, param_opt5[0],param_opt5[1],param_opt5[2],param_opt5[3])
plt.plot(temps, valeurs_calculees_10, linestyle="None", linewidth=1.0, marker='+', ms=7, color='blue',label="Données calculées")
plt.ylabel("Débit (m3/s)")
plt.xlabel("Temps (heure)")
plt.legend()
plt.grid()

plt.show()
