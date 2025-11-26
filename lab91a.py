import numpy as np
import matplotlib.pyplot as plt
# --- 1. Définition des Paramètres du Circuit ---
# Valeurs nominales
R = 1500  # Résistance en Ohms (1.5 kΩ)
C = 0.22e-6  # Capacité en Farads (0.22 μF)
V_g_peak_to_peak = 8  # Tension crête à crête de la source (V)

# --- 2. Analyse Théorique ---

# i) Expression Analytique de la Fonction de Transfert H(jω) = V_R / V_g
# C'est un filtre passe-haut puisque la sortie est prise aux bornes de la résistance R.
# H(jω) = R / (R + 1/(jωC)) = (jωRC) / (1 + jωRC)

# Fréquence de coupure théorique (en radians/seconde)
omega_c = 1 / (R * C)
# Fréquence de coupure théorique (en Hertz)
f_c = omega_c / (2 * np.pi)

print(f"**Fréquence de coupure théorique (f_c):** {f_c:.2f} Hz\n")

# Définir la plage de fréquences (échelle logarithmique pour le tracé)
f_min = 1
f_max = 10 * f_c  # Allons jusqu'à 10 fois la fréquence de coupure pour une bonne visualisation
f = np.logspace(np.log10(f_min), np.log10(f_max), 500)
omega = 2 * np.pi * f

# Calcul de l'expression complexe de la fonction de transfert H(jω)
# H(jω) = (j * omega * R * C) / (1 + j * omega * R * C)
H = (1j * omega * R * C) / (1 + 1j * omega * R * C)

# Calcul du module et de la phase
module_H = np.abs(H)
phase_H = np.angle(H, deg=True)  # Phase en degrés

# Calcul du module en décibel (dB)
module_H_dB = 20 * np.log10(module_H)

# iii) Valeur du module à haute fréquence (f → ∞)
# Pour un passe-haut, |H(jω)| → 1 lorsque ω → ∞.
# La valeur en dB à haute fréquence est donc : 20 * log10(1) = 0 dB.
# La fréquence de coupure est le point où le module diminue de -3 dB par rapport à cette valeur,
# soit -3 dB.
valeur_a_coupe_dB = 0 - 3

# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(10, 8))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
plt.semilogx(f, module_H_dB, color='blue')
plt.title('Diagramme de Bode - Module $\\left(|V_R / V_g|\\right)$')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--")

# Ajout de la fréquence de coupure sur le module
plt.axvline(x=f_c, color='red', linestyle='--', label=f'Fréquence de coupure $f_c$ = {f_c:.2f} Hz')
plt.axhline(y=valeur_a_coupe_dB, color='green', linestyle=':', label='-3 dB')
plt.legend()

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange')
plt.title('Diagramme de Bode - Phase $\\left(\\phi_{V_R} - \\phi_{V_g}\\right)$')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--")

# Ajout de la fréquence de coupure sur la phase
plt.axvline(x=f_c, color='red', linestyle='--', label=f'Fréquence de coupure $f_c$ = {f_c:.2f} Hz')
plt.axhline(y=45, color='green', linestyle=':', label='$45^{\\circ}$')
plt.legend()

plt.tight_layout()
plt.show()