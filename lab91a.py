import numpy as np
import matplotlib.pyplot as plt

# --- 1. Définition des Paramètres du Circuit ---
# Valeurs nominales (inchangées)
R = 1500  # Résistance en Ohms (1.5 kΩ)
C = 0.22e-6  # Capacité en Farads (0.22 μF)
V_g_peak_to_peak = 8  # Tension crête à crête de la source (V)

# --- 2. Analyse Théorique ---

# Fréquence de coupure théorique (en radians/seconde)
omega_c = 1 / (R * C)
# Fréquence de coupure théorique (en Hertz)
f_c = omega_c / (2 * np.pi)

# --- MODIFICATION 1: Ajustement de la plage de fréquences pour centrer f_c ---
# Pour centrer f_c sur une échelle logarithmique, on prend un intervalle symétrique en décades,
# par exemple, de 2 décades en dessous à 2 décades au-dessus.
f_min = 0.01 * f_c  # 2 décades en dessous de f_c
f_max = 100 * f_c   # 2 décades au-dessus de f_c
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

# Valeur du module à la fréquence de coupure (-3 dB)
valeur_a_coupe_dB = -3.0

# --- MODIFICATION 2: Calcul de l'asymptote ---

# Asymptote basse fréquence (f << f_c) : Pente de +20 dB/décade, passe par 0 dB à f_c.
# A_LF(f) = 20 * log10(f/f_c)
asymptote_lf_dB = 20 * np.log10(f / f_c)

# Asymptote haute fréquence (f >> f_c) : 0 dB.
asymptote_hf_dB = np.zeros_like(f)

# L'asymptote totale est la ligne brisée :
# - Pente de +20 dB/décade pour f < f_c
# - Pente de 0 dB pour f > f_c
asymptote_dB = np.where(f < f_c, asymptote_lf_dB, asymptote_hf_dB)


# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(12, 10))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
plt.semilogx(f, module_H_dB, color='blue', label='Module réel $|H(j\omega)|$')
# Ajout de l'asymptote
plt.semilogx(f, asymptote_dB, color='black', linestyle='-', linewidth=1.5, label='Asymptote (Approximation)')

plt.title('Diagramme de Bode - Filtre Passe-Haut RC')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--", alpha=0.7)

# Ajout de la fréquence de coupure sur le module
plt.axvline(x=f_c, color='red', linestyle='--', label=f'Fréquence de coupure $f_c$ = {f_c:.2f} Hz')
plt.axhline(y=valeur_a_coupe_dB, color='green', linestyle=':', label='-3 dB')
plt.axhline(y=0, color='gray', linestyle=':', label='Asymptote haute fréquence (0 dB)') # Ligne de 0 dB
plt.legend(loc='lower left')
plt.ylim(np.min(module_H_dB) - 5, 5) # Ajuster les limites pour mieux voir l'asymptote basse fréquence

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange', label='Phase réelle $\\angle H(j\omega)$')
plt.title('Diagramme de Bode - Phase')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--", alpha=0.7)

# Ajout de la fréquence de coupure sur la phase
plt.axvline(x=f_c, color='red', linestyle='--', label=f'$f_c$ = {f_c:.2f} Hz')
plt.axhline(y=45, color='green', linestyle=':', label='$45^{\\circ}$ (à $f_c$)')

# Asymptotes de phase: 90 degrés à basse fréquence, 0 degré à haute fréquence.
plt.axhline(y=90, color='gray', linestyle=':', label='Asymptote basse fréquence ($90^{\\circ}$)')
plt.axhline(y=0, color='gray', linestyle=':', label='Asymptote haute fréquence ($0^{\\circ}$)')

plt.legend(loc='lower left')
plt.yticks(np.arange(0, 91, 15)) # Ajuster les ticks pour mieux visualiser

plt.tight_layout()
# plt.savefig('diagramme_de_bode_filtre_rc_modifie.png')

plt.show()