import numpy as np
import matplotlib.pyplot as plt

# --- 1. Définition des Paramètres du Circuit ---
R = 1500  # Résistance en Ohms (1.5 kΩ)
C = 0.22e-6  # Capacité en Farads (0.22 μF)

# --- 2. Analyse Théorique ---
omega_c = 1 / (R * C)
f_c = omega_c / (2 * np.pi)

# Définir la plage de fréquences (centrée sur f_c pour le tracé)
f_min = 0.01 * f_c  
f_max = 100 * f_c   
f = np.logspace(np.log10(f_min), np.log10(f_max), 500)
omega = 2 * np.pi * f

# Calcul de H(jω)
H = (1j * omega * R * C) / (1 + 1j * omega * R * C)
module_H_dB = 20 * np.log10(np.abs(H))
phase_H = np.angle(H, deg=True) 
valeur_a_coupe_dB = -3.0

# Calcul de l'asymptote du module (déjà présente)
asymptote_lf_dB = 20 * np.log10(f / f_c)
asymptote_hf_dB = np.zeros_like(f)
asymptote_dB = np.where(f < f_c, asymptote_lf_dB, asymptote_hf_dB)


# --- Calcul de l'asymptote de la phase ---
f1_phase = 0.1 * f_c  # 1 décade sous fc
f2_phase = 10 * f_c   # 1 décade au-dessus de fc

# Région centrale (entre f1 et f2) : 90 - 45 * log10(f / f1)
phase_asymptote = 90 - 45 * (np.log10(f) - np.log10(f1_phase))

# Région 1: f < f1_phase -> 90 degrés
phase_asymptote = np.where(f < f1_phase, 90, phase_asymptote)

# Région 3: f > f2_phase -> 0 degrés
phase_asymptote = np.where(f > f2_phase, 0, phase_asymptote)


# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(12, 10))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
# Utilisation de chaînes brutes r'' pour éviter les avertissements LaTeX
plt.semilogx(f, module_H_dB, color='blue', label=r'Module $|H(j\omega)|$')
plt.semilogx(f, asymptote_dB, color='blue', linestyle='--', linewidth=1.5, label='Asymptotes')

plt.title('Diagramme de Bode module circuit RC passe-haut')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--", alpha=0.7)

plt.axvline(x=f_c, color='red', linestyle='--', label=f'Fréquence de coupure $f_c$ = {f_c:.2f} Hz')
plt.axhline(y=valeur_a_coupe_dB, color='green', linestyle=':', label='-3 dB')
plt.legend(loc='lower left')
plt.ylim(np.min(module_H_dB) - 5, 5)

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange', label=r'Phase $\angle H(j\omega)$')

# Ajout de l'asymptote de la phase (la ligne noire brisée)
plt.semilogx(f, phase_asymptote, color='orange', linestyle='--', linewidth=1.5, label='Asymptote')

plt.title('Diagramme de Bode phase circuit RC passe-haut')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--", alpha=0.7)

plt.axvline(x=f_c, color='red', linestyle='--', label=f'$f_c$ = {f_c:.2f} Hz')
plt.axhline(y=45, color='green', linestyle=':', label=r'$45^{\circ}$ (à $f_c$)')

# Lignes de référence des plateaux

plt.legend(loc='lower left')
plt.yticks(np.arange(0, 91, 15))

plt.tight_layout()
plt.show() # Pour maintenir la fenêtre ouverte