import numpy as np
import matplotlib.pyplot as plt

# --- 1. Paramètres du Circuit RLC (inchangés, R est la résistance interne RL) ---
L = 36.28e-3  
C = 0.1e-6 
R = 100    # Résistance interne RL de l'inductance, utilisée pour le calcul du Q

# --- 2. Analyse Théorique ---
omega_0 = 1 / np.sqrt(L * C)
f_0 = omega_0 / (2 * np.pi) # Fréquence de résonance fn ou fr (selon la définition)

# Le facteur de qualité est calculé avec R = RL
Q = (1/R) * np.sqrt(L/C) 

# --- NOUVEAU: Le gain max asymptotique est 0 dB pour le passe-bas (H=Vc/Vg) ---
gain_asymptote_fo = 0 
print(f"**Fréquence de résonance théorique (f_0):** {f_0:.2f} Hz")
print(f"**Facteur de qualité (Q):** {Q:.2f}")
print(f"**Gain asymptotique à f0 (G_croisement):** {gain_asymptote_fo:.2f} dB\n")

# Plage de fréquences
f_min = f_0 / 100 
f_max = f_0 * 100 
f = np.logspace(np.log10(f_min), np.log10(f_max), 500)
omega = 2 * np.pi * f

# --- CORRECTION DE LA FONCTION DE TRANSFERT (H = Vc / Vg) ---
# H(jω) = (1/(jωC)) / (R + jωL + 1/(jωC))
Z_total = R + 1j * omega * L + (1 / (1j * omega * C))
Z_C = 1 / (1j * omega * C)
H = Z_C / Z_total
module_H_dB = 20 * np.log10(np.abs(H))
phase_H = np.angle(H, deg=True) 

# --- CORRECTION DES ASYMPTOTES pour le PASSE-BAS (Ordre 2) ---

# 1. Asymptote du Module 
# Basse fréquence (f < f0): 0 dB
# Haute fréquence (f >= f0): Pente -40 dB/décade, passant par 0 dB à f0
asymptote_dB = np.where(f < f_0, 
                        0, # 0 dB pour f < f0
                        -40 * np.log10(f / f_0)) # -40 dB/décade pour f >= f0

# 2. Asymptote de la Phase (Transition Verticale / Échelon 0° à -180°)
# La phase saute de 0° à -180° à la fréquence de résonance f_0
phase_asymptote = np.where(f < f_0, 
                           0, # 0° pour f < f0
                           -180) # -180° pour f >= f0


# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(12, 10))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
plt.semilogx(f, module_H_dB, color='blue', label=r'Courbe réelle $|V_C / V_g|$')
plt.semilogx(f, asymptote_dB, color='red', linestyle='--', linewidth=1.5, label='Asymptotes (0 dB / -40 dB/décade)')

plt.title('Diagramme de Bode - Module du Filtre RLC Passe-Bas')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--", alpha=0.7)

# Lignes de référence
plt.axvline(x=f_0, color='red', linestyle=':', label=f'$f_0$ = {f_0:.2f} Hz')
plt.axhline(y=0, color='gray', linestyle=':', label='Gain asymptotique (0 dB)')
# Note: L'atténuation à -3 dB se produit à la fréquence de coupure f_c

plt.legend(loc='lower left')
# Ajustement des limites Y. La surtension réelle sera d'environ 20*log10(Q) = 14.3 dB.
plt.ylim([-50, 20]) 

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange', label=r'Phase réelle $\phi(V_C) - \phi(V_g)$')
plt.semilogx(f, phase_asymptote, color='red', linestyle='--', linewidth=1.5, label='Asymptotes (0° / Transition Verticale / -180°)')

plt.title('Diagramme de Bode - Phase du Filtre RLC Passe-Bas')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--", alpha=0.7)

# Lignes de référence de phase
plt.axvline(x=f_0, color='red', linestyle=':', label=f'$f_0$ = {f_0:.2f} Hz')
plt.axhline(y=0, color='gray', linestyle=':', label='$0^{\\circ}$ (Asymptote B.F.)')
plt.axhline(y=-180, color='gray', linestyle=':', label='$-180^{\\circ}$ (Asymptote H.F.)')
plt.axhline(y=-90, color='magenta', linestyle='-.', alpha=0.8, label=r'$-90^{\circ}$ (à $f_0$)')

plt.legend(loc='lower left')
plt.yticks(np.arange(-180, 1, 30))

plt.tight_layout()
plt.show()