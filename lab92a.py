import numpy as np
import matplotlib.pyplot as plt

# --- 1. Paramètres du Circuit RLC (inchangés) ---
L = 27e-3  
C = 0.1e-6 
R = 100    

# --- 2. Analyse Théorique ---
omega_0 = 1 / np.sqrt(L * C)
f_0 = omega_0 / (2 * np.pi) 
Q = (1/R) * np.sqrt(L/C) # Q ≈ 5.2

# Calcul du décalage de l'asymptote à la résonance
gain_asymptote_fo = 20 * np.log10(1 / Q)
print(f"**Gain asymptotique à f0 (G_croisement):** {gain_asymptote_fo:.2f} dB")

# Plage de fréquences
f_min = f_0 / 100 
f_max = f_0 * 100 
f = np.logspace(np.log10(f_min), np.log10(f_max), 500)
omega = 2 * np.pi * f

# Fonction de transfert H(jω)
Z_total = R + 1j * omega * L + (1 / (1j * omega * C))
H = R / Z_total
module_H_dB = 20 * np.log10(np.abs(H))
phase_H = np.angle(H, deg=True) 

# --- Asymptotes Corrigées ---

# 1. Asymptote du Module (Passe par G_croisement à f0)
asymptote_dB = np.where(f < f_0, 
                        20 * np.log10(f / f_0) + gain_asymptote_fo, # Pente +20 dB/décade décalée
                        -20 * np.log10(f / f_0) + gain_asymptote_fo) # Pente -20 dB/décade décalée

# 2. Asymptote de la Phase (Échelon +90° à -90°) - Inchagée car correcte
phase_asymptote = np.where(f < f_0, 90, -90) 


# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(12, 10))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
plt.semilogx(f, module_H_dB, color='blue', label=r'Courbe réelle $|H(f)|$')
plt.semilogx(f, asymptote_dB, color='blue', linestyle='--', linewidth=1.5, label=f'Asymptote (Croisement à {gain_asymptote_fo:.2f} dB)')

plt.title('Diagramme de Bode - Module du Filtre RLC Passe-Bande')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--", alpha=0.7)

# Lignes de référence
plt.axvline(x=f_0, color='red', linestyle='--', label=f'$f_0$ = {f_0:.2f} Hz')
plt.axhline(y=0, color='gray', linestyle=':', label='Gain max Réel (0 dB)')
plt.axhline(y=gain_asymptote_fo, color='purple', linestyle=':', label=f'Point de croisement Asymptote ({gain_asymptote_fo:.2f} dB)')

plt.legend(loc='lower left')
plt.ylim([np.min(module_H_dB)-10, 20]) 

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange', label=r'Phase réelle $\phi(f)$')
plt.semilogx(f, phase_asymptote, color='orange', linestyle='--', linewidth=1.5, label='Asymptote de phase (Échelon +90°/-90°)')

plt.title('Diagramme de Bode - Phase du Filtre RLC Passe-Bande')
plt.xlabel(r'Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--", alpha=0.7)

plt.axvline(x=f_0, color='red', linestyle='--', label=f'$f_0$ = {f_0:.2f} Hz')
plt.axhline(y=0, color='blue', linestyle=':', label='$0^{\\circ}$ (à $f_0$)')
#plt.axhline(y=-90, color='gray', linestyle=':', label='$-90^{\circ}$')
#plt.axhline(y=90, color='gray', linestyle=':', label='$+90^{\circ}$')

plt.legend(loc='lower left')
plt.yticks(np.arange(-90, 91, 30))

plt.tight_layout()
plt.show()