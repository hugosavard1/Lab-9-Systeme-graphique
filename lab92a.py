import numpy as np
import matplotlib.pyplot as plt

# --- 1. Définition des Paramètres du Circuit RLC Série ---
L = 27e-3  # Inductance en Henry (27 mH)
C = 0.1e-6 # Capacité en Farads (0.1 μF)
R = 100    # Résistance en Ohms (Ajoutée pour l'amortissement et la forme du filtre)
V_g_peak_to_peak = 0.2 # Tension crête à crête de la source (V)

print(f"**Paramètres du circuit RLC :**")
print(f"R = {R} Ω")
print(f"L = {L*1e3:.0f} mH")
print(f"C = {C*1e6:.1f} μF\n")

# --- 2. Analyse Théorique ---

# Pulsation de résonance (naturelle)
omega_0 = 1 / np.sqrt(L * C)
f_0 = omega_0 / (2 * np.pi) # Fréquence de résonance en Hz

# Facteur de qualité Q (pour un RLC série, avec la sortie aux bornes de R)
Q = (1/R) * np.sqrt(L/C)

# Largeur de bande à -3dB (bande passante)
delta_omega = R / L
delta_f = delta_omega / (2 * np.pi)

# Fréquences de coupure théoriques (approximatives pour Q élevé)
f_c1 = f_0 * np.sqrt(1 + (1/(2*Q))**2) - f_0 / (2*Q) # Fréquence de coupure inférieure
f_c2 = f_0 * np.sqrt(1 + (1/(2*Q))**2) + f_0 / (2*Q) # Fréquence de coupure supérieure

# Pour des Q élevés, on peut aussi utiliser f_c1 = f_0 - delta_f/2 et f_c2 = f_0 + delta_f/2

print(f"**Fréquence de résonance (f_0):** {f_0:.2f} Hz")
print(f"**Facteur de qualité (Q):** {Q:.2f}")
print(f"**Largeur de bande (-3dB) (Δf):** {delta_f:.2f} Hz")
print(f"**Fréquence de coupure inférieure (f_c1):** {f_c1:.2f} Hz")
print(f"**Fréquence de coupure supérieure (f_c2):** {f_c2:.2f} Hz\n")

# Définir la plage de fréquences (centrée autour de la résonance)
f_min = f_0 / 10 # 1/10 de la fréquence de résonance
f_max = f_0 * 10  # 10 fois la fréquence de résonance
f = np.logspace(np.log10(f_min), np.log10(f_max), 500)
omega = 2 * np.pi * f

# Fonction de transfert pour un filtre RLC série avec sortie aux bornes de R (Passe-bande)
# H(jω) = V_R / V_g = R / (R + jωL + 1/(jωC))
# On peut simplifier en divisant par R: H(jω) = 1 / (1 + jω(L/R) + 1/(jωRC))
# Ou mieux, en termes de Q et ω0: H(jω) = 1 / (1 + jQ(ω/ω0 - ω0/ω))

Z_total = R + 1j * omega * L + (1 / (1j * omega * C))
H = R / Z_total

# Calcul du module et de la phase
module_H = np.abs(H)
phase_H = np.angle(H, deg=True) # Phase en degrés

# Calcul du module en décibel (dB)
module_H_dB = 20 * np.log10(module_H)

# --- 3. Tracé de la Réponse en Fréquence (Diagrammes de Bode) ---

plt.figure(figsize=(10, 8))

# --- Tracé du Module (dB) vs Fréquence (Log) ---
plt.subplot(2, 1, 1)
plt.semilogx(f, module_H_dB, color='blue', label='Courbe réelle $|H(f)|$')
plt.title('Diagramme de Bode - Module $\\left(|V_R / V_g|\\right)$ du Filtre RLC Passe-Bande')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Module (dB)')
plt.grid(which="both", ls="--")

# Ajout de la fréquence de résonance
plt.axvline(x=f_0, color='red', linestyle='--', label=f'Fréquence de résonance $f_0$ = {f_0:.2f} Hz')
# Ajout des fréquences de coupure
plt.axvline(x=f_c1, color='green', linestyle=':', label=f'$f_{{c1}}$ = {f_c1:.2f} Hz (-3 dB)')
plt.axvline(x=f_c2, color='green', linestyle=':', label=f'$f_{{c2}}$ = {f_c2:.2f} Hz (-3 dB)')
plt.axhline(y=-3, color='purple', linestyle=':', label='-3 dB') # Ligne des -3 dB
plt.legend()
plt.ylim([-40, 5]) # Ajustement de l'échelle pour une meilleure vue

# --- Tracé de la Phase (Degrés) vs Fréquence (Log) ---
plt.subplot(2, 1, 2)
plt.semilogx(f, phase_H, color='orange')
plt.title('Diagramme de Bode - Phase $\\left(\\phi_{V_R} - \\phi_{V_g}\\right)$ du Filtre RLC Passe-Bande')
plt.xlabel('Fréquence $f$ (Hz)')
plt.ylabel('Déphasage (Degrés)')
plt.grid(which="both", ls="--")

# Ajout de la fréquence de résonance (où la phase est de 0 degrés)
plt.axvline(x=f_0, color='red', linestyle='--', label=f'Fréquence de résonance $f_0$ = {f_0:.2f} Hz')
plt.axhline(y=0, color='blue', linestyle=':', label='$0^{\\circ}$') # Ligne des 0 degrés à la résonance
# Ajout des fréquences de coupure (où la phase est de +/- 45 degrés)
plt.axvline(x=f_c1, color='green', linestyle=':', label=f'$f_{{c1}}$ = {f_c1:.2f} Hz (-45°)')
plt.axvline(x=f_c2, color='green', linestyle=':', label=f'$f_{{c2}}$ = {f_c2:.2f} Hz (+45°)')
plt.axhline(y=-45, color='gray', linestyle=':')
plt.axhline(y=45, color='gray', linestyle=':')
plt.legend()

plt.tight_layout()
plt.show()