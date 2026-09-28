"""
Résolution d'un circuit électrique à 10 résistances.

Les lois de Kirchhoff (lois des nœuds et des mailles) donnent un système
linéaire R . I = U de 10 équations. On le résout par élimination de Gauss
avec pivot partiel, puis rétro-substitution, codées sans bibliothèque externe.

Projet L3 Physique, CY Cergy Paris Université (septembre 2024).

Utilisation :
    python circuit.py
"""

import sys


def valeur_absolue(x):
    return x if x >= 0 else -x


# Seuil en dessous duquel un pivot est considéré comme nul
precision = 1e-6

# Valeurs des résistances (en ohms)
R1, R2, R3, R4, R5, R6, R7, R8, R9, R10 = 2400, 8000, 5300, 7900, 1900, 8900, 1400, 6000, 9700, 6800

# Tension appliquée (en volts)
V = 50

# Système R . I = U issu des lois de Kirchhoff :
# - 5 premières lignes : lois des nœuds (coefficients +1 / -1)
# - 5 dernières lignes : lois des mailles (coefficients = résistances)
R = [[-1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
     [0, -1, 0, 0, 1, 0, 1, 0, 0, 0],
     [0, 0, 1, 0, 1, -1, 0, -1, 0, 0],
     [0, 0, 0, 1, 0, 1, 0, 0, -1, 0],
     [0, 0, 0, 0, 0, 0, 1, 1, 0, -1],
     [R1, -R2, R3, 0, -R5, 0, 0, 0, 0, 0],
     [0, 0, -R3, R4, 0, -R6, 0, 0, 0, 0],
     [0, 0, 0, 0, R5, 0, -R7, R8, 0, 0],
     [0, 0, 0, 0, 0, R6, 0, -R8, R9, -R10],
     [R1, 0, 0, R4, 0, 0, 0, 0, R9, 0]]

U = [0, 0, 0, 0, 0, 0, 0, 0, 0, V]

# Copies du système d'origine, pour vérifier la solution à la fin
R_origine = [ligne[:] for ligne in R]
U_origine = U[:]

# ---------------------------------------------------------------------------
# Élimination de Gauss avec pivot partiel
# ---------------------------------------------------------------------------

n = len(U)
for k in range(n):
    # Recherche de la ligne avec le plus grand pivot (en valeur absolue)
    Max = valeur_absolue(R[k][k])
    p = k
    for i in range(k + 1, n):
        val = valeur_absolue(R[i][k])
        if val > Max:
            Max = val
            p = i

    if Max < precision:
        print("La matrice est singulière ou proche de l'être.")
        sys.exit(1)

    # Échange des lignes k et p
    R[k], R[p] = R[p], R[k]
    U[k], U[p] = U[p], U[k]

    # Élimination sous le pivot
    for i in range(k + 1, n):
        h = R[i][k] / R[k][k]
        for j in range(k, n):
            R[i][j] -= h * R[k][j]
        U[i] -= h * U[k]

# ---------------------------------------------------------------------------
# Rétro-substitution
# ---------------------------------------------------------------------------

I = [0.0] * n
I[-1] = U[-1] / R[-1][-1]
for i in range(n - 2, -1, -1):
    s = 0
    for j in range(i + 1, n):
        s += R[i][j] * I[j]
    I[i] = (U[i] - s) / R[i][i]

# ---------------------------------------------------------------------------
# Résultats
# ---------------------------------------------------------------------------

print("\nIntensités dans chaque résistance (une valeur négative signifie que")
print("le courant circule dans le sens opposé au sens choisi sur le schéma) :")
for i in range(n):
    print(f"  I{i + 1:<2} = {I[i]:+.3e} A")

# Courant fourni par le générateur
z = I[0] + I[1]
print(f"\nCourant du générateur : I = I1 + I2 = {z:.3e} A")

# Résistance équivalente vue par le générateur
R_total = V / z
print(f"Résistance équivalente : {R_total:.1f} ohms")

# Vérification : résidu max |R . I - U| sur le système d'origine
residu = max(valeur_absolue(sum(R_origine[i][j] * I[j] for j in range(n)) - U_origine[i])
             for i in range(n))
print(f"Vérification : résidu max |R.I - U| = {residu:.1e}")

# Vérification physique : bilan de puissance (conservation de l'énergie).
# La puissance fournie par le générateur doit être entièrement dissipée
# par effet Joule dans les résistances.
resistances = [R1, R2, R3, R4, R5, R6, R7, R8, R9, R10]
P_joule = sum(r * i ** 2 for r, i in zip(resistances, I))
print(f"Bilan de puissance : V.I = {V * z * 1e3:.3f} mW, somme des R.I² = {P_joule * 1e3:.3f} mW")
