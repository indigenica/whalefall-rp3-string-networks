# Auditoría independiente U4–U11 · Estimadores de red (R_Γ, R_ξ, drenaje, X, q)

*11 de octubre de 2026, noche (Madrid) · auditor independiente (Claude, sesión en la nube).*
- *Plan y criterios:* `audit/u4_u11/PLAN.md`, fijados en el commit `84fcf65` antes de escribir código.
- *Código:* `a_lib.py` (estimadores reimplementados desde los preregistros), `a0_units.py` y `a1_network.py`.
- *Salida:* `a1_network.json`. Se ejecuta en unos 15 s: `python3 audit/u4_u11/a1_network.py`.
- *Independencia:* no se ha leído ni importado `network/*_analyze.py`, `network/u9_extra.py`, `network/u10a_extra.py` ni `lya/*.py`. Solo los preregistros, los informes, los JSON de resultados y las redes de `network/results/u4/runs/`.

## Resumen en cuatro líneas

1. **R_Γ, R_ξ, el drenaje semilla a semilla, ρ_M/ρ_a y las pendientes dξ/dℓ se reproducen** en U6, U9, U10a, U10b y U11, casi siempre hasta la cuarta cifra. R_Γ no depende del esquema de derivada: tres esquemas razonables dan 1.562–1.565 en U11, con σ = 0.094.
2. **El intercambio X por semilla no coincide** con lo publicado (|Δf_X| ≤ 0.004, de una emisión de ~1–2 %). Es una cuestión de definición que el preregistro no fija. **No cambia nada:** R_Γ^corr sale 1.555 frente a 1.544 (0.13σ), con el mismo veredicto.
3. **q se reproduce bit a bit con el convenio del autor, pero ese convenio decide el veredicto de U11.** La ventana k_fís ∈ [10H, m/4] se puede evaluar en el instante final del par (lo que hace el autor), en su media o en el instante inicial. Con el primero, q + 2σ = 0.993 < 1 (Q_MENOR_QUE_1). **Con los otros dos, q + 2σ = 1.004 y 1.023: Q_COMPATIBLE_CON_1.** Según el criterio del plan, esto es un FAIL de robustez del veredicto q de U11.
4. **La σ de bootstrap de q infravalora mucho su incertidumbre** donde hay pocas semillas: en U10a (dos semillas) la σ publicada es 0.005–0.012, y el convenio mueve q entre 0.48 y 0.79. En U9 (densas), entre 0.64 y 1.20.

**La puerta δA de U7c se reproduce** en sus dos filas decisivas (§7, añadida después, cuando CLASS ya funcionaba en la sesión). κ(x₀) sale a 0.04 %, y δA un +2.6 % por encima en filas y referencias por igual, así que ninguna fila cambia de lado del umbral.

## 1. Convenios leídos de los datos

Comprobados directamente sobre las 1380 filas de las 46 redes, sin el código del autor (`a0_units.py`):
- H = 1/(aτ) y ℓ = ln(m_r/H), con m_r = √m2 = 0.490, exactos a precisión de máquina;
- ρ_a y ρ_M de las máscaras son densidades físicas; la suma de `spectrum.dE` es ρ_a/2 (la parte cinética) a 3×10⁻⁷ para ℓ ≥ 4;
- ε_lat = Γ/(8πH³f²ℓ · 0.24ℓ) con f² = 0.5 reproduce las medias por ℓ publicadas (`rp3_mean`, `u1_mean`) a 10⁻⁴.

## 2. Tabla de comprobaciones (resumen)

64 PASA, 1 MATIZADO y 7 DISCREPANCIA. Las 7 discrepancias son las 4 de f_X y 3 de q con el convenio propio del auditor. La tabla completa está en `a1_network.json`.

| Estudio | Comprobación | Auditado | Publicado | Estado |
|---|---|---|---|---|
| U11 | R_Γ drenada | 1.5617 | 1.5617 ± 0.094 | PASA |
| U11 | R_Γ, otros dos esquemas de derivada | 1.5645 / 1.5644 | 1.5617 | PASA |
| U11 | R_ξ drenada | 1.6072 | 1.6072 ± 0.058 | PASA |
| U11 | R_Γ^corr | 1.555 | 1.544 ± 0.087 | PASA |
| U11 | drenada/U6 por semilla (6 RP³ + 4 U(1)) | 0.748–0.864 / 0.666–0.762 | idénticos a 10⁻⁴ | PASA |
| U11 | ρ_M/ρ_a al final (4 casos) | 0.142 / 0.186 / 0.898 / 0.596 | 0.142 / 0.187 / 0.897 / 0.596 | PASA |
| U11 | f_X U(1), semillas 61–64 | −0.002 / −0.009 / −0.001 / −0.019 | −0.006 / −0.013 / −0.003 / −0.023 | DISCREPANCIA (sin efecto) |
| U11 | q_RP3(ℓ = 5), convenio del autor | 0.873 | 0.873 ± 0.059 | PASA |
| U11 | q_RP3(ℓ = 5), media geométrica del par | 0.892 | 0.873 | PASA (0.3σ), **cambia el veredicto** |
| U6 | R_Γ, R_Γ integral, R_ξ | 1.3893 / 1.3719 / 1.6166 | iguales | PASA |
| U6 | q_RP3, convenio del autor / media | 1.026 / 1.109 | 1.026 ± 0.118 | PASA / MATIZADO |
| U9 | R_Γ y R_ξ densas | 1.6884 / 1.5809 | iguales | PASA |
| U9 | dξ/dℓ RP³ densas y poco densas | 0.312, 0.308 / 0.313 | iguales | PASA |
| U10a | R_Γ drenadas, densas / poco densas | 1.7554 / 1.6308 | iguales | PASA |
| U10a | brecha G | 0.0736 | 0.0736 | PASA |
| U10a | q_RP3, convenio del autor | 0.700 / 0.593 | iguales | PASA |
| U10a | q_RP3, media geométrica | 0.648 / 0.479 | 0.700 ± 0.012 / 0.593 ± 0.005 | DISCREPANCIA |
| U10b | pendientes, ξ final, R_Γ y q de las lanzadas tarde | −0.200 / 1.730 / 1.302 / 0.186 | iguales | PASA |

## 3. q: reproducible, pero frágil

**Qué hace el estimador** (ESPEC_U2 §3, ESPEC_U6). Entre dos instantáneas con espectro, F(k) = a₂⁴ dE₂(k) − a₁⁴ dE₁(k) es la energía Goldstone nueva en cada bin de k comóvil. Se ajusta F ∝ k^(−q) en k_fís = k/a ∈ [10H, m/4] con al menos 8 bins. Después, q(ℓ) = q0 + q1(ℓ − 5) para ℓ ≥ 4.6.

**Lo que el preregistro no fija** es en qué instante del par se evalúan a y H para convertir esa ventana física en bins comóviles. El autor usa el instante final. Es una elección legítima, pero no la única.

**Cómo cambia q0(ℓ = 5)** según el convenio:

| | final del par (autor) | media geométrica | inicial del par | k en el borde del bin | ajuste ponderado por F |
|---|---|---|---|---|---|
| U11 (6 RP³) | **0.873 ± 0.060** | 0.892 ± 0.058 | 0.905 ± 0.059 | 0.909 ± 0.066 | 0.938 ± 0.039 |
| U6 (6 RP³) | **1.026 ± 0.116** | 1.109 ± 0.127 | 0.904 ± 0.046 | 0.944 ± 0.038 | 0.903 ± 0.041 |
| U9 densas (2 RP³) | **0.738 ± 0.067** | 0.806 ± 0.146 | 1.166 ± 0.233 | 1.196 ± 0.177 | 0.637 ± 0.119 |
| U10a densas (2) | **0.700 ± 0.012** | 0.648 ± 0.054 | 0.621 ± 0.045 | 0.786 ± 0.086 | 0.716 ± 0.065 |
| U10a poco densas (2) | **0.593 ± 0.005** | 0.479 ± 0.008 | 0.558 ± 0.054 | 0.599 ± 0.013 | 0.793 ± 0.023 |
| U10b lanzadas tarde (3) | **0.186 ± 0.041** | 0.172 ± 0.028 | 0.166 ± 0.029 | 0.177 ± 0.086 | 0.152 ± 0.042 |

Las tres primeras columnas son la ambigüedad más limpia, la del instante. Las dos últimas son variantes del ajuste, menos neutras: con bins lineales en k, usar el borde superior en lugar del centro deforma los primeros bins.

**Lectura.**
- **El espectro no es una ley de potencias limpia** en la ventana: si lo fuera, desplazar unos pocos bins no movería q. Que lo haga indica curvatura, o ruido en los primeros bins, que son pocos (9–19 en U11).
- **El veredicto preregistrado de U11, Q_MENOR_QUE_1, depende del convenio.** Con el del autor, q + 2σ = 0.993. Con la media geométrica, 1.004; con el inicial, 1.023. INFORME_U11 ya decía que el margen era de 0.01 y que q es el estimador más inestable. La auditoría lo cuantifica: **la sistemática del convenio (≈ 0.03 entre instantes, ≈ 0.065 incluyendo las variantes del ajuste) es del orden de la σ estadística (0.06)** y no está contabilizada.
- **En U9 y U10a la σ de bootstrap no sirve como error de q.** Con dos semillas, el bootstrap sobre semillas da 0.005–0.07, y el convenio mueve q entre 0.1 y 0.5. Las cifras «q ≈ 0.6–0.7» de U10a y «q ≈ 0.74» de U9 deberían leerse como «q entre ≈ 0.5 y 1.2».
- **R_Γ y R_ξ no tienen este problema**: son integrales del espectro, no pendientes.

**Recomendación** (no cambia ningún veredicto por sí misma; eso le toca al autor, con una enmienda):
1. dar q_RP3(ℓ = 5) en U11 como **0.87 ± 0.06 (estad.) ± 0.03–0.07 (convenio)** y decir que es compatible con 1 a 2σ si se suma la sistemática;
2. en el artículo y en el contrato, no presentar Q_MENOR_QUE_1 como un resultado robusto: q < 1 es lo que prefieren los datos con este estimador, no algo establecido;
3. la tarea 6 (U12) ya recorre q ∈ {0.6, 0.75, 0.87, 1.0}, así que cubre este rango.

## 4. El intercambio X

ESPEC_U10a define la emisión limpia como Γ − ⟨X⟩, «con ⟨X⟩ interpolado al mismo ℓ», pero no fija cómo se promedian las ~8 muestras de X de cada instantánea ni cómo se forma f_X. El auditor promedia por instantánea y toma el cociente de medias en la ventana. Sale |f_X| ≤ 2 % en las cuatro semillas, como en el informe, pero con valores distintos por semilla (diferencias de 0.002–0.004). El efecto en R_Γ^corr es de +0.011 (0.13σ). **La conclusión de U11 («X es casi nulo una vez drenado el baño») se mantiene.**

## 5. Detalles menores

- **La ventana de R_Γ pierde su último punto.** La derivada centrada necesita una instantánea posterior, y en ℓ_fin − 0.05 no la hay: el último punto de la rejilla sale NaN, también en los JSON publicados (`per_ell`). En la práctica la ventana es [ℓ_fin − 1.3, ℓ_fin − 0.19]. No cambia nada, pero conviene decirlo.
- **R_Γ publicado es el cociente de las medias en la ventana**, no la media de los cocientes por ℓ. Las dos dan lo mismo a 10⁻³ en U11.
- **La σ de bootstrap con dos semillas** (U9, U10a) es una cota inferior poco fiable, también para R_Γ (U10a poco densas: 0.019).

## 6. Pendiente

- **La puerta de forma de U7c** (el ajuste a la familia térmica de 149 modelos) no se ha reimplementado; solo δA (§7).
- **El kernel de la red y el drenaje en sí** (que el rozamiento actúe solo fuera de los núcleos y se apague en ℓ = 3.8) no se han auditado: necesitarían ejecutar redes. Una red de 128³ con y sin drenaje cabe en esta máquina.

## 7. La puerta δA de U7c sobre dos filas (addenda 1 del plan)

*`a3_lya_gate.py` → `a3_lya_gate.json`. CLASS v3.4.1 con el parche, compilado en la sesión, sin caché. Tarda unos 80 minutos con dos núcleos compartidos.*

- **Reimplementado:** κ(x₀) de U7a, con una integral analítica propia en era de radiación pura; el estimador δA, con cuadratura propia en ln k y la cola más allá de 1000 h/Mpc extrapolada como ley de potencias; las referencias y el umbral estricto.
- **No reimplementado:** la distribución en momento de la fila (S6) y los diccionarios de parámetros de CLASS.

| Comprobación | Auditado | Publicado | Estado |
|---|---|---|---|
| ν_S4, x₀ = 10 / 30 | 3.0841 / 1.0392 | 3.0853 / 1.0396 | PASA (0.04 %) |
| κ, x₀ = 10 / 30 | 0.5382 / 1.5973 | 0.5380 / 1.5967 | PASA |
| umbral estricto (mezcla de 1 keV) | 0.2404 | 0.2341 | PASA (+2.7 %) |
| B2κ, R_Γ = 2.09, x₀ = 30, K_sw = 0: f_em / δA | 0.4363 / 0.2596 | 0.4362 / 0.2530 | PASA (+2.6 %) |
| B2κ, R_Γ = 1.69, x₀ = 30, K_sw = 0: f_em / δA | 0.3849 / 0.2246 | 0.3848 / 0.2188 | PASA (+2.6 %) |

**Las dos filas quedan del mismo lado del umbral:** la de 2.09 fuera, la de 1.69 dentro. Todos los δA del auditor salen un 1.3–2.7 % más altos, filas y referencias por igual. Es la diferencia de cuadratura (la cola y la rejilla en k), y se cancela en el cociente fila/umbral. La κ analítica de la era de radiación pura coincide con la numérica del autor (que sigue g*(T)) a 0.04 %.
