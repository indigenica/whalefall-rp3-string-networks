# U4 · CLASS para Lyman-α y redes grandes para el atractor

*8 de octubre de 2026 · informe parcial: U4a completo, U4b fase 1 completa, fase 2 (RP³) en marcha.*
- *Preregistro:* `docs/preregistros/ESPEC_U4.md`, con la enmienda 1 (cuadratura en momentos), escrita antes de calcular las filas.
- *Resultados:* `lya/results/u4a/U4A_RESULT.json` (el log, en `u4a_class.log`) y `network/results/u4/U4B_RESULT.json`, con las redes en `runs/`.
- *Máquina:* WSL de indi-xeon (44 núcleos, 52 GB). CLASS v3.4.1 (commit 814ce59) con `env/class_qm_trapz_log.patch`.

## 0. En una frase

**Con CLASS, el sector oscuro supera Lyman-α en la misma ventana que con nuestro resolvedor: ramas B2 y B3 con espectro log-plano.** El veredicto preregistrado de U4a es **CIERRE_SEGÚN_RAMA**, con todas las validaciones superadas.

Las redes de 768³ llegan por primera vez al atractor (**U4b fase 1: GATES_PASSED**). Para que el cierre sea firme falta medir R_Γ en ese atractor (fase 2). Queda además abierta una pregunta nueva: ξ crece despacio con ℓ (dξ/dℓ ≈ 0.06–0.18), y hay que ver si ese crecimiento, extrapolado, es compatible con Lyα (U5).

## 1. U4a · Lyman-α con CLASS

**Validaciones** (todas superadas; ninguna ejecución de CLASS falló):

| | Resultado | Criterio |
|---|---|---|
| C0, lectura del fichero de distribución (Fermi–Dirac frente a analítica) | 0.37 % | < 0.5 % |
| C0b, cuadratura del parche (20 nodos) | δA a 0.006 %; T² a 0.04 % en 0.5–20 h/Mpc | < 0.5 % |
| C1, precisión (90 nodos, ℓ_max = 50, tolerancias ×0.1), dos filas decisivas | δA cambia un 0.88 % y un 0.14 % | < 2 % |
| C2 (informativa), frente a `lya/boltz.py` en esas filas | T² a 0.25 % y 0.63 % en 0.5–20 h/Mpc; δA de boltz un 9–12 % mayor | — |

C2 era la validación que le faltó a U3a. Nuestro resolvedor reproduce la forma de CLASS a menos del 1 % y es algo más severo en δA. **U3a queda confirmada.**

**Umbrales con CLASS.**
- El umbral estricto, mixta de 1 keV con f = 0.14, da δA = 0.234.
- El amplio, térmica de 5.3 keV, da 0.336.
- Referencias térmicas: 2 keV → 0.682; 3.5 keV → 0.475.

**Filas (62).** Primaria = puerta de forma estricta. Todas las filas son representables (rms < 0.1); en las permitidas, rms ≤ 0.008.

| Rama · espectro | δA | f_em | m_G (eV) | Permitidas |
|---|---|---|---|---|
| B2 log-plano (−1σ / central / +1σ) | 0.114–0.200 | 0.14–0.58 | 0.85–5.5×10⁻¹⁵ | **12/12** |
| B3 log-plano | 0.070–0.133 | 0.08–0.43 | 1.3–6.9×10⁻¹⁵ | **12/12** |
| B2 y B3, espectro IR | 0.61–0.85 | 0.88–0.99 | ~10⁻¹⁸–10⁻¹⁷ | 0/24 |
| B1 (los dos espectros) | 0.65–0.99 | 0.74–1.00 | ~10⁻²⁰–10⁻¹⁶ | 0/14 |

- La fila permitida más próxima al umbral es B2 +1σ, x₀ = 10, K_sw = 0: δA = 0.200, a un 15 % del umbral.
- La radiación oscura en las filas permitidas es ΔN_eff(CMB) ≤ 0.006.

**Veredicto preregistrado: CIERRE_SEGÚN_RAMA** (B2 y B3 permitidas a ±1σ, B1 excluida). Coincide con los umbrales estrictos y amplios y con los criterios de forma y de δA. La puerta de red viene de U4b fase 1, que pasó.

**Matiz sobre la puerta de forma.**
- En las 24 filas permitidas, el mejor ajuste cae en el borde inferior de la malla de mezclas térmicas: f = 0.05 con m ≈ 0.5–1 keV.
- La forma pide una fracción caliente todavía menor. Es lo esperable: la distribución tiene una cola dN/dln q ∝ 1/q, es decir, un pequeño porcentaje de partículas muy rápidas sobre un grueso frío.
- El límite de Garcia-Gallego et al. (f ≤ 0.14 m^1.1) se aplica aquí por debajo de 1 keV, fuera del rango en que se ajustó.
- El veredicto no depende de esto, porque la puerta δA pasa por sí sola en las 24 filas. Aun así, U5b amplía la malla.

## 2. U4b fase 1 · el control U(1) a 768³

Cuatro redes U(1) de 768³, con ℓ_fin = 5.95 y ℓ_ref = 5.65:

| Red | Horas | ξ(3.0) | ξ(4.0) | ξ(4.5) | ξ(5.0) | **ξ(5.65)** | ξ(5.95) |
|---|---|---|---|---|---|---|---|
| ξ_i = 0.25, semilla 61 | 5.2 | 0.331 | 0.409 | 0.453 | 0.507 | **0.568** | 0.613 |
| ξ_i = 2.0, semilla 61 | 5.8 | 1.621 | 0.519 | 0.538 | 0.553 | **0.613** | 0.610 |
| ξ_i = 0.5, semilla 61 | 6.9 | 0.432 | 0.448 | 0.478 | 0.521 | **0.565** | 0.629 |
| ξ_i = 0.5, semilla 62 | 5.4 | 0.394 | 0.408 | 0.467 | 0.546 | **0.668** | 0.752 |

- **G1′:** ξ_U1(ℓ_ref) = 0.62 (0.565 y 0.668), dentro de [0.5, 2.5]. **Pasa.**
- **G2′:** memoria de la condición inicial, |ξ(2.0) − ξ(0.25)|/media = **7.6 %**, por debajo del 20 %. **Pasa.**
  - La tendencia con N es 31 % a 256³ (U3b, ℓ = 4.5), 13 % a 384³ (prueba en la nube, ℓ = 4.96) y 7.6 % a 768³.
- **Estado: GATES_PASSED.** Es la primera vez que el control U(1) llega al atractor en este proyecto.
- **Reproducibilidad:** en ℓ = 4.5, ξ = 0.45–0.48, lo mismo que U2 y U3b a 256³ (0.44–0.45). El tamaño de la caja no cambia ξ a un ℓ dado.
- **G3′ (informativa):** dξ/dℓ en [4.95, 5.65] = **0.062, 0.181, 0.093 y 0.085** (media 0.105; dispersión entre semillas grande).
  - No es la pendiente 0 de Hindmarsh et al. ni la de 0.24–0.25 de GHV y Buschmann et al.
  - Nuestro ℓ_max = 6 está por debajo del estado del arte (ℓ ≈ 8–9). Este dato es un control, no una contribución a esa disputa.

## 3. Qué queda abierto

1. **U4b fase 2** (en marcha): dos redes RP³ y dos U(1) a 576³, R_Γ en [ℓ_ref − 1, ℓ_fin − 0.05] y la pendiente espectral.
   - Si |ΔR_Γ| ≤ 0.3 frente a 1.69, el cierre según rama es firme.
   - Si no, se repite U4a con el R_Γ nuevo (ESPEC_U4).
2. **La pendiente de ξ (U5a, por preregistrar antes de calcular).** B2 supone ξ constante (1.2 en ℓ ~ 70) y B1 supone ξ ≈ 0.24 ℓ (≈ 17).
   - Extrapolada linealmente, la pendiente medida da ξ(70) ≈ 7, entre las dos.
   - Hay que pasar esa rama por S6 y Lyα con CLASS.
   - Si la rama queda excluida, la supervivencia del sector oscuro pasa a depender de que ξ deje de crecer entre ℓ ≈ 6 y 70, que es justo la disputa abierta en la literatura.
3. **Malla de forma ampliada (U5b):** m < 0.5 keV y f < 0.05.
4. **Isocurvatura de ruido blanco (U5c).** En el escenario post-inflacionario, las cuerdas dejan fluctuaciones de isocurvatura con espectro de ruido blanco por debajo de k_* ≈ a H en H = m_G (Gorghetto, Trifinopoulos y Valogiannis 2025, arXiv:2511.04734). Hay que comprobar que, con m_G ~ 10⁻¹⁵ eV, quedan muy por debajo de las escalas de Lyα.
5. **Literatura nueva.** Yamada (arXiv:2607.21701, julio de 2026) simula redes de cuerdas ℤ_N con vacío PSU(N); para N = 2, PSU(2) = SO(3) ≅ RP³, nuestro caso. Encuentra escalado y una densidad que crece con N² − 1, pero no mide la emisión. Nuestro R_Γ le es complementario y hay que compararlos.
