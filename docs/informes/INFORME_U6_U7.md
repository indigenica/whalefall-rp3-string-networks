# U6 y U7 · La red RP³ con más estadística, y la emisión calibrada con las redes más grandes publicadas

*9 de octubre de 2026.*
- *Preregistros:* `docs/preregistros/ESPEC_U6.md` y `ESPEC_U7.md`.
- *Resultados:* `network/results/u4/U6_RESULT.json` (con las redes en `runs/`) y `lya/results/u7/U7_RESULT.json` (el log, en `u7_class.log`).

## 0. En una frase

**Con la emisión calibrada con Correia, Hindmarsh et al. (arXiv:2512.13653), el sector oscuro de la ballena pasa Lyα en las 24 filas: CALIBRADA_PERMITIDA.** La red RP³, medida ahora con seis semillas, da R_Γ = 1.39 ± 0.06 y una pendiente espectral q = 1.03 ± 0.12, compatible con el q = 1 de ese trabajo. La duda de U5 sobre q > 1 venía de una estimación pobre.

## 1. U6 · RP³ a 576³: seis semillas RP³ y cuatro U(1), espectros en cada instantánea

| | U4b (2 + 2 semillas) | **U6 (6 + 4 semillas)** |
|---|---|---|
| R_Γ (ventana [4.36, 5.61]) | 1.40 ± 0.13 | **1.389 ± 0.062** (bootstrap); integral 1.372 ± 0.077 |
| R_ξ | 1.65 ± 0.10 | **1.617 ± 0.053** |
| q_RP3 en ℓ = 5 | 1.10–1.22 (ℓ ≈ 5.4, dos semillas, espectros cada tres instantáneas) | **1.03 ± 0.12** |
| dq/dℓ (RP³) | ≈ 0.7 ± 0.2 (prueba con dos semillas) | **−0.12 ± 0.29**: SIN_TENDENCIA_SIGNIFICATIVA; q(ℓ_fin) = 0.95 |

- **R_Γ se confirma** y su error se reduce a la mitad.
- **La pendiente espectral de RP³ es compatible con 1 y no crece.** Es lo que concluyen Correia et al. para U(1) (q = 1).
  - La lectura de U5 (q ≈ 1.1–1.2) salía de dos semillas con pocos espectros.
  - Aun así, el error (±0.12) todavía no permite afirmar q ≤ 1.025, que es el q_crit de U5f.
- **El estimador de q para U(1) no es fiable** (q = −0.18 ± 0.26, con dispersión grande): ese canal tiene muy pocos bins y está dominado por el ruido. No se usa.
- **dξ/dℓ en [ℓ_ref − 0.7, ℓ_ref]:**
  - U(1): 0.094–0.115, lo mismo que a 768³;
  - **RP³: 0.26–0.37**, con ξ_RP3 ≈ 0.9–1.0.
  - La red RP³ crece más deprisa que la U(1), y R_ξ sube a lo largo de la ventana. Aun así, R_Γ, el cociente de emisión que entra en la cadena, no muestra tendencia.
  - Si es el transitorio «desde abajo» que describen Correia et al. (nuestras redes empiezan con ξ_i = 0.5) se puede comprobar barato: redes RP³ con densidad inicial alta (ξ_i = 2). Queda propuesto como U9.

## 2. U7 · Calibración de la emisión (arXiv:2512.13653)

**U7a.** Nuestra emisión de U(1) (B2, log-plano) escala: ν es constante en ℓ ≈ 100–113.

| x₀ | ν_S4 | κ = 1.66/ν_S4 |
|---|---|---|
| 10 | 3.08 | 0.54 ± 0.06 |
| 30 | 1.04 | 1.60 ± 0.16 |

El valor medido (1.66 ± 0.17) cae entre los dos cortes; equivale a x₀ ≈ 18.

**U7b** (CLASS, mismas puertas que U4a y U5; R_Γ = 1.40 ± 0.13):

| Familia | x₀ | f_em | δA (estricto 0.234) | m_G (eV) | Permitidas |
|---|---|---|---|---|---|
| B2κ (emisión calibrada con ν_ax) | 10 | 0.19–0.36 | 0.079–0.103 | 1.5–4.5×10⁻¹⁵ | 6/6 |
| B2κ | 30 | 0.19–0.36 | 0.155–0.203 | 1.5–4.5×10⁻¹⁵ | 6/6 |
| B2ζ (ξ* = 1.49) | 10 | 0.35–0.56 | 0.157–0.194 | 0.9–2.2×10⁻¹⁵ | 6/6 |
| B2ζ | 30 | 0.15–0.31 | 0.126–0.168 | 1.6–5.2×10⁻¹⁵ | 6/6 |

- **Veredicto U7: CALIBRADA_PERMITIDA.** Pasan las 24 filas, incluso con R a +1σ.
- La fila más justa es B2κ, x₀ = 30, R a +1σ, con δA = 0.203, un 13 % por debajo del umbral.
- La forma corresponde a un 1–3 % de una componente térmica de 0.2–0.26 keV, igual que en U5b.

## 3. Lectura de conjunto (U4–U7)

- **La predicción de la ballena 5 queda así:** una materia oscura de pseudo-Goldstones con m_G ≈ 1–5×10⁻¹⁵ eV, con un 20–35 % emitido por la red de cuerdas RP³. En Lyα eso se ve como un 1–3 % de componente muy caliente (~0.25 keV) y una supresión δA ≈ 0.08–0.20.
- **Su condición de vida** sigue siendo la de U5: densidad de cuerdas constante y q ≤ 1 en ℓ ≈ 70.
  - Es lo que concluye la simulación más grande publicada (12288³).
  - Nuestras redes RP³ dan ahora q = 1.03 ± 0.12, sin tendencia.
  - La única tensión que queda en nuestros datos es el crecimiento de ξ_RP3 en la ventana (U9).
- **Pendiente:**
  - la auditoría independiente de U4–U7;
  - U9, RP³ con densidad inicial alta;
  - las cotas de reliquias calientes para el ~2 % a 0.25 keV;
  - y un escrito sobre la red RP³ (R_Γ, R_ξ, q) frente a U(1), comparado con Yamada (arXiv:2607.21701).
