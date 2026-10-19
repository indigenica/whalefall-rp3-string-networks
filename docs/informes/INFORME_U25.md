# U25 · El paso de tiempo de producción queda dentro de lo preregistrado (< 2 % en R_Γ, < 1 % en R_ξ)

*19 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U25.md` (commit `55070bd`).
  - Enmienda 1 (18 oct, ≈16:50): el drenaje es de primer orden en dτ, así que el sesgo es b = 2·|Δ|. Es anterior a cualquier red de U25.
  - Nota del 18 oct, 20:53: sin veredicto con una sola semilla.
- *Código:* `network/lattice.py` (`tau_grid`), `network/u25_produce.py` y `network/u25_analyze.py`.
- *Resultado:* `network/results/u25/U25_RESULT.json`; las redes, en `network/results/u25/runs/`.
- *Ejecución:* en la nube, 6 redes de 384³ con dτ = 0.125 (≈ 1.5–1.7 h cada RP³ y ≈ 14 min cada U(1)). Las parejas con dτ = 0.25 son las D0 de U24.

## 0. En una frase

**Veredicto preregistrado: PASO_ROBUSTO.**
- Con dτ/2, R_Γ cambia un −0.10 % ± 0.14 % y R_ξ un −0.16 % ± 0.09 %.
- Con el sesgo de primer orden de la enmienda 1, el paso de producción (dτ = 0.25) sesga R_Γ en b = 0.2 % y R_ξ en 0.3 %.
- **b + 2σ_b vale 0.8 % en R_Γ y 0.7 % en R_ξ,** frente a los umbrales del 2 % y del 1 %.

## 1. Las cantidades

Son 3 semillas. Δ_X = X(dτ = 0.125)/X(0.25) − 1 en la ventana [4.36, 5.21], con bootstrap emparejado.

| | Δ | σ (bootstrap) | σ entre semillas (informativo) | b = 2·\|Δ\| | b + 2σ_b |
|---|---|---|---|---|---|
| **R_Γ** | −0.10 % | 0.14 % | 0.18 % | 0.20 % | **0.78 %** (umbral 2 %) |
| **R_ξ** | −0.16 % | 0.09 % | 0.12 % | 0.31 % | **0.68 %** (umbral 1 %) |
| ε_RP³ | +0.04 % | 0.12 % | | 0.09 % | |
| ε_U(1) | +0.15 % | 0.05 % | | 0.29 % | |
| ξ_RP³ | +0.06 % | 0.07 % | | 0.13 % | |
| ξ_U(1) | +0.22 % | 0.03 % | | 0.44 % | |

- **También se cumple con el error entre semillas,** algo mayor: b + 2σ_b vale 0.93 % en R_Γ y 0.78 % en R_ξ.
- **El único cambio claramente distinto de cero** es ξ_U(1) (+0.2 %, unas 7σ). Es pequeño y se compensa en parte en el cociente.

## 2. Qué significa

- **El paso de tiempo de producción no limita el artículo.** Queda dentro de lo preregistrado: < 2 % en R_Γ y < 1 % en R_ξ. b + 2σ_b vale 0.8 % y 0.7 %, o 0.9 % y 0.8 % con el error entre semillas, muy por debajo de sus errores (±0.12 y ±0.06).
- **El 0.016 de la prueba de 128³ de U2** se hizo sin drenaje y con otro protocolo. U25 lo sustituye para el protocolo drenado de producción.
- **La semilla 61 sola** daba Δ_RΓ = −0.4 %, y las tres juntas, −0.1 %. La dispersión entre semillas es pequeña porque las parejas comparten las condiciones iniciales.

## 3. Lo que U25 no resuelve

- **No es 576³,** aunque con m_r·a·dx = 1 al final en los dos tamaños se espera el mismo efecto.
- **No prueba dτ = 0.0625,** por coste. Con un sesgo medido de 0.2 %, no parece necesario.
