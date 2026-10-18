# U6 · La red RP³ con más estadística: R_Γ, R_ξ y la pendiente espectral q(ℓ) · Preregistro

*9 de octubre de 2026, ≈09:00 (Madrid). Escrito antes de producir las redes de U6.*
- Se conocen U4 y U5: DEPENDE_DEL_ESCALADO y DEPENDE_DE_q.
- El análisis `network/u6_analyze.py` se fijó antes de los datos y se probó sobre las dos semillas RP³ existentes (`network/results/u4/U6_dryrun_2seeds.json`). Esa prueba no cuenta como resultado.

## Motivo

U5 deja el sector oscuro de la ballena colgado de dos números en el ℓ físico (≈ 70):
- **el crecimiento de ξ:** con la extrapolación lineal de nuestra pendiente, toda la rama B4 queda excluida; la pendiente crítica es < 0.062 y la medida, 0.133 ± 0.030;
- **la pendiente espectral q:** q_crit = 1.025 para B2.

Ninguna red alcanza ℓ ≈ 70. U6 no puede decidir eso, y no da veredicto sobre la ballena. Su objetivo es dejar bien medido lo que sí es nuestro y nuevo: la emisión de una red RP³ (≅ SO(3)) frente a una U(1) y la evolución de su espectro en la ventana accesible. Ese resultado se puede publicar por sí mismo y comparar con Yamada (arXiv:2607.21701).

## Producción

`u4b_produce.py 4`, a 576³:
- RP³ con ξ_i = 0.5, semillas 63–66;
- U(1) con ξ_i = 0.5, semillas 63 y 64;
- espectros en cada una de las 30 instantáneas, no cada tres.

Mismo kernel, máscaras y paso que U4b.

## Análisis (`network/u6_analyze.py`)

Usa las semillas 61–66 de RP³ y 61–64 de U(1). Las 61 y 62 tienen espectros cada tres instantáneas.
- **R_Γ y R_ξ** en la ventana de U4b, [ℓ_fin − 1.3, ℓ_fin − 0.05], con los dos estimadores de U2/U4b y bootstrap sobre semillas.
- **q(ℓ):** pendiente de la emisión instantánea F(k) ∝ k^(−q) entre dos instantáneas consecutivas, en k_fís ∈ [10 H, m/4], con al menos 8 bins de k.
  - Ajuste q(ℓ) = q0 + q1 (ℓ − 5) para ℓ ≥ 4.6, con las semillas juntas y bootstrap sobre semillas.
  - Se informa q en ℓ_fin y la tendencia: CRECE si q1 > 2σ, DECRECE si q1 < −2σ, y SIN_TENDENCIA_SIGNIFICATIVA en otro caso. Lo mismo para U(1), como comparación con la literatura.
- **Pendientes dξ/dℓ** de RP³ y U(1) en [ℓ_ref − 0.7, ℓ_ref].

**Lectura (informativa):**
- Si q_RP3 CRECE y ya supera 1 al final de la ventana, la vía de escape «q ≤ 1 en el ℓ físico» exigiría que la tendencia se invirtiera, en contra de lo medido.
- Si q_RP3 no tiene tendencia y es ≤ 1, la vía queda abierta.
