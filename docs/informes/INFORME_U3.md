# U3 · Las dos pruebas que deciden el sector oscuro

*6 de octubre de 2026*
- *Preregistro:* `docs/preregistros/ESPEC_U3.md`, con cuatro enmiendas fechadas. La cuarta es posterior a la auditoría.
- *Resultados:* `lya/results/U3A_RESULT.json` y `network/results/u3/U3B_RESULT.json`.
- *Figuras:* `lya/results/u3a_lya.png` y `network/figures/u3b_control.png`.
- *Auditoría independiente:* `audit/u3/AUDIT_U3.md`.

## 0. En una frase

**La red no se puede llevar al atractor con N = 256 (U3b: GATES_FAILED), así que el veredicto preregistrado sigue siendo INCONCLUSO. Pero la prueba de Lyman-α con una transferencia de Boltzmann cambia el diagnóstico del sector oscuro:**
- con un espectro de emisión log-plano, las ramas de escalado constante (B2) y de la tendencia medida (B3) pasan Lyα, incluso con R_Γ a +1σ;
- el crecimiento logarítmico (B1) y cualquier espectro dominado por el infrarrojo quedan excluidos.

Lectura condicional a la red: **CIERRE_SEGÚN_RAMA**, con el umbral estricto y con el amplio. El resolvedor no queda validado según su preregistro (W1).

## 1. U3b · ¿Se puede llevar la red al atractor? No, con N = 256

**Montaje.** Antes de pasar a cuerdas físicas, la red evoluciona con un núcleo de anchura comóvil fija, de 3 celdas («crecimiento de núcleo», como en Hindmarsh et al. 2021). Así se le da más tiempo conforme para relajarse antes de la ventana de análisis.

**Resultado** (control U(1), tres redes de 256³):

| | ξ_U1(ℓ = 4.5) |
|---|---|
| F1, ξ_i = 0.5 (semillas 41 y 42) | 0.52 y 0.35 (media 0.44) |
| F2, ξ_i = 2.0 (semilla 43) | 0.60 |
| U2 sin preparación (referencia) | 0.45 |

- **G1** (ξ_U1 ≥ 0.5): **falla**.
- **G2** (memoria de la condición inicial < 20 %): **falla**, con un 31 %.
- Veredicto: **GATES_FAILED**, confirmado por la auditoría. R_Γ sigue siendo el de U2, 1.69 ± 0.29.

**Por qué.** Con cualquier preparación, la red pasa por un colapso cuando el cociente horizonte/núcleo es ≈ e^2.5–e^3 (`network/figures/u3b_control.png`). Después sube despacio. Con N = 256 la ventana (ℓ ≤ 4.85) no basta para llegar al atractor. **Hace falta más rango dinámico: N ≳ 1024, es decir ℓ ≳ 6.2.**

Las redes RP³ de F se detuvieron (enmienda 3) cuando las dos puertas ya estaban decididas por el control U(1).

## 2. U3a · Lyman-α con una transferencia de Boltzmann

**Herramienta.**
- CLASS no se pudo instalar (PyPI y GitHub bloqueados).
- `lya/boltz.py` es un resolvedor lineal propio en gauge síncrono (Ma & Bertschinger 1995), que trata la distribución de momentos completa de las partículas emitidas.
- La puerta es el criterio δA de Murgia et al. 2017.
- La auditoría lo comprobó ecuación a ecuación y lo reprodujo con código independiente (≤ 5×10⁻⁴ en la transferencia; δA 0.13781 frente a 0.13783).

**Validación:**

| Prueba | Resultado | Estado |
|---|---|---|
| W1, ajuste de Viel et al. 2005 | k_{1/2} un 11–18 % mayor; δA un 4–10 % menor | **FALLA** (preregistrada) |
| W1′, ajuste δA(m) de Dienes et al. 2021 | retirada (enmienda 4): se construye con el ajuste de Viel, no es independiente | — |
| W2, límite frío | 2×10⁻¹⁰ | PASA |
| W3, convergencia en las filas decisivas: nodos, ℓ_max, paso, k_max = 1500 | ≤ 1.1 % | PASA (tolerancia del 3 %, enmienda 4) |

**Conclusión:** U3a no está validada según su preregistro.

**Por qué el resultado sigue siendo informativo.** Los umbrales se calculan con el mismo código que las filas. Según la auditoría, las simplificaciones del resolvedor (bariones, fotones sin acoplamiento fuerte, RSA) mueven el cociente fila/umbral ≤ 0.3 %.

**Umbrales** (mismo código, k_max = 1000 h/Mpc):
- **estricto:** δA ≤ 0.241, el de la mezcla en el borde de Garcia-Gallego et al. 2025 a 1 keV. Con su tabla, más permisiva, sería mayor.
- **amplio:** δA ≤ 0.351, el del WDM térmico de 5.3 keV.

**Resultado sobre las filas de S6** (K_sw = 0; ε medida en U2; tres ramas; R a −1σ / central / +1σ):

| Rama | ir_delta | energy_logflat |
|---|---|---|
| B1 (GHV, ξ ∝ ℓ) | excluida (δA ≈ 0.99) | excluida (δA ≈ 0.67–0.76) |
| B2 (Hindmarsh, ξ\* = 1.2) | excluida (δA ≈ 0.76–0.84) | **permitida (δA = 0.15–0.22)** |
| B3 (tendencia medida) | excluida (δA ≈ 0.62–0.72) | **permitida (δA = 0.09–0.15)** |

Las 12 filas permitidas con K_sw = 0 dan:
- m_G = 1.7–6.9×10⁻¹⁵ eV;
- una fracción emitida (tibia) f_em = 0.15–0.58;
- δA = 0.09–0.22;
- ΔN_eff = 0.002–0.006.

La fila más ajustada (B2 a +1σ, x0 = 10) queda a un 90 % del umbral estricto.

**Veredicto.**
- **Preregistrado: INCONCLUSO.** ESPEC_U2 §6 lo fuerza cuando falla el control U(1) de la red, y falla en U2 y en U3b.
- **Condicional a la red:** CIERRE_SEGÚN_RAMA, con el umbral estricto y con el amplio. B2 y B3 quedan permitidas a −1σ, central y +1σ, y B1 queda excluida en todos los casos.

**Lectura.**
- La transferencia de Boltzmann es menos severa que el indicador λ_fs de U1, porque la materia oscura es una mezcla: parte fría (desalineamiento) y parte tibia (emisión).
- Si el espectro de emisión es **log-plano** y la red escala con ξ constante, o sigue la tendencia medida, el sector oscuro sobrevive. Entonces deja una supresión de pequeña escala δA ≈ 0.1–0.2, medible con datos de Lyα de alta resolución.
- Si el espectro está dominado por el infrarrojo (la imagen de GHN a ℓ grande), el sector oscuro queda falsado en las tres ramas.
- La pendiente medida en la red RP³ es q = 0.76–0.97, cercana a log-plana. Pero está medida a ℓ ≈ 5 y no se puede extrapolar con fiabilidad.

**Advertencias.**
- δA es poco fiable para distribuciones multimodales (Dienes et al. 2021), y estas lo son: frío + tibio. Además depende de la convención de k_max (enmienda 4).
- No es CLASS.
- Falta confirmarlo con `classy` y, mejor aún, con simulaciones de Lyα para esta distribución concreta.

## 3. Estado tras U3

| Pieza | Antes de U3 | Después de U3 |
|---|---|---|
| Red, ε | R_Γ = 1.69 ± 0.29, INCONCLUSO | igual: la preparación de núcleo no basta con N = 256 |
| Lyα | indicador λ_fs: las tres ramas excluidas | Boltzmann + δA: **B2 y B3 sobreviven con el espectro log-plano**; B1 excluida (sin validación preregistrada) |
| Qué decide | — | (1) cómo escalan las redes a ℓ ≈ 70; (2) el espectro de emisión (q ≤ 1 frente a q > 1); (3) Lyα para la mezcla, con CLASS y simulaciones |

## 4. Reproducir

```
make network-u3     # U3b (~2–3 h)
make lya            # U3a (~1 h con k_max = 1000)
```
