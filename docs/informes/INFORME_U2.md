# U2 · Cuánto emite de verdad la red RP³, y qué implica para el cierre

*6 de octubre de 2026 · network/*
- *Preregistro:* `ESPEC_U2.md`, con 5 enmiendas escritas antes de la producción.
- *Resultados:* `network/results/u2/U2_RESULT.json` y `results/s6_measured_closure.json`. *Figura:* `network/figures/u2_red.png`.
- *Auditoría independiente:* `audit/u2/AUDIT_U2.md`.

## 0. En una frase

**La red de cuerdas ℤ₂ de RP³ emite Goldstones libres 1.7 veces más que una red U(1) con la misma f, la misma masa de núcleo y en el mismo instante: R_Γ = 1.69 ± 0.29.** U1 pedía lo contrario, una red RP³ entre 20 y 1 000 veces menos eficiente que la U(1).
- **Veredicto preregistrado: INCONCLUSO.** El control U(1) queda justo por debajo del rango de la literatura (ξ_U1 = 0.45 frente a ≥ 0.5): la ventana de análisis aún no está en el régimen de escalado.
- **Con las puertas de U1, ninguna de las tres extrapolaciones de la literatura permite cerrar,** ni siquiera a −1σ.
- **Con un tratamiento de Lyman-α más fino (exploratorio),** que admite una fracción de materia oscura tibia, la ballena sobrevive si las redes escalan con ξ constante (Hindmarsh et al.) y el espectro es log-plano. Con el crecimiento logarítmico de GHV queda falsada.

## 1. Qué se ha hecho

**Montaje.**
- **Teoría:** la del campo Q de la acción común. Matriz simétrica 4×4 (10 componentes), vacío RP³ y cuerdas ℤ₂ con núcleo físico, en era de radiación.
- **Control:** una U(1) con la misma f² = ½ y la misma masa radial, m² = 0.24, en el mismo código C (leapfrog conforme, OpenMP).
- **Estimador:** el de GHV: radiación libre = 2 × energía cinética tangente al vacío fuera de los núcleos.
- **Producción:** 20 redes, ~3 h de CPU.
  - Conjunto A: 256³, 2 semillas RP³ y 4 U(1).
  - Conjunto B: 192³.
  - Conjunto C: núcleo más resuelto.
  - Conjunto D: paso de tiempo.
  - Conjunto E: condiciones iniciales con ξ_i = 0.25 y 1.0.
- **Ventana de análisis:** ℓ = ln(m/H) ∈ [3, 4.85].

**Validaciones.**
- **V2 (espectro de masas):** reproduce {0×3, 0.24×6, 1.84} con error ≤ 6×10⁻⁵.
- **V3 (detector de cuerdas):** exacto.
- **V4 (estimador sin cuerdas):** a⁴ρ constante al 2.9 %, con 0 cuerdas.
- **V1 (energía):** falló tal como estaba preregistrado. Es un desfase constante ∝ dτ² de la energía sombra del leapfrog, sin deriva secular (< 10⁻³). Se mantuvo dτ = 0.25 y se añadió el conjunto D (enmienda 1). El efecto sobre R_Γ es 0.016.

**Auditoría independiente** (código propio):
- fuerza, ecuaciones de movimiento, hessiano, detector y longitud de cuerda: PASS;
- **mismo estimador en los dos modelos**: un mismo paquete de Goldstone da exactamente 1.0000000;
- R_Γ recalculado exactamente.

**Bug encontrado:** la prueba de tendencia de R_Γ con ℓ ajustaba sobre un NaN, así que la rama B3 preregistrada no se evaluaba. Corregido: la tendencia es −0.48 por unidad de ℓ (significativa con el bootstrap del código) y **B3 se evalúa ahora**. No cambia el veredicto.

## 2. El número

| Magnitud | Valor |
|---|---|
| **R_Γ** (emisión RP³ / U(1), mismo ℓ) | **1.69 ± 0.29** (estadístico 0.04; sistemático dominante: condición inicial, 0.27) |
| Estimador integral (enmienda 5) | 1.67 ± 0.04 |
| R_ξ (densidad de cuerdas) | 1.56 ± 0.09 |
| η (emisión por unidad de energía de cuerda) | RP³ 0.60, U(1) 0.54 |
| ξ en ℓ = 4.5 | RP³ 0.74, U(1) **0.45** |
| Radiación masiva / Goldstone (final) | RP³ 0.98, U(1) 0.72 |
| Pendiente espectral q (RP³, en 10H < k < m/4) | 0.76–0.97, cercana a log-plana |

| Sistemático | ΔR_Γ |
|---|---|
| Caja (B) | 0.057 |
| Resolución del núcleo (C) | 0.024 |
| Paso de tiempo (D) | 0.016 |
| Condición inicial (E: 0.08 / 0.27) | dominante |
| Máscara (r_m = 1.5 / 4) | ≤ 0.005 |

**Lectura física.**
- La red RP³ es más densa, 1.56 veces. Los defectos ℤ₂ se aniquilan por pares en lugar de sumar devanados, y aun así sobrevive más longitud.
- Por unidad de energía de cuerda, emite casi igual que la U(1). No hay ningún mecanismo que la haga «silenciosa».
- Según la auditoría, si se restara la energía del baño masivo que acaba convirtiéndose en Goldstones, R_Γ subiría un ~9 %, es decir, más lejos del cierre.

## 3. Por qué INCONCLUSO

La puerta de control (ESPEC §6) pedía ξ_U1(ℓ = 4.5) ∈ [0.5, 2.5]. Las cuatro semillas de A dan 0.43–0.50 (media 0.45); B y E100 dan 0.52–0.54.

La figura muestra la causa, y la auditoría la confirma:
1. La red empieza sobredensa.
2. Colapsa hacia ℓ ≈ 3 y se pasa de largo: ξ ≈ 0.36 en U(1).
3. En la ventana vuelve a subir despacio (dξ/dℓ ≈ 0.07–0.20, frente a ≈ 0.24 en la literatura).

Es decir, **todavía no ha llegado al atractor de escalado**. El cociente se ha medido sobre redes que se están aproximando a ese estado, no sobre redes ya en escalado. Su dependencia de la condición inicial, 0.27, ya está dentro de σ.

## 4. Cierre con la ε medida (S6)

ε en el ℓ físico (≈ 50–85), según la rama:

| Rama | Hipótesis de la literatura | ε_RP³ |
|---|---|---|
| B1 | GHV: ξ ≈ 0.24ℓ | R_Γ ≈ 1.7 |
| B2 | Hindmarsh et al.: ξ\* = 1.20 | R_Γ·1.2/(0.24ℓ) ≈ 0.10–0.34 |
| B3 | Tendencia medida, R = c0 + c1/ℓ | ≈ 0.05–0.28 |

U1 pedía ε ≲ 10⁻³–6×10⁻².

**Puertas preregistradas (P1–P5 de U1):** las tres ramas quedan **excluidas** en los 12 nodos con K_sw = 0, incluso a −1σ y siempre por P5, la tibieza de lo emitido. Sin la puerta de control, el veredicto sería FALSADA_EN_FORMA_ACTUAL.

**Para pasar haría falta** que la red emitiera muchas menos veces de lo medido (con R central y K_sw = 0):
- 200–60 veces menos en B1;
- 14–3.5 veces menos en B2;
- 8–2 veces menos en B3.

**Exploratorio 1: Lyman-α para materia oscura mixta** (no preregistrado).
- Se sustituye la franquicia fija del 10 % de P5 por la cota de Garcia-Gallego et al. 2025 (arXiv:2504.06367): f_WDM < 0.14 (m/keV)^1.1. La masa equivalente se toma del mismo integrador de flujo libre.
- **B2 y B3 sobreviven con el espectro log-plano**, incluso a ±1σ. B1 sigue excluida.
- Predicción en ese caso:
  - m_G ≈ 2–7×10⁻¹⁵ eV;
  - entre el 15 % y el 54 % de la materia oscura es tibia, con un λ_fs como el de un WDM térmico de 2.2–4.4 keV;
  - ΔN_eff ≈ 0.002–0.006.
- **Una supresión de pequeña escala en el bosque Lyα, en el borde de las cotas actuales, sería entonces una predicción falsable.**

**Exploratorio 2: Goldstones sin masa** (sin m_G, con la materia oscura de otro origen). La red produce radiación oscura:
- **B1:** ΔN_eff(CMB) ≈ 0.22–0.31, al borde de Planck (0.29). CMB-S4 lo decidiría.
- **B2:** ≈ 0.015.

## 5. Estado y siguiente paso

- **El eslabón 8 deja de ser NOT_CALIBRATED y pasa a medido en cociente:** R_Γ = 1.69 ± 0.29 a ℓ ≈ 3–5. Falta el régimen de escalado.
- **La ballena 5 en su forma actual queda en tensión fuerte.** Solo sobrevive si se cumplen tres condiciones a la vez:
  1. las redes escalan con ξ constante (la disputa GHV / Hindmarsh, abierta en la literatura);
  2. el espectro es log-plano;
  3. Lyα admite una fracción tibia.
- **Siguiente prueba** (cierra INCONCLUSO, que es lo que pide la puerta):
  - (a) Redes que lleguen al atractor: fase de preparación con cuerdas gruesas (PRS) o con crecimiento del núcleo, y después cuerdas físicas; N ≥ 384 si la máquina WSL lo permite.
  - (b) La puerta Lyα mixta, preregistrada con la distribución f(p) real y una transferencia de Boltzmann (CLASS), en vez del indicador λ_fs.
  - (c) El espectro q a mayor ℓ.

## 6. Reproducir

```
cd network
python validate.py                                       # V1–V4 (~15 s)
python produce.py D; python produce.py A_U1; python produce.py A_RP3a
python produce.py E; python produce.py B; python produce.py C   # ~3 h de CPU con 2 núcleos
python analyze.py                                        # U2_RESULT.json + figura
cd ..; python -B chain/s6_measured.py                      # cierre con ε medida (~4 min)
```

Necesita gcc con OpenMP. Cola WSL: `make network-u2`.
