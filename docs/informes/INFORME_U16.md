# U16 · R_ξ y R_Γ no dependen de la resolución del núcleo

*13 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U16.md` (commit `c9a640c`, antes de producir ninguna red). El análisis, `network/u16_analyze.py`, se fijó en `4cb9c29`, antes de completar el primer bloque.
- *Producción:* `network/u16_produce.py`, 22 redes en la nube (4 vCPU), ~6 h en total. *Resultado:* `network/results/u16/U16_RESULT.json`.
- *Incidencia:* el contenedor se reinició durante la RP³ de volumen de la semilla 62. Se relanzó el mismo comando: las 20 redes hechas se saltaron y la interrumpida se repitió desde el principio (`produce.log`). No cambia nada.

## 0. En una frase

**Veredicto preregistrado:**
- **Resolución:** RESOLUCION_ROBUSTA para R_ξ y para R_Γ.
- **Volumen:** VOLUMEN_ROBUSTO para R_Γ e INDETERMINADO para R_ξ, por falta de potencia, sin efecto visible.

Con el mismo campo inicial en las tres mallas, pasar de m_r·a·dx = 1 a 0.6 cambia R_ξ en −0.26 % y R_Γ en +1.6 %. El criterio de resolución de todas nuestras redes (m_r·a·dx = 1 al final) no sesga los cocientes más allá de un par de por ciento.

## 1. Brazo de resolución

ℓ_fin = ln 96 = 4.56, L·H·a = 2 al final, ventana [3.26, 4.51]. Hay tres semillas, con RP³ y U(1) emparejadas en las tres resoluciones.

| r = m_r a_fin dx | N | R_ξ | R_Γ | R_Γ integral | ξ_RP³ | ξ_U(1) | ε_RP³ | ε_U(1) |
|---|---|---|---|---|---|---|---|---|
| 1 | 192 | 1.567 | 1.605 | 1.381 | 0.641 | 0.409 | 0.345 | 0.215 |
| 0.75 | 256 | 1.566 | 1.661 | 1.379 | 0.638 | 0.407 | 0.353 | 0.213 |
| 0.6 | 320 | 1.563 | 1.631 | 1.373 | 0.639 | 0.409 | 0.354 | 0.217 |

**Diferencias emparejadas Δ = X(0.6) − X(1)**, media de las tres semillas ± error estándar:

| | Δ | relativo | por semilla | extrapolación al continuo (r² / r) |
|---|---|---|---|---|
| R_ξ | −0.0040 ± 0.0007 | −0.26 % | −0.0039, −0.0052, −0.0029 | 1.562 / 1.558 |
| R_Γ | +0.026 ± 0.002 | +1.6 % | +0.023, +0.029, +0.026 | 1.667 / 1.696 |
| R_Γ integral | −0.009 ± 0.005 | −0.6 % | −0.012, +0.002, −0.016 | 1.370 / 1.362 |
| ξ_RP³ | −0.0019 ± 0.0007 | −0.30 % | | |
| ε_RP³ | +0.0089 ± 0.0006 | +2.6 % | | |
| ε_U(1) | +0.0020 ± 0.0002 | +0.9 % | | |

**Lectura:**
- **Las condiciones emparejadas funcionan.** Las tres semillas dan diferencias del mismo signo y tamaño, y el error de Δ es 50–100 veces menor que la dispersión entre semillas (R_Γ va de 1.52 a 1.68 según la semilla).
- **La densidad de cuerdas no depende de la resolución:** ξ cambia menos de un 0.3 % en los dos modelos.
- **La emisión sube un poco al resolver mejor**, más en RP³ (+2.6 %) que en U(1) (+0.9 %). R_Γ sube un 1.6 %, de forma significativa pero por debajo del 5 % del criterio.
- **No es monótono:** el máximo está en r = 0.75. Por eso las extrapolaciones al continuo (+4 % en r², +6 % lineal) no son fiables, y no se usan para corregir nada. Si se quiere una cota, el sesgo de R_Γ con r = 1 está entre 0 y −6 %, del lado de subestimar la emisión de RP³.
- **El estimador integral** no muestra el efecto (−0.6 % ± 0.3 %).

## 2. Brazo de volumen

Con r = 1 y ℓ_fin = 4.56, se compara L = 3.33 τ_fin (N = 320, dos semillas, condiciones normales) con L = 2 τ_fin (N = 192, las tres semillas del brazo anterior).

| | L = 2 τ_fin | L = 3.33 τ_fin | Δ | σ | Clase |
|---|---|---|---|---|---|
| R_ξ | 1.567 | 1.558 | −0.009 (−0.6 %) | 0.087 (5.5 %) | INDETERMINADO |
| R_Γ | 1.605 | 1.615 | +0.011 (+0.7 %) | 0.078 (4.9 %) | VOLUMEN_ROBUSTO |
| R_Γ integral | 1.381 | 1.441 | +0.061 | — | (informativo) |

- Por semilla, la caja grande da ξ_RP³ = 0.603 y 0.623, frente a 0.586–0.683 de la pequeña. La caja grande dispersa menos, como se espera de 4.6 veces más volumen de Hubble.
- **No hay efecto visible del volumen**, pero con dos y tres semillas no emparejadas la prueba solo descarta efectos mayores que ~10 % (2σ).

## 3. Qué significa

- **El +16 % de R_ξ del conjunto C de U2** (auditoría U2, hallazgo 8; m_r·a·dx ≤ 0.6 con N fijo) no viene de la resolución: aquí, a L·H fijo, R_ξ cambia un 0.26 %. Tampoco se ve en el brazo de volumen, aunque ese brazo solo descarta ≳ 10 %. Lo más probable es que fuera estadístico (C tenía una o dos semillas) o una mezcla de volumen y ventana. **Queda sin explicar del todo, y se dice.**
- **Para el artículo de la red:** se puede afirmar que R_ξ y R_Γ son robustos frente a la resolución del núcleo a ≤ 2 % (de m_r·a·dx = 1 a 0.6), con condiciones iniciales emparejadas. Para el volumen la prueba es más débil: con R_Γ se descartan efectos mayores que ~10 %, y para R_ξ el resultado es INDETERMINADO.
- **Para el sector oscuro**, el veredicto es RESOLUCION_ROBUSTA, así que según ESPEC_U16 no se corrige nada y U11 no se reclasifica. A título informativo:
  - con el cambio medido (+1.6 %), U11 queda en R_Γ + 1σ ≈ 1.585 + 0.09 = 1.675 ≤ 1.69, todavía en LIMPIA_ZONA_PERMITIDA;
  - con las extrapolaciones al continuo, que no son fiables (+4 % y +6 %), quedaría en 1.72–1.75, en LIMPIA_INTERMEDIA según el criterio de U11. Es la zona entre las filas de U7c que pasan (R_Γ = 1.69) y la esquina que falla (2.09 con x₀ = 30).
  - Para fijar R_Γ en el continuo a mejor que un 2 % haría falta una cuarta resolución (r ≈ 0.5, N = 384, dentro de lo que cabe en la nube).

## 4. Lo que U16 no resuelve

- **ℓ ≤ 4.56:** la ventana está en la aproximación al escalado y sin drenaje. Se supone que el sesgo de resolución depende de r y no de ℓ.
- **No mide q** (menos de 8 bins en la ventana).
- **El volumen**, con poca potencia: para bajar σ a ~2 % harían falta unas 6 semillas de cada caja, unas 4 h más en esta máquina, o pocas en WSL.

## 5. Enmienda 1 · cuarta resolución, r = 0.5 (384³)

*Preregistrada en ESPEC_U16, enmienda 1 (13 de octubre, 09:35), antes de lanzarla. Informativa: el veredicto de §0 no cambia. Seis redes, ~2.3 h, con el mismo campo inicial rellenado desde 320³ (comprobado a 3×10⁻⁸).*

| r = m_r a_fin dx | 1 | 0.75 | 0.6 | 0.5 |
|---|---|---|---|---|
| R_ξ | 1.567 | 1.566 | 1.563 | 1.564 |
| R_Γ (local) | 1.605 | 1.661 | 1.631 | 1.669 |
| R_Γ (integral) | 1.381 | 1.379 | 1.373 | 1.374 |
| ξ_RP³ / ξ_U(1) | 0.641 / 0.409 | 0.638 / 0.407 | 0.639 / 0.409 | 0.638 / 0.408 |
| ε_RP³ / ε_U(1) | 0.345 / 0.215 | 0.353 / 0.213 | 0.354 / 0.217 | 0.359 / 0.215 |

Por semilla, R_Γ (local) en r = 1, 0.75, 0.6 y 0.5:
- semilla 61: 1.629, 1.675, 1.652, 1.687;
- semilla 62: 1.516, 1.569, 1.545, 1.580;
- semilla 63: 1.679, 1.752, 1.705, 1.751.

**Diferencias emparejadas Δ = X(0.5) − X(1):**
- R_ξ: −0.0031 ± 0.0012 (−0.20 %);
- R_Γ local: **+0.065 ± 0.004 (+4.0 %)**;
- R_Γ integral: −0.007 ± 0.006 (−0.5 %);
- ε_RP³: +4.0 %; ε_U(1): 0.0 %.

La clasificación de ESPEC_U16 aplicada a r = 0.5, a título informativo, sigue siendo RESOLUCION_ROBUSTA para los dos (|Δ| < 5 %), pero para R_Γ local por poco.

**Lectura:**
- **La densidad de cuerdas está convergida:** R_ξ y ξ de cada modelo cambian ≤ 0.5 % en todo el recorrido, de r = 1 a 0.5.
- **La emisión local de RP³ sube al resolver mejor**, +4 % entre r = 1 y 0.5. La de U(1) no cambia. Encima hay un zigzag de ±2 % (bajada en r = 0.6) que se repite igual en las tres semillas: no es ruido, sino algo sistemático del estimador o del muestreo.
- **El estimador integral no lo ve** (−0.5 %). El integral cuenta la energía total emitida en la ventana, Δ(a⁴ρ_a); el local promedia derivadas centradas entre instantáneas. Que solo cambie el local apunta a que al menos parte del efecto está en cómo cae la derivada respecto a las instantáneas, que se colocan en los mismos ℓ pero con pasos dτ distintos en unidades físicas. No se ha separado.
- **Extrapolación al continuo con cuatro puntos:** R_Γ local = 1.679 (en r²) y 1.714 (lineal), es decir +4.6 % y +6.8 % sobre r = 1. R_ξ = 1.562 y 1.560. Con el zigzag, estas extrapolaciones solo dan el orden de magnitud.
- **Para U11, a título informativo** (no se reclasifica: la enmienda no lleva veredicto):
  - con +4 %, R_Γ = 1.62 ± 0.09, y +1σ = 1.71 supera por poco el 1.69 de LIMPIA_ZONA_PERMITIDA (sería LIMPIA_INTERMEDIA);
  - con el estimador integral no cambia nada.
  - **Para el artículo:** R_ξ es robusto a resolución a ≤ 0.5 %; R_Γ local puede llevar un sesgo de −4 a −7 % a m_r·a·dx = 1 en RP³, que el estimador integral no muestra. Merece una frase.

## Nota · 14 de octubre de 2026: U18 revisa el efecto en R_Γ

U18 (`docs/informes/INFORME_U18.md`) repitió la semilla 61 con 270 instantáneas en lugar de 30. Con muestreo denso, los estimadores local e integral de R_Γ coinciden a < 1 %, y el cambio con la resolución pasa de +1.6 % a −2.0 % (local) y −0.8 % (integral). **El aumento de R_Γ local de §1 y §5 se debe, al menos en buena parte, al muestreo con 30 instantáneas, no a la resolución del núcleo.** Lo que sí se mantiene es que R_Γ depende de m_r·a·dx en ≤ ~2 % y R_ξ en ≤ 0.5 %. El veredicto RESOLUCION_ROBUSTA no cambia.
