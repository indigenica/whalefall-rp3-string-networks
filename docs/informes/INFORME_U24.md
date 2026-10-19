# U24 · El drenaje del baño masivo no cambia R_Γ (pero sí toca a las cuerdas)

*19 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U24.md` (commit `55070bd`, 18 oct, 14:20) y su enmienda 1 (18 oct, 20:53), escrita tras la semilla 61 y antes de ver las semillas 62 y 63.
- *Código:* `network/lattice.py` (`damp_rate`), `network/u24_produce.py` y `network/u24_analyze.py`.
- *Resultado:* `network/results/u24/U24_RESULT.json`; las redes, en `network/results/u24/runs/`.
- *Ejecución:* en la nube con 4 hilos, 30 redes de 384³ en 16.9 h de cálculo. Hubo dos reinicios del contenedor, que obligaron a rehacer dos veces una misma red (`produce.out`).

## 0. En una frase

**Veredicto preregistrado: DRENAJE_ROBUSTO.**
- R_Γ cambia como mucho un **1.9 % ± 1.4 %** al apagar el drenaje antes o después, o al hacerlo la mitad o el doble de fuerte.
- El umbral era |δ| ≤ 5 % con σ ≤ 3 % en las cuatro variantes.

**La dependencia de R_Γ con el drenaje queda acotada en el 5 % preregistrado:** el máximo es 1.9 % ± 1.4 %, y su borde a 2σ, 4.7 % (5.5 % con el error entre semillas). Pero el drenaje sí toca a la red: con el doble de intensidad, R_ξ baja un 5.6 % ± 1.5 %, y la emisión por unidad de longitud (R_Γ/R_ξ) sube ≈ 8 %.

## 1. Las cuatro variantes frente a la referencia D0

Son 3 semillas (61, 62, 63) por variante y modelo. δ = X^v/X^D0 − 1 en la ventana [4.36, 5.21], con bootstrap emparejado por semilla (el preregistrado). Entre paréntesis va el error entre semillas (ddof = 1, entre √3), que es informativo (enmienda 1).

| Variante | ℓ_d | Γ_d | δR_Γ | δR_ξ | δε_RP³ | δε_U(1) |
|---|---|---|---|---|---|---|
| Dearly | 3.4 | m_r | −0.2 % ± 0.3 % (0.4) | +0.6 % ± 0.2 % | −1.9 % | −1.8 % |
| Dlate | 4.1 | m_r | +0.3 % ± 0.4 % (0.6) | −0.1 % ± 0.1 % | +2.7 % | +2.3 % |
| Dweak | 3.8 | m_r/2 | **+1.9 % ± 1.4 %** (1.8) | +2.5 % ± 2.2 % | −6.7 % | −8.4 % |
| Dstrong | 3.8 | 2 m_r | **+1.8 % ± 1.0 %** (1.1) | **−5.6 % ± 1.5 %** | +7.8 % | +6.0 % |

- **El veredicto también se cumple con el error entre semillas**, que es algo mayor: como mucho 1.8 %, frente al umbral del 3 %.
- **El momento en que se apaga el drenaje no importa para el cociente** (Dearly y Dlate, |δ| ≤ 0.3 %).
- **Para las emisiones de cada modelo sí importa.** Apagarlo después sube ε_RP³ y ε_U(1) un +2.7 % y un +2.3 % (6σ y 8σ), y apagarlo antes las baja ≈ 1.9 %; ξ se mueve ≈ ±2 % en el mismo sentido.
- **Es la lectura informativa que pedía el ESPEC, pero con el signo contrario:** si quedara baño, Dlate bajaría ε. Así que el momento del drenaje actúa sobre las cuerdas, no sobre un baño residual.
- **La intensidad sí mueve las emisiones de cada modelo un 6–8 %,** pero en el mismo sentido en RP³ y en U(1), así que el cociente cambia menos de un 2 %.

## 2. Lo que el drenaje hace a las cuerdas (informativo)

La expectativa del preregistro, «R_ξ cambia < 1 % como en U11», **no se cumple con la intensidad:**
- con el doble de intensidad, R_ξ baja un 5.6 % ± 1.5 % (3.7σ);
- con la mitad, sube un 2.5 % ± 2.2 %.

- **No es que las emisiones sigan a ξ** en cada modelo. Con el doble, ε_RP³ sube +7.8 % con ξ_RP³ +2.0 %, y ε_U(1) sube +6.0 % con ξ_U(1) +8.1 %.
- **R_Γ se mueve poco porque las emisiones de los dos modelos cambian en el mismo sentido.** La frase contraria de la enmienda 1, punto 5, escrita tras la semilla 61, no se confirma.
- **El mecanismo más probable** sigue siendo que el umbral de «cerca del vacío» del drenaje incluya la parte externa de los núcleos de cuerdas en movimiento.

**Consecuencia para el artículo:**
- R_Γ es robusto frente al drenaje, dentro del 5 % preregistrado.
- Ya no se puede usar «R_ξ no cambia» como prueba de que el drenaje no toca a las cuerdas.
- R_ξ lleva una sensibilidad al drenaje de ≈ 6 % por un factor 2 en la intensidad.
- La emisión por unidad de longitud (R_Γ/R_ξ = 0.97 desde abajo) lleva una de ≈ 8 % (+7.9 % ± 2.1 %, las tres semillas en el mismo sentido). El «13 %» del artículo la lleva dentro, y así se dice en la limitación (ii).

## 3. La referencia a 384³ frente a 576³ (enmienda 1, punto 4)

| | R_Γ en [4.36, 5.21] |
|---|---|
| D0, 384³, 3 semillas | **1.51 ± 0.14** (bootstrap de semillas; error entre semillas, 0.17; por semilla: 1.20, 1.62, 1.78) |
| U11, 576³, 6 + 4 semillas, misma ventana | 1.56 ± 0.09 |

- La diferencia es de 0.3σ: **no se detecta efecto de tamaño finito,** aunque con este error no se pueden excluir efectos de hasta ≈ 30 %. Lo más probable es que el 1.20 de la semilla 61 sola sea dispersión: a 384³ las tres semillas dispersan ≈ 20 %.
- R_ξ(D0) = 1.54, frente a 1.58 en U11 en la misma ventana (1.61 en la ventana completa).

## 4. Notas de registro

- **La ventana de ε tiene 9 puntos** ([4.36, 5.12]), y la de ξ, 10. Es lo mismo que en U11, y la enmienda 1 lo declara.
- **El bootstrap con 3 semillas** da un σ algo menor que el error entre semillas (factor ≈ √(2/3)). El veredicto se cumple con los dos.
- **Coste real:** 53–66 min por red RP³ y 9–11 min por red U(1); 16.9 h en total, frente a las 10 h estimadas.
- **El analizador** quedó fijado en `109ce51`, antes de la primera red terminada. Con una semilla da SIN_VEREDICTO (enmienda 1).

## 5. Lo que U24 no resuelve

- **No es 576³:** la ventana es más corta (0.85 en ℓ en vez de 1.25) y solo hay 3 semillas.
- **No prueba otros umbrales de «cerca del vacío»,** que es justo lo que parece tocar a las cuerdas.
- **El paso de tiempo** lo mide U25.
