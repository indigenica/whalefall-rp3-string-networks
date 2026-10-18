# U18 · Los estimadores de R_Γ coinciden con instantáneas densas; con 30 se equivocan en ambos sentidos

*14 de octubre de 2026 (Madrid).*
- *Preregistro:* `docs/preregistros/ESPEC_U18.md` (commit `c7ff41b`).
- *Producción:* `network/u18_produce.py`, ~3.5 h en la nube. *Análisis:* `network/u18_analyze.py`. *Resultado:* `network/results/u18/U18_RESULT.json`.
- *Corrección de código:* la primera ejecución del análisis falló por un nombre de clave mal construido en la validación V1 (`R_Gamma_30` en lugar de `R_Gamma_local_30`). Se corrigió sin cambiar ningún criterio (commit `ac72485`).

## 0. En una frase

**Veredicto preregistrado: INDETERMINADO.** Con instantáneas densas, el cambio del estimador local con la resolución es −0.031: ni ≤ 0.010 (MUESTREO) ni ≥ +0.015 (FISICO).

Lo que sí queda claro es más importante que el veredicto:
- **Con 270 instantáneas, el estimador local y el integral coinciden a 0.4–0.8 %. Con las 30 de siempre discrepaban un 13–15 %.**
- a⁴ρ_a lleva una fluctuación resoluble (rms 3–4 %) que 30 instantáneas muestrean mal.
- Los dos estimadores de U16 estaban afectados, y el efecto de resolución que se veía (+1.6 % a r = 0.6) cambia de signo con muestreo denso.

## 1. Validación

**V1:** con solo las 30 instantáneas de U16, R_ξ, R_Γ local y R_Γ integral reproducen bit a bit los de la semilla 61 de U16, en las dos resoluciones. La evolución no depende de cómo se trocee ni de las medidas.

## 2. Resultados (semilla 61)

| | r = 1, 30 inst. | r = 1, densas | r = 0.6, 30 inst. | r = 0.6, densas |
|---|---|---|---|---|
| R_ξ | 1.500 | 1.501 | 1.496 | 1.495 |
| R_Γ local | 1.629 | **1.607** | 1.653 | **1.576** |
| R_Γ integral | 1.416 | **1.601** | 1.404 | **1.588** |
| \|local − integral\| | 0.213 | **0.006** | 0.249 | **0.012** |

**Cambios con la resolución**, Δ = X(0.6) − X(1):

| Estimador | 30 instantáneas (= U16) | densas |
|---|---|---|
| R_Γ local | +0.023 | −0.031 |
| R_Γ integral | −0.012 | −0.013 |
| R_ξ | −0.004 | −0.006 |

**Fluctuación de a⁴ρ_a** en la ventana, respecto a un ajuste cúbico en ℓ:
- RP³: rms 4.2 % (r = 1) y 4.0 % (r = 0.6);
- U(1): 3.2 % y 2.8 %;
- unas 13 veces el ruido entre instantáneas vecinas, así que es resoluble;
- periodo dominante ≈ 0.31 en ℓ, el mismo en los dos modelos y las dos resoluciones.

**El estimador local «suavizado»** (regresión en ±0.15 en ℓ) da 1.80 y 1.74. Esa ventana abarca justo un periodo de la fluctuación, y una regresión lineal sobre un periodo completo de una oscilación tiene un sesgo grande. Es un defecto del diseño de ese estimador auxiliar: no se usa.

## 3. Lectura

- **La hipótesis de ESPEC_U18 era aliasing de una oscilación rápida** (periodo de 0.04–0.07 en ℓ). Lo que aparece es otra cosa: una fluctuación más lenta, con periodo ≈ 0.31 en ℓ, presente en los dos modelos. No se ha identificado su origen. Con 30 instantáneas (una cada Δℓ ≈ 0.12) se muestrea con dos o tres puntos por periodo y sesga los dos estimadores: el integral, que depende solo de los extremos de la ventana, hasta un 13 %, y el local unos pocos por ciento.
- **El efecto de resolución de U16 no es robusto.** Con muestreo denso, el estimador local da −2.0 % entre r = 1 y 0.6 (en U16, +1.6 %), y el integral −0.8 %. El +4 % de r = 0.5 (U16, enmienda 1) es, con toda probabilidad, del mismo origen; no se ha repetido con muestreo denso.
- **R_ξ es inmune:** es una medida instantánea, no una derivada.
- **Para los R_Γ publicados** (U6, U11 y U16, todos con unas 30 instantáneas): en esta semilla, el local con 30 instantáneas se aleja del denso un +1.4 % (r = 1) y un +4.9 % (r = 0.6). Las redes de 576³ tienen su propia fase de la fluctuación, así que no se puede corregir a posteriori. **Conviene añadir una sistemática de muestreo de ~±5 % a R_Γ local**, y no usar el estimador integral con muestreo escaso. Para U11 (1.56 ± 0.09), eso lleva el error a ~±0.12. No se reclasifica nada.

## 4. Lo que U18 no resuelve

- **Una sola semilla:** las cifras de §3 son de una realización.
- **El origen de la fluctuación de periodo ≈ 0.31 en ℓ.** Podría ser física (una respiración coherente de la red o del baño tras la condición inicial) o de la medida (máscaras). Habría que verla en redes con drenaje.
- **r = 0.5 con muestreo denso** (384³, ~4 h) queda pendiente.
- **Siguiente paso recomendado:** producir las redes futuras con ≥ 100 instantáneas en la ventana, o con medidas baratas de ρ_a entre instantáneas, y repetir R_Γ de U11 en WSL con muestreo denso.
