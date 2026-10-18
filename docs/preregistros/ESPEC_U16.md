# U16 · Resolución y volumen de la red: R_ξ y R_Γ hacia el continuo · Preregistro

*12 de octubre de 2026, 22:05 (Madrid). Escrito antes de producir ninguna red de U16. Se ejecuta en la nube (4 vCPU, 16 GB).*

## Motivo

Todas las redes de U4–U11 terminan con m_r·a·dx = 1: al final, el núcleo de la cuerda mide una celda. Es el criterio habitual en la literatura de cuerdas físicas, pero es justo lo primero que pregunta un árbitro: **¿cuánto cambian R_ξ y R_Γ si se resuelve mejor el núcleo?**

Lo que se sabe antes:
- **U2, conjunto C** (auditoría U2, hallazgo 8): con m_r·a·dx ≤ 0.6 a 256³, R_ξ subió de 1.52 a 1.75 (+16 %) y R_Γ cambió solo un 1.4 %. Pero en C, al bajar la resolución con N fijo, también se encogió la caja en unidades del horizonte para cada ℓ, así que el efecto mezcla **resolución y volumen**.
- El borrador para PRD dice que R_ξ ≈ 1.6 es robusto en todos los protocolos, pero no ensaya la resolución.

U16 separa los dos efectos:
- **brazo de resolución:** L·H fijo para cada ℓ y m_r·a_fin·dx = r variable;
- **brazo de volumen:** r = 1 fijo y la caja más grande.

## Convenios (los de `network/lattice.py`)

- a = sτ con s = r/(m_r τ_fin), así que m_r·a·dx = r en τ_fin y ℓ = ln(r τ²/τ_fin).
- dx = 1, dτ = 0.25 (el cociente dτ/dx no cambia; en unidades físicas el paso de tiempo baja con r).
- La caja mide N, y el horizonte comóvil vale τ. Con τ_fin = N/2, L·H·a = N/τ = 2 al final, igual que en U4–U11, y ℓ_fin = ln(rN/2).

## Producción (`network/u16_produce.py`)

Protocolo de U6: sin drenaje, ξ_i = 0.5, arranque en ℓ_i = 1 y 30 instantáneas.

**Brazo de resolución**, ℓ_fin = ln 96 = 4.564 en las tres resoluciones:

| N | τ_fin | r = m_r a_fin dx |
|---|---|---|
| 192 | 96 | 1 |
| 256 | 128 | 0.75 |
| 320 | 160 | 0.6 |

- **Condiciones iniciales emparejadas.** Para cada semilla se genera un único campo físico: ruido blanco en 320³ filtrado con la longitud de correlación de la red de 320³, y truncado en Fourier a 256³ y 192³. El filtro deja en el corte de 192³ una amplitud < 10⁻¹⁰, así que las tres redes empiezan con el mismo campo y solo se diferencian por la discretización. En RP³, la dirección n es el autovector principal del campo matricial truncado, y en U(1) la fase del campo truncado, como en `lattice.initial`.
- Las comparaciones se hacen **por semilla**, emparejadas, lo que quita la mayor parte del ruido entre semillas.
- **Semillas 61, 62 y 63**, en RP³ y en U(1): 18 redes.

**Brazo de volumen**, r = 1 y ℓ_fin = 4.564:
- N = 320 con τ_fin = 96, es decir, L = 3.33 τ_fin en lugar de 2 τ_fin: una caja con 4.6 veces más volumen de Hubble.
- Condiciones iniciales normales de `lattice.initial` (no emparejables con 192³), semillas 61 y 62, en RP³ y en U(1): 4 redes.

**Coste estimado** con el banco de pruebas de esta máquina (paso RP³ de 0.46 s a 192³ con 4 hilos, ∝ N³): unas 5–6 h en total. Orden: semilla por semilla en el brazo de resolución (192, 256, 320; RP³ y después U(1)) y luego el brazo de volumen.

## Análisis (`network/u16_analyze.py`, con los estimadores de U6: `analyze.window_ratio` e `int_ratio`)

**Ventana:** [ℓ_fin − 1.3, ℓ_fin − 0.05] = [3.26, 4.51], con la rejilla de 10 puntos de U6.

**Medidas**, para cada (N, r) y cada semilla:
- ξ y ε_lat de RP³ y de U(1) promediados en la ventana;
- R_ξ, R_Γ (diferencias centradas) y R_Γ integral.

**Cambio con la resolución:** Δ_X = X(r = 0.6) − X(r = 1), por semilla (las tres semillas de RP³ emparejadas con las tres de U(1) de la misma semilla); σ_Δ es el error estándar de la media de las tres diferencias.

**Extrapolación al continuo** (informativa): ajuste de X(r) = X₀ + X₂ r² y de X(r) = X₀ + X₁ r a los tres valores medios. Se dan los dos X₀.

**Clasificación**, por separado para R_ξ y para R_Γ:
- **RESOLUCION_ROBUSTA:** |Δ_X| < 0.05·X(r = 1), y además σ_Δ < 0.05·X(r = 1) (la prueba tiene potencia para ver un 5 %).
- **SESGO_DE_RESOLUCION:** |Δ_X| ≥ 0.05·X(r = 1) y |Δ_X| ≥ 2σ_Δ.
- **INDETERMINADO:** en otro caso.

**Cambio con el volumen:** Δ_V = X(L = 3.33 τ_fin) − X(L = 2 τ_fin) a r = 1, con las semillas no emparejadas: σ_V sale de la dispersión entre semillas de cada brazo. Clasificación igual, con los nombres **VOLUMEN_ROBUSTO**, **SESGO_DE_VOLUMEN** e **INDETERMINADO**.

**Consecuencias, fijadas ahora:**
- Si R_ξ o R_Γ salen con SESGO_DE_RESOLUCION, el valor a 576³ (r = 1) lleva ese sesgo. El informe dará el valor corregido al continuo con la extrapolación en r², y la diferencia entre las dos extrapolaciones como sistemática. El artículo para PRD deberá incluirlo.
- **Para el sector oscuro** solo importa R_Γ. Si el R_Γ de U11 corregido al continuo pasa de 1.69 + 1σ, U11 dejaría de estar en LIMPIA_ZONA_PERMITIDA según su propio criterio. Eso se informará, pero U11 no se reclasifica sin una enmienda suya.

## Lo que U16 no resuelve

- **ℓ ≤ 4.56.** La ventana está en la aproximación al escalado, no en el atractor, y antes del drenaje de U10a/U11. Se supone que el sesgo de resolución depende sobre todo de r y no de ℓ; no se comprueba.
- **No mide q:** a ℓ ≈ 4.5 la ventana [10H, m_r/4] tiene menos de 8 bins.
- **Tres semillas.** Con condiciones emparejadas basta para un efecto del 5 %, si las diferencias por semilla son coherentes; si no, el veredicto será INDETERMINADO.

## Enmienda 1 · 13 de octubre de 2026, 09:35 (Madrid): cuarta resolución, r = 0.5

*Escrita después de leer el resultado de U16 (`docs/informes/INFORME_U16.md`): RESOLUCION_ROBUSTA, con R_Γ cambiando +1.6 % de forma no monótona (máximo en r = 0.75). El equipo pidió la cuarta resolución para ver si ese cambio se estabiliza.*

- **Nueva fila del brazo de resolución:** N = 384, τ_fin = 192, r = 0.5 (ℓ_fin = ln 96, L·H·a = 2 al final, como las demás), semillas 61, 62 y 63, RP³ y U(1): 6 redes, ~5 h en la nube.
- **Condiciones iniciales emparejadas:** el mismo campo físico de las otras tres resoluciones. Se genera en 320³ como siempre y se rellena con ceros en Fourier hasta 384³. Es exacto porque el campo está limitado en banda: el filtro vale < 10⁻¹⁰ mucho antes del corte de 320³.
- **El veredicto preregistrado no cambia** (se decide con Δ = X(0.6) − X(1), y ya está dado). Con la nueva fila se informa, sin veredicto nuevo:
  - Δ₀.₅ = X(0.5) − X(1) por semilla, con la misma clasificación aplicada a título informativo;
  - las extrapolaciones al continuo con cuatro puntos (en r² y en r);
  - si R_Γ(r) es monótona entre 0.75 y 0.5.
- **RAM:** una RP³ de 384³ necesita ≈ 12 GB, de los 15.6 GB de la máquina; las redes van de una en una.
