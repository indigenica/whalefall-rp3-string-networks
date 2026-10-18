# U9 · ¿Es transitorio el crecimiento de la red RP³? Condiciones iniciales densas · Preregistro

*9 de octubre de 2026, ≈17:35 (Madrid). Escrito antes de producir las redes de U9.*

## Motivo

- En U6, la densidad de cuerdas de RP³ crece deprisa en la ventana: dξ/dℓ = 0.26–0.37 en [4.66, 5.36], con ξ_RP3 ≈ 0.9–1.0. La de U(1) crece a 0.09–0.12.
- R_Γ no tiene tendencia, pero R_ξ sube.
- Correia, Hindmarsh et al. (arXiv:2512.13653) atribuyen este tipo de crecimiento a redes que empiezan con poca densidad y se acercan al atractor desde abajo. Todas nuestras redes RP³ empezaron con ξ_i = 0.5.
- U9 repite RP³ con condiciones iniciales densas (ξ_i = 2), igual que la puerta G2′ de U4b hizo con U(1).

## Producción

`u4b_produce.py 5`, a 576³, con el mismo kernel y los espectros en cada instantánea, como en U6:
- RP³ con ξ_i = 2.0, semillas 61 y 62;
- U(1) con ξ_i = 2.0, semillas 61 y 62.

## Análisis (`network/u9_analyze.py`, fijado antes de los datos)

ℓ_ref = ℓ_fin − 0.3 = 5.36; ventana de pendiente [ℓ_ref − 0.7, ℓ_ref].
- **M_RP3**, memoria de la condición inicial de RP³ en ℓ_ref: |ξ̄_denso − ξ̄_disperso|/media. El disperso es la media de las seis semillas con ξ_i = 0.5 de U6; el denso, la media de las dos de U9.
- **Pendientes dξ/dℓ de RP³:** s_denso (media de 2) frente a s_disperso (media de 6).
- **Clasificación:**
  - **TRANSITORIO_DESDE_ABAJO:** M_RP3 < 20 % y s_denso ≤ 0.5 · s_disperso;
  - **CRECIMIENTO_GENUINO:** M_RP3 < 20 % y s_denso ≥ 0.7 · s_disperso;
  - **INDETERMINADO:** en otro caso, incluido M_RP3 ≥ 20 %.
- **Informativos:**
  - R_Γ y R_ξ de las redes densas (RP³ denso frente a U(1) denso), en la ventana de U4b, comparados con U6 (1.389 ± 0.062 y 1.62);
  - la memoria de U(1) a 576³;
  - q_RP3 de las redes densas.

## Lectura

- **TRANSITORIO:** la medida de R_Γ de U6 corresponde al régimen que escala, y la red RP³ se comporta como dicen Correia et al. para U(1).
- **CRECIMIENTO_GENUINO:** R_ξ seguiría subiendo con ℓ. Habría que ver si R_Γ también lo hace a ℓ mayores; a 576³ no se puede.

---

## Enmienda 1 (9 de octubre, ≈20:45 Madrid) · análisis informativos añadidos tras una revisión externa

**Qué se había visto al escribirla.** Solo la primera red RP³ densa (semilla 61), comparada con las seis poco densas de U6, en una lectura provisional:
- densidad en ℓ_ref: 1.08 frente a 0.99, memoria del 8.5 %;
- pendiente dξ/dℓ: 0.312 frente a 0.313 de media;
- emisión de la red densa 1.9 veces la de las poco densas en la ventana de U4b, con q ≈ 0.64;
- radiación masiva fuera de los núcleos, ρ_M/ρ_a ≈ 1.4 en la red densa frente a ≈ 1.1 en la poco densa al comienzo de la ventana.

No se había visto la segunda RP³ densa (semilla 62) ni ninguna U(1) densa.

**Qué no cambia.** La clasificación preregistrada (TRANSITORIO_DESDE_ABAJO / CRECIMIENTO_GENUINO / INDETERMINADO), sus umbrales y `network/u9_analyze.py` quedan como estaban.

**Definición de las etiquetas, para los informes.** Las etiquetas describen solo la ventana simulada.
- CRECIMIENTO_GENUINO significa: *crecimiento persistente dentro de la ventana numérica estudiada (ℓ ≤ 5.66), independiente de la densidad inicial de cuerdas ensayada*. No significa crecimiento demostrado hasta ℓ ≈ 70.
- TRANSITORIO_DESDE_ABAJO significa que la pendiente depende de la densidad inicial en esa ventana. No demuestra que ξ sea constante a ℓ grande.
- Ninguna de las dos refuta ni confirma la lectura de Correia et al. a ℓ grande.

**Análisis informativos nuevos** (`network/u9_extra.py`, salida `results/u4/U9_EXTRA.json`; no dan veredicto ni cambian la clasificación):

1. **Distribución por semilla.** ξ(ℓ_ref) y la pendiente de cada red densa, frente al rango y la media ± desviación de las seis poco densas.
2. **¿Desde arriba?** Para cada red densa: el ℓ del mínimo de ξ tras el desplome inicial (ℓ > 2.5) y si después ξ decrece en algún tramo de la ventana. Si nunca decrece, el informe dice que la prueba **no** ensayó una aproximación desde arriba.
3. **Historias de energía**, por red y en función de ℓ, fuera de los núcleos (máscara rm2.5): a⁴ρ_a, a⁴ρ_M y a³ρ_M.
4. **Cota de conversión masiva → Goldstone.**
   - C(ℓ) = parte negativa de H · a⁻³ d(a³ρ_M)/d ln τ, con la misma derivada local de tres instantes que Γ en `analyze.derived`.
   - Es una cota superior de la energía que la radiación masiva cede al sector Goldstone, porque atribuye a conversión todo el corrimiento al rojo por encima del no relativista.
   - Γ_corr = Γ − C es una cota inferior de la emisión de las cuerdas.
   - Se informan ε_lat y ε_corr, y los cocientes denso/poco denso para RP³ y para U(1), así como R_Γ con ε_corr.
5. **Fin de la relajación** ℓ_relax: el primer ℓ a partir del cual la media densa de ρ_M/ρ_a queda a menos del 10 % de la media poco densa del mismo modelo, y sigue así hasta el final. Si no ocurre antes de ℓ_fin − 0.3, se informa «sin relajación en la ventana».
6. **Barrido de la ventana.**
   - R_Γ (RP³ densa / U(1) densa, y RP³ poco densa / U(1) poco densa como referencia), R_Γ con ε_corr y q_RP3 (densa y poco densa).
   - Comienzos de ventana ℓ_ini ∈ {4.36, 4.66, 4.96, 5.16} y, si existe, ℓ_relax; final siempre en ℓ_fin − 0.05.
7. **Modelos de ξ(ℓ)** en [4.36, ℓ_fin], sobre la media de cada conjunto (RP³ poco densa ×6, RP³ densa, U(1) 576³ poco densa ×4, U(1) 576³ densa, U(1) 768³ con ξ_i = 0.5 ×4).
   - Se comparan una recta, ξ = a + bℓ, y una aproximación a una constante, ξ = ξ∞ − c/ℓ. Los dos modelos tienen dos parámetros.
   - Los errores por punto son la dispersión entre semillas dividida por √n; en los conjuntos de una o dos redes se usa la dispersión relativa de las seis RP³ poco densas.
   - Se informa Δχ² = χ²_sat − χ²_lin y el ξ∞ del modelo de saturación.
   - Lectura: |Δχ²| < 4 → INDISTINGUIBLE; Δχ² ≥ 4 → LINEAL_PREFERIDO; Δχ² ≤ −4 → SATURACION_PREFERIDA.
   - Se avisa de que con Δℓ ≈ 1.3 lo esperable es INDISTINGUIBLE, y de que un ξ∞ ajustado no es una predicción para ℓ ≈ 70.

**Prueba nueva propuesta (no preregistrada aquí):** U10 haría que una red RP³ densa se acerque de verdad al atractor desde arriba. Para eso, la red se prepararía con crecimiento de núcleo, de modo que no se desplome hacia ℓ ≈ 3, y después se pasaría a cuerdas físicas. Solo se preregistrará si las lecturas 2 y 6 lo justifican.

### Enmienda 1b (9 de octubre, ≈21:05 Madrid) · ajustes de `u9_extra.py` tras probarlo con la semilla densa 61

Al probar el script con los datos de la enmienda 1 (la misma semilla densa 61 y las redes de U6), la cota de conversión del punto 4 salió muy holgada. Atribuye a conversión el 40 % de la emisión incluso en las redes poco densas, que no tienen relajación inicial. Se añaden dos cosas, sin quitar nada:

- **Cota inferior de la conversión** con a⁴ρ_M: solo cuentan las pérdidas que no explica el corrimiento al rojo relativista. Se informan las dos cotas. ε_corr_lo y ε_corr acotan la emisión de las cuerdas por arriba y por abajo.
- **Fiabilidad de q en el barrido de ventana:** se marca «pocos puntos» cuando hay menos de tres estimaciones espectrales por semilla.

**Resultado de la prueba** (provisional, una sola semilla densa; no es el análisis final): la cota inferior es 0 y la superior está entre el 40 % (poco densas) y el 110 % (densa). Con los datos guardados (energía masiva integrada fuera de los núcleos), **no se puede separar** la conversión de la emisión de las cuerdas. Lo que sí se ve es que la red densa conserva un exceso de radiación masiva del 18 % al final de la ventana, es decir, que la relajación no termina dentro de ella. La separación fina necesitaría un diagnóstico nuevo en el kernel (el flujo de energía entre sectores o el espectro de la radiación masiva). Eso sería una prueba nueva, no un análisis de U9.
