# U10a · ¿De dónde sale el exceso de emisión de las redes densas? Drenaje de la radiación masiva y medida directa del intercambio · Preregistro

*9 de octubre de 2026, ≈23:58 (Madrid). Escrito antes de implementar nada de U10a y antes de cualquier ejecución.*

## Motivo

**Lo que dejó U9** (`docs/informes/INFORME_U9.md`). Las redes que empiezan densas (ξ_i = 2) emiten más que las poco densas: ×1.85 en RP³ y ×1.52 en U(1). Tienen además un espectro más azul (q ≈ 0.74 frente a 1.03). Por eso R_Γ sale 1.69 en lugar de 1.39. Al final de la ventana conservan un 20–29 % más de radiación masiva fuera de los núcleos.

**La hipótesis.** Parte de lo que el estimador cuenta como emisión de las cuerdas es radiación masiva de la condición inicial que se convierte en Goldstones. Con los datos guardados no se pudo separar.

**Por qué importa.** U7c mostró que R_Γ decide si el sector oscuro pasa Lyα: pasa con 1.39 y con 1.69, y falla con 2.09 y x₀ = 30.

U10a ataca la hipótesis de dos maneras independientes: una intervención (A) y una medida directa (B).

## A · Intervención: drenar la radiación masiva inicial

- **Qué se hace.** Hasta ℓ_d = 3.8, a la parte masiva de la velocidad se le aplica un rozamiento de ritmo físico γ = m_r (amortiguamiento crítico de los modos masivos).
  - La parte masiva es la componente normal al vacío: para RP³, la velocidad menos su proyección tangente w nᵀ + n wᵀ; para U(1), la componente radial.
  - El rozamiento solo actúa **fuera de los núcleos**, donde el parámetro de orden está cerca del vacío: hueco λ₁ − λ₂ > 0.6 en RP³ y |φ|/v > 0.6 en U(1). La dinámica de las cuerdas no se frena directamente.
  - Desde ℓ_d la evolución es la física de siempre.
- **Por qué ℓ_d = 3.8.** Las redes densas de U9 tocan fondo hacia ℓ ≈ 3.3–3.6, y la ventana de medida empieza en 4.36. Así se drena el baño que deja el desplome y queda medio e-fold de evolución libre antes de medir.
- **Lectura.** Si el exceso viene del baño inicial, las redes densas drenadas deben dar la misma R_Γ y la misma q que las poco densas drenadas.

## B · Medida directa del intercambio masivo → Goldstone

**La identidad.** Para U(1), la energía Goldstone es e_G = ½ r²(θ̇² + |∇θ|²), con r = |φ|. Cumple ∂_t e_G = X + ∇·(r² θ̇ ∇θ), donde

X = 2 (ṙ/r) (e_G,grad − e_G,cin)

es la transferencia local desde el modo radial (masivo) al sector Goldstone. El término de divergencia es el flujo de energía Goldstone. Integrado sobre la región fuera de los núcleos, es lo que sale de las cuerdas, es decir, su emisión.

**Para RP³ se usa el análogo plano a plano**, en la base propia de Q:
- e_G = Σ_k [u_k² + Σ_i g_ki²], con u_k = e_kᵀ Q̇ e₁ y g_ki = e_kᵀ ∂_i Q e₁, para k = 2…4;
- X = Σ_k 2 (λ̇₁ − λ̇_k)/(λ₁ − λ_k) · (Σ_i g_ki² − u_k²), con λ̇_k = e_kᵀ Q̇ e_k.

Es exacta en el límite abeliano. Omite los acoplamientos no abelianos (las rotaciones entre e₂…e₄ y la curvatura de la conexión). **Su validez se comprueba en V5**, no se supone.

**Muestreo.**
- Después de cada instantánea, se evalúan ⟨X⟩, ⟨e_G,cin⟩ y ⟨e_G,grad⟩ cada 4 pasos durante 32 pasos.
- La región es la de las celdas a más de (4/m + 8 dτ_fís) de cualquier plaqueta perforada en la instantánea; así se cubre lo que se mueven las cuerdas en ese intervalo.
- ⟨X⟩ es la media de las muestras válidas, con fracción ≥ 0.5.
- **Emisión limpia:** Γ_cuerdas = Γ − ⟨X⟩. Γ es la del estimador de siempre (rm2.5), y ⟨X⟩ se interpola al mismo ℓ.

**También se guarda el espectro de la velocidad masiva** fuera de los núcleos en cada instantánea, como diagnóstico informativo.

## V5 · Validación de X (antes de la producción, en la nube)

- **Montaje.** Caja sin cuerdas en Minkowski (a = 1), N = 64, dτ = 0.1, 600 pasos. Vacío más perturbaciones gaussianas pequeñas en todas las componentes: amplitud 0.15 y l_c = 4 celdas, con el vacío de RP³ en n = (1, 0, 0, 0). Se comprueba que el número de plaquetas perforadas es 0.
- **Medida.** Cada 2 pasos se mide la energía Goldstone total E_G(t) (cinética + gradiente) y la transferencia ∫X dt.
- **Criterio, para cada modelo:** max_t |ΔE_G(t) − ∫₀ᵗ X| ≤ 0.15 · max_t |ΔE_G(t)|.
- **Si falla en un modelo**, X queda como diagnóstico cualitativo en ese modelo, y el veredicto de U10a se basa solo en A.

## Producción

`u4b_produce.py 6` en la máquina WSL, a 576³, con espectros en cada instantánea (ESPEC_U9) y X en todas las redes. En este orden:

1. RP³ drenada, ξ_i = 2, semilla 61.
2. U(1) drenadas, ξ_i = 2 y ξ_i = 0.5, semilla 61 (en pareja).
3. RP³ drenada, ξ_i = 0.5, semilla 61.
4. RP³ **sin drenar**, ξ_i = 2, semilla 61, con X. Repite U9 s61 bit a bit, más X; es una comprobación de reproducibilidad.
5. U(1) sin drenar, ξ_i = 2, semilla 61, con X.
6. RP³ drenada, ξ_i = 2 y ξ_i = 0.5, semilla 62.
7. U(1) drenadas, semilla 62 (en pareja).

**Tiempo estimado:** ≈ 2 h por red RP³ y ≈ 0.7 h por red U(1), unas 13 h en total. El análisis se puede hacer con la semilla 61 (pasos 1–5) como **lectura provisional**. El veredicto exige las dos semillas drenadas.

## Análisis (`network/u10a_analyze.py`, fijado antes de los datos)

**Ventana y medidas.**
- La ventana de R_Γ es la de U4b, [ℓ_fin − 1.3, ℓ_fin − 0.05], con el mismo estimador que U6 y U9.
- R_Γ^dren,denso = RP³ drenada densa / U(1) drenada densa; R_Γ^dren,poco = lo mismo con las poco densas.
- q_RP3 como en U6, para cada protocolo.
- Exceso masivo: (ρ_M/ρ_a)_denso / (ρ_M/ρ_a)_poco − 1 al final de la ventana, en las drenadas.
- Brecha de protocolo: G = |R_Γ^dren,denso − R_Γ^dren,poco| / media. En U9, sin drenar, G = 0.195.

**Veredicto A:**
- **RELAJACION_EXPLICA:** G < 0.08, |q_denso − q_poco| < 0.15 y exceso masivo < 10 %. La R_Γ limpia es la media de las dos R_Γ drenadas.
- **RELAJACION_NO_EXPLICA:** exceso masivo < 10 % y G ≥ 0.15. El drenaje quita el baño, pero la brecha sigue.
- **INDETERMINADO:** cualquier otro caso, incluido que el drenaje no quite el exceso masivo.

**B, informativo** (solo si V5 pasa en ese modelo):
- f_X = ⟨X⟩/Γ en la ventana, para cada red;
- R_Γ^limpia = cociente de (Γ − X) para cada protocolo;
- comprobación cruzada: (Γ − X) de la RP³ densa sin drenar frente a Γ de la RP³ densa drenada.

**Control:** la RP³ sin drenar con X debe reproducir U9 s61 en ξ(ℓ) con un error relativo ≤ 10⁻⁵.

## Lo que U10a no resuelve

- El comportamiento a ℓ ≈ 70.
- Si la radiación masiva que **emiten las cuerdas** (no la de la condición inicial) debe contarse en la emisión. En el universo real esa radiación acaba en Goldstones, así que es emisión legítima. El drenaje solo actúa antes de ℓ_d y no la toca dentro de la ventana.
- U10b (aproximación desde arriba) sigue pendiente.

---

## Enmienda 1 (10 de octubre, ≈01:00 Madrid) · resultado de V5 y ajustes de implementación, antes de la producción

**V5 tal como estaba preregistrada** (`network/validate_u10a.py`, resultado en `network/results/u4/U10A_V5.json`):

| Modelo | Error relativo de ∫X | Plaquetas perforadas | Veredicto |
|---|---|---|---|
| U(1) | 0.123 | 0 → 0 | **PASA** |
| RP³ | 0.376 | 0 → 230 | **FALLA** |

En RP³, con amplitud 0.15 se forman cuerdas durante la evolución, así que la prueba no cumple su propia condición de caja sin cuerdas.

**Pruebas informativas** (no preregistradas, no se usan), RP³ sin cuerdas:
- con amplitud 0.08, error 0.34;
- con amplitud 0.05, error 0.47.

Con amplitudes pequeñas el error crece porque domina la discretización. X usa diferencias centradas, mientras la dinámica usa la acción de la red con el laplaciano de 7 puntos. Para U(1) existe una identidad exacta en la red, con energías de enlace, pero no se adopta aquí.

**Consecuencia**, como estaba preregistrada: **en RP³, X queda como diagnóstico cualitativo.** El veredicto de U10a se basa solo en A, el drenaje. En U(1), X es cuantitativo.

**Comprobación del drenaje** (también en `U10A_V5.json`): con factor 0, la velocidad masiva pasa de rms 0.078 a 0.006 en RP³ (las celdas cerca de núcleos no se tocan) y de 0.096 a 2×10⁻⁹ en U(1). La velocidad tangente cambia menos de 10⁻⁷.

**Ajustes de implementación** (no cambian el diseño):
1. **X solo donde sirve.** Como en RP³ X ya no es cuantitativo, no se muestrea en las redes RP³ drenadas: ahorra ≈ 40 % de tiempo. Se mantiene en todas las U(1) y en la RP³ de control sin drenar (paso 4).
2. **El drenaje se aplica cada 4 pasos** con el factor acumulado, exp(−γ a 4 dτ). Es el mismo ritmo, y el índice del paso se calcula a partir de τ, así que no depende de cómo se trocee la evolución.
3. **Radio de la región de X**, aclarado: 4/(m a) + 0.5 · 32 · dτ celdas (cuatro anchuras de núcleo más lo que recorre una cuerda a v ≤ 0.5 durante el muestreo). Sustituye a la expresión ambigua «4/m + 8 dτ_fís».
4. **Exceso masivo del veredicto A:** el máximo en valor absoluto entre RP³ y U(1). Queda fijado en `network/u10a_analyze.py` antes de los datos.
5. **La evolución sin drenar y con X es idéntica bit a bit** a la de siempre: se comprobó con N = 64 y con N = 96. `make test` pasa, incluida la reproducción de una red de U2.

**Exploración previa** (N = 128, τ_fin = 64, una semilla, ξ_i = 2; no cuenta para nada). Con el drenaje hasta justo después del desplome (ℓ_d = 3.6):
- ρ_M/ρ_a al final baja de 1.5 a 0.18 en RP³ y de 0.99 a 0.19 en U(1);
- la ρ_a de RP³ drenada (1.01×10⁻³) se acerca a la de la poco densa sin drenar (1.00×10⁻³), frente a 1.94×10⁻³ sin drenar.

Con ℓ_d = 2.8 (antes del desplome) el drenaje apenas tiene efecto, porque el baño se produce en el desplome. Esto confirma la elección de ℓ_d = 3.8 en la producción.

**Tiempo estimado de la producción:** ≈ 2 h por RP³ drenada, ≈ 2.4 h para la RP³ de control con X y ≈ 0.5–0.7 h por U(1). En total, 12–13 h.

---

## Enmienda 2 (10 de octubre, ≈19:45 Madrid) · desviaciones de la producción, escrita DESPUÉS de ver los datos

**Tamaño de la caja.** `u4b_produce.py 6` tomó N del planificador, que lo calcula según la RAM libre. Al lanzar había 41.9 GB libres y eligió **N = 512**, no 576. Consecuencias:
1. **Escala de la ventana:** ℓ_fin = 5.545 y la ventana de R_Γ es [4.25, 5.50]; con 576 habría sido [4.36, 5.61].
2. **El veredicto A no cambia de diseño**, porque todas sus redes están a 512³.
3. **El control** «la RP³ sin drenar reproduce U9 s61 bit a bit» **no se puede hacer** (otro N); se informa como «no aplicable».
4. **Comparaciones con las redes sin drenar poco densas:** solo existen a 576³ (U6 y U9). Se hacen a la misma ventana absoluta de ℓ, que da R_Γ = 1.43 (U6) y 1.67 (U9), frente a 1.39 y 1.69 en su ventana propia. Son informativas y entre tamaños distintos.

**Números adicionales calculados después de los datos:** están en `network/u10a_extra.py` → `U10A_EXTRA.json`. No forman parte del veredicto.
