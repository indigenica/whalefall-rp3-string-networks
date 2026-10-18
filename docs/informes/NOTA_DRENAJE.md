# Nota · Por qué es legítimo drenar el baño masivo inicial

*11 de octubre de 2026, noche. Tarea 7 de `docs/PLAN_NUBE_2026-10-11.md`.*
- *Comprobación numérica (informativa, sin veredicto):* `network/t7_drain_fade.py` → `network/results/u4/T7_DRAIN_FADE.json`.
- *Datos:* las redes de U6 y sus gemelas drenadas de U11 (576³, mismas semillas).

## 0. En una frase

**El drenaje quita un artefacto de la condición inicial de la red, no una parte de la producción cosmológica.** En el universo real, la radiación masiva que deje la transición decae a nivel árbol hacia ℓ ≈ 11–12, mucho antes de ℓ ≈ 70, y lo que produce son Goldstones duros, que cuentan como radiación oscura y no como materia oscura. En la propia red, la contribución del baño ya se está apagando dentro de la ventana.

## 1. Qué hace el drenaje y qué no hace

- **Hace:** hasta ℓ_d = 3.8, frena la parte masiva de la velocidad con un rozamiento de ritmo físico m_r, solo donde el campo está cerca del vacío (ESPEC_U10a, A). Desde ℓ_d la evolución es la física de siempre.
- **No hace:** no toca la radiación masiva que emiten las cuerdas después de ℓ_d, ni su conversión en Goldstones. Esa parte es física y sigue dentro de la R_Γ drenada. Por eso ρ_M/ρ_a de las redes drenadas **crece** a lo largo de la ventana (de 0.10 a 0.14 en RP³ y de 0.13 a 0.19 en U(1)): es la radiación masiva de la propia red, que se va acumulando.
- **Para qué sirve:** separar en el laboratorio dos cosas que el estimador de emisión mezcla. Una es la energía que las cuerdas ponen en Goldstones. La otra es la que pasa del baño masivo de la condición inicial (un campo aleatorio con longitud de correlación l_c, lejos del vacío) a los Goldstones, por interacción no lineal. La primera escala con la red; la segunda es un transitorio de la preparación.

## 2. Por qué en el universo real no llega a ℓ ≈ 70

**(a) Decae a nivel árbol mucho antes.** Las anchuras de los modos masivos de Q en el vacío (S3, `results/s3_sectors.json`, en las mismas unidades que la red) son:

| Modo | Γ/m | ℓ_dec = ln(m/Γ) |
|---|---|---|
| quintete (×5) | 6.0×10⁻⁶ | 12.0 |
| singlete | 2.4×10⁻⁵ | 10.6 |
| traza (pesado) | 2.5×10⁻⁵ | 10.6 |

El decaimiento ocurre cuando H ≈ Γ, es decir, hacia ℓ ≈ 11–12. En ℓ = 70, la fracción que sobrevive es ~exp(−Γ/2H) = exp(−3×10⁻⁶ · e⁷⁰), que es cero a todos los efectos. Un plasma solo acelera el decaimiento, y una población con mucha ocupación, como la de la red, convierte más deprisa todavía. Las anchuras son de vacío y a nivel árbol (lo dice `docs/modelo/MODELO_ACCION_COMUN.md`): sirven como cota del tiempo de vida, no como descripción de la cascada.

**(b) Lo que produce no es materia oscura.** Un modo de masa m_r decae en Goldstones con p ≈ m_r/2. En ℓ_dec ≈ 12 eso es x = p/H ≈ e¹²/2 ≈ 10⁵. Son Goldstones ultrarrelativistas que siguen siéndolo hasta hoy, porque m_G ≈ 10⁻¹⁵ eV. Cuentan como radiación oscura. La parte de la emisión que fija la materia oscura es la blanda: la de x ≈ x₀ = 10–30, emitida en las últimas e-foldings antes de que H ≈ m_G.

**(c) Su huella es ΔN_eff, que se acota aparte.** Si el baño lleva una fracción ε_b de la energía del plasma cuando decae (g* = 106.75), la radiación oscura resultante hoy da ΔN_eff ≈ 2.9 ε_b, tras la dilución por el recalentamiento del plasma del Modelo Estándar. Para ΔN_eff < 0.3 basta con ε_b ≲ 0.1. Este presupuesto depende de la historia térmica de la transición de Q, que la red no modela; es un control de radiación oscura, no de materia oscura.

## 3. En la propia red, el baño ya se está apagando

Si la contribución del baño fuera un transitorio, el cociente de emisión entre la red drenada y su gemela sin drenar debería subir hacia 1 a lo largo de la ventana. Eso es lo que se ve (ventana [4.36, 5.61], semillas de U11):

| | drenada/sin drenar (inicio → final) | pendiente por unidad de ℓ | ρ_M/ρ_a sin drenar | ρ_M/ρ_a drenada |
|---|---|---|---|---|
| RP³ (6 semillas) | 0.72 → 0.88 | **+0.28 ± 0.08** (3.7σ) | 1.11 → 0.90 | 0.10 → 0.14 |
| U(1) (4 semillas) | 0.77 → 0.79 | +0.04 ± 0.03 | 0.76 → 0.60 | 0.13 → 0.19 |

- **En RP³ el baño se apaga de forma significativa** dentro de la ventana. A ese ritmo, la diferencia entre drenada y sin drenar desaparecería en menos de una unidad de ℓ más allá del final de la ventana.
- **En U(1) la tendencia es la misma, pero no es significativa** con cuatro semillas.
- **Las dos poblaciones masivas van en sentidos opuestos:** la del baño (sin drenar) baja, la de la propia red (drenada) sube. Se acercan, como corresponde a un transitorio que se olvida.

Esto es un indicio, no una prueba: la ventana es corta y, en la red, la conversión es clásica y no lineal, con un ritmo que no es el de la ecuación de §2(a).

## 4. Qué queda abierto

- **El sesgo propio del rozamiento.** Puede afectar a la red, aunque U10a y U11 muestran que la densidad de cuerdas casi no cambia (R_ξ = 1.607 drenada frente a 1.617 sin drenar).
- **ℓ_d = 3.8 es una elección.** Con un ℓ_d más tardío se drenaría también radiación masiva de las cuerdas. Con uno más temprano, quedaría baño sin drenar. No se ha hecho un barrido de ℓ_d.
- **ε_b** (§2c) no está calculado: depende de cómo se ordena Q en la cosmología de la ballena (S3: «ordering switch» durante la inflación de bolsillo), y no de la red.

## 5. Consecuencia

Para el sector oscuro, la magnitud adecuada es la R_Γ drenada (U11: 1.56 ± 0.09). La R_Γ sin drenar (U6: 1.39) está sesgada a la baja, porque el baño inflaba más la emisión de U(1) (28 %) que la de RP³ (19 %). Esto no significa que el drenaje represente toda la producción cosmológica: en la cosmología, la energía del estado inicial acaba como radiación oscura dura y se cuenta en ΔN_eff, no en la materia oscura.
