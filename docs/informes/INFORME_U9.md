# U9 · Redes RP³ que empiezan densas: qué es robusto y qué depende de cómo empieza la red

*9 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U9.md`, con las enmiendas 1 y 1b, que añaden análisis informativos a raíz de una revisión externa.
- *Resultados:* `network/results/u4/U9_RESULT.json` (clasificación preregistrada) y `U9_EXTRA.json` (enmienda 1). Las redes están en `runs/*_xi2.json` y el log en `produce_phase5.log`.
- *Producción:* dos RP³ y dos U(1) a 576³ con ξ_i = 2 (semillas 61 y 62), frente a las seis RP³ y cuatro U(1) con ξ_i = 0.5 de U6.

## 0. En una frase

**La densidad de cuerdas de RP³ es robusta y la emisión no.**
- Empiece densa o poco densa, la red RP³ acaba con la misma densidad (memoria del 4.6 %) y crece al mismo ritmo dentro de la ventana simulada (dξ/dℓ = 0.310 frente a 0.313). El cociente de densidad con U(1) es el mismo: R_ξ = 1.58 frente a 1.62.
- En cambio, el cociente de emisión y la pendiente espectral dependen de cómo empieza la red: R_Γ = 1.69 frente a 1.39, y q = 0.74 frente a 1.03.

**Veredicto preregistrado: CRECIMIENTO_GENUINO.** Según la enmienda 1, significa: *crecimiento persistente dentro de la ventana numérica estudiada (ℓ ≤ 5.66), independiente de la densidad inicial de cuerdas ensayada*. No demuestra que la red siga creciendo hasta ℓ ≈ 70, ni refuta la lectura de Correia et al. a ℓ grande (§2).

## 1. Densidad de cuerdas: universal en la ventana

| | RP³ ξ_i = 0.5 (6 semillas) | RP³ ξ_i = 2, semilla 61 | RP³ ξ_i = 2, semilla 62 |
|---|---|---|---|
| ξ(ℓ_ref = 5.36) | 0.990 ± 0.055 | 1.078 | 0.995 |
| dξ/dℓ en [4.66, 5.36] | 0.313 ± 0.046 (0.26–0.37) | 0.312 (−0.03 σ) | 0.308 (−0.12 σ) |

- **Memoria de la condición inicial en ℓ_ref: 4.6 %** para RP³ y 9.6 % para U(1). Las dos están por debajo del límite preregistrado del 20 %.
- **Las redes densas caen dentro de la distribución de las poco densas**, en densidad y en pendiente.
- **R_ξ(denso) = 1.58 ± 0.09, frente a 1.62 ± 0.05 en U6.** El cociente de densidad entre RP³ y U(1) no depende de cómo empieza la red. Es el resultado más sólido de la red hasta ahora.
- **Recta frente a saturación** (punto 7 de la enmienda): INDISTINGUIBLE en los cinco conjuntos (|Δχ²| ≤ 3.2), como se había anunciado. Los ajustes de saturación dan ξ∞ ≈ 2.3–2.4 para RP³ y ≈ 1.05–1.22 para U(1). Este último está, por casualidad o no, en el valor de Hindmarsh et al. (1.19). **No son predicciones para ℓ ≈ 70**: con Δℓ ≈ 1.3 la ventana no distingue entre una recta y una curva que se aplana.

## 2. ¿Se acercó alguna red «desde arriba»? No

Las cuatro redes densas (dos RP³ y dos U(1)) se desploman en ℓ ≈ 3.3–4.3, cuando el horizonte es solo unas 25–75 veces más grueso que la cuerda. Quedan por debajo de su valor final y, dentro de la ventana, ξ **sube en los 8 tramos de las cuatro**.

Por tanto, U9 no ensayó lo que describen Correia et al.: una red que se acerca al régimen de escalado bajando desde arriba. Lo que muestra es más limitado: **el crecimiento en ℓ ≤ 6 no es memoria de haber empezado poco denso**, porque una red que empieza muy densa, tras su desplome, sigue la misma trayectoria. Eso debilita una interpretación concreta de su argumento («nuestras redes crecen porque empezaron poco densas»), pero no lo refuta a ℓ grande.

## 3. Emisión: depende de cómo empieza la red

| | ξ_i = 0.5 (U6) | ξ_i = 2 (U9) |
|---|---|---|
| R_Γ en la ventana preregistrada [4.36, 5.61] | **1.39 ± 0.06** | **1.69** (bootstrap ± 0.03 con dos semillas por modelo, poco fiable; por ℓ va de 1.0 a 2.9) |
| R_Γ empezando la ventana en 4.66 / 4.96 / 5.16 | 1.41 / 1.29 / 1.22 | 1.50 / 1.74 / 2.09 |
| q_RP3 (ℓ = 5) | **1.03 ± 0.12** | **0.74 ± 0.07** |
| q_RP3 empezando en 4.66 / 4.96 | 0.97 / 0.95 | 0.69 / 0.73 |
| Emisión densa / poco densa | — | RP³: 1.85; U(1): 1.52 |

- **El exceso de emisión de las redes densas no se cancela en el cociente.** RP³ gana un 85 % y U(1) un 52 %. Por eso R_Γ sube de 1.39 a 1.69.
- **Al retrasar la ventana, los dos protocolos no convergen.** R_Γ baja con las redes poco densas y sube con las densas, aunque en las ventanas tardías hay pocos puntos.
- **La pendiente es más baja con las redes densas:** sobra energía en las ondas cortas.
- **Qué lo causa** (puntos 3–5 de la enmienda): las redes densas cargan más radiación masiva fuera de los núcleos. Al final de la ventana tienen todavía un 20 % (RP³) y un 29 % (U(1)) más que las poco densas, así que **la relajación no termina dentro de la ventana**.
  - La lectura más simple: parte de lo que el estimador llama «emisión de las cuerdas» es radiación masiva de la condición inicial que se va convirtiendo en Goldstones a escalas cercanas a la masa. Eso explicaría a la vez el exceso de Γ y el q más bajo.
  - **Con los datos guardados no se puede separar.** La cota de conversión queda entre el 0 y el 113 % de la emisión (enmienda 1b), y las emisiones «corregidas» no tienen sentido.
- **Exploratorio, post hoc, no usado:** con los espectros instantáneos guardados, el exceso de RP³ denso es menor (×1.1–1.2) si solo se cuentan las ondas por debajo de m/4. Es compatible con la lectura anterior, pero la emisión por bandas sale ruidosa e incoherente entre bandas con dos semillas, y no se le da peso.

## 4. Qué cambia

**La medida de la red** (eslabón 8; artículo §4):
- **Robusto:** la red RP³ es ≈ 1.6 veces más densa que la U(1), con independencia de la densidad inicial ensayada.
- **Robusto:** la red RP³ radia más que la U(1) con los dos protocolos (R_Γ > 1). No hay red silenciosa.
- **Robusto en la ventana:** la densidad de RP³ crece unas tres veces más deprisa que la de U(1).
- **Con sistemática de protocolo:** R_Γ ≈ 1.4–1.7 en la ventana preregistrada y 1.2–2.1 según dónde empiece. q ≈ 0.7–1.0. Ya no se puede publicar un único R_Γ = 1.39 ± 0.06 como «la» emisión de la red en escalado.

**El sector oscuro** (eslabón 10):
- U4a pasó con R_Γ hasta 1.98, con la emisión sin calibrar. U7 (calibrada) se calculó solo hasta R_Γ = 1.53. **No sabemos si las filas calibradas pasan con R_Γ ≈ 1.7–2.1.** Lo decide **U7c** (`docs/preregistros/ESPEC_U7c.md`, preregistrado antes de calcular: 24 filas de CLASS, 1–2 h en WSL).
- **Las dos sistemáticas empujan en sentidos opuestos.** Más R_Γ hace la materia oscura más tibia. Menos q (más energía en las ondas cortas) hace la emisión irrelevante: con q < 1, la materia oscura es fría de desalineamiento y pasa Lyα **sin la firma tibia** (U5f).
- La condición de vida de U5 (ξ constante y q ≤ 1.025 a ℓ ≈ 70) sigue siendo la misma, y sigue fuera del alcance de cualquier red.

**Los tres escenarios** (de la revisión externa):
- **A, la red crece en nuestras cajas y luego se estabiliza:** B2 sigue siendo candidata. U9 no lo excluye.
- **B, sigue creciendo hasta ℓ ≈ 70:** la rama B2 queda excluida por Lyα (U5a). La medida de la red sigue siendo de interés por sí misma.
- **C, la emisión depende de cómo empieza la red:** **U9 muestra que, a 576³, es así.** Antes de dar una calibración definitiva hay que estudiar la relajación y reforzar el estimador.

## 5. Pruebas siguientes posibles

1. **U7c** (preregistrado): Lyα con R_Γ = 1.39, 1.69 y 2.09.
2. **U10a:** un diagnóstico nuevo en el kernel que mida directamente el flujo de energía entre la radiación masiva y la Goldstone, o el espectro de la masiva, para separar la relajación de la emisión de las cuerdas. Supone tocar el código C y repetir una pareja de cada modelo y protocolo (≈ 5 h).
3. **U10b:** una red densa que se acerque de verdad desde arriba. Se prepararía con crecimiento de núcleo para que no se desplome hacia ℓ ≈ 3, y después se pasaría a cuerdas físicas. Es la prueba directa del argumento de Correia et al.
4. **Para el artículo sobre la red (PRD):** la pregunta interesante que deja U9 es *por qué la red RP³ es 1.6 veces más densa y crece unas tres veces más deprisa que la U(1)*, con independencia de cómo empiece. ¿Tiene que ver con los tres Goldstones y la topología ℤ₂? Es una pregunta de física de campos que tiene sentido aunque la materia oscura de la ballena no sobreviva.
