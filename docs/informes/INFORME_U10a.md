# U10a · Drenar la radiación masiva de la condición inicial: cuánto de la «emisión» era de las cuerdas

*10 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U10a.md`, con la enmienda 1 (V5 y ajustes, antes de la producción) y la enmienda 2 (desviaciones, después de los datos).
- *Resultados:* `network/results/u4/U10A_RESULT.json` (veredicto preregistrado) y `U10A_EXTRA.json` (números añadidos después de ver los datos). La validación está en `U10A_V5.json`.
- *Producción:* diez redes a **512³** (el planificador eligió 512 por la RAM libre; enmienda 2). Son ocho drenadas (RP³ y U(1), ξ_i = 2 y 0.5, semillas 61 y 62) y dos sin drenar con X (ξ_i = 2, semilla 61). Unas 13 h en WSL.

## 0. En una frase

**Más de la mitad de lo que medíamos como emisión de las redes densas era radiación masiva de la condición inicial convirtiéndose en Goldstones.**
- **Lo muestran dos métodos independientes, que coinciden en U(1):** el drenaje deja el 43 % de la emisión, y la medida directa del intercambio dice que el 45 % es de las cuerdas.
- **Las redes poco densas también llevaban baño:** el drenaje les quita entre un 23 % (RP³) y un 32 % (U(1)) de la emisión.
- **Sin el baño, los dos protocolos casi coinciden:** R_Γ = 1.76 ± 0.09 (densas) y 1.63 ± 0.02 (poco densas). La brecha baja de 0.195 a 0.074, y **la R_Γ limpia queda en ≈ 1.6–1.8**.
- **La pendiente espectral baja** a q ≈ 0.6–0.7.
- **El veredicto preregistrado es INDETERMINADO,** porque una de sus tres condiciones (el exceso masivo relativo) no se cumple (§2).

## 1. Qué hace el drenaje

Radiación masiva fuera de los núcleos, ρ_M/ρ_a, en ℓ = 4.3 / 4.8 / 5.3:

| | RP³ | U(1) |
|---|---|---|
| Densa sin drenar (512³) | 1.43 / 1.23 / 1.09 | 0.98 / 0.90 / 0.80 |
| Densa drenada | 0.19 / 0.24 / 0.26 | 0.20 / 0.25 / 0.29 |
| Poco densa sin drenar (U6, 576³) | 1.13 / 0.99 / 0.92 | 0.78 / 0.66 / 0.63 |
| Poco densa drenada | 0.10 / 0.12 / 0.14 | 0.13 / 0.15 / 0.18 |

El drenaje quita el 80–90 % del baño. Lo que queda va creciendo poco a poco: es la radiación masiva que emiten las propias cuerdas, que es física.

## 2. Veredicto A, preregistrado

| Condición | Valor | Umbral | ¿Cumple? |
|---|---|---|---|
| Brecha de protocolo G | **0.074** (U9 sin drenar: 0.195) | < 0.08 | sí |
| \|q_denso − q_poco\| | **0.107** | < 0.15 | sí |
| Exceso masivo relativo E_M (máximo entre modelos) | **0.87** (RP³: 0.27 frente a 0.14; U(1): 0.29 frente a 0.19) | < 0.10 | **no** |

**Veredicto: INDETERMINADO.** Las dos condiciones sobre R_Γ y q se cumplen. Falla la tercera: tras el drenaje, la red densa conserva el doble de radiación masiva **relativa** que la poco densa. En términos absolutos es de cinco a diez veces menos que sin drenar. El criterio era relativo y se informa tal cual.

## 3. Números informativos (después de los datos; `U10A_EXTRA.json`)

**Emisión que queda tras el drenaje:**

| | RP³ | U(1) |
|---|---|---|
| Densa (frente a la misma red sin drenar, 512³) | **0.38** | **0.43** |
| Poco densa (frente a U6, 576³, misma ventana de ℓ) | 0.77 | 0.68 |

**La medida directa del intercambio, X:**
- **U(1), cuantitativa** (V5 pasó):
  - en la densa sin drenar, X es el **55 %** de la emisión, así que las cuerdas aportan el 45 %; **el drenaje dice 43 %**. Dos métodos independientes coinciden;
  - en las drenadas, X es negativo: −0.3 en las densas y ≈ −0.03 en las poco densas. Tras el drenaje es el sector Goldstone el que cede algo de energía al masivo, que vuelve a llenarse. Con esa corrección, la emisión de las cuerdas en U(1) densa drenada sería un ≈ 30 % mayor que la medida.
- **RP³, cualitativa** (V5 falló): en la densa sin drenar, X ≈ 58 % de la emisión. El drenaje deja el 38 %. Es coherente en orden de magnitud.

**R_Γ por protocolo** (ventana [4.25, 5.50]):

| | R_Γ | R_ξ |
|---|---|---|
| Drenadas, densas (2 semillas: 1.89 y 1.62) | 1.76 ± 0.09 | 1.39 |
| Drenadas, poco densas (2 semillas: 1.67 y 1.59) | 1.63 ± 0.02 | 1.59 |
| Densa sin drenar, 512³ (1 semilla) | 2.03 | 1.34 |
| U6, poco densas sin drenar (576³, misma ventana de ℓ) | 1.43 ± 0.06 | 1.59 |
| U9, densas sin drenar (576³, misma ventana de ℓ) | 1.67 ± 0.03 | 1.56 |

- **Las R_Γ drenadas caen entre 1.59 y 1.89 en las cuatro parejas.** El baño inflaba más la emisión de U(1) que la de RP³ en las redes poco densas, por eso la R_Γ de U6 (1.39) era **baja**.
- **Con la corrección por X de las U(1) drenadas**, solo disponible en U(1), la R_Γ densa drenada bajaría hacia ≈ 1.35 y la poco densa apenas se movería (≈ 1.58). No sabemos si RP³ tiene una corrección parecida (allí X es cualitativa).
- **Estimación limpia, con su incertidumbre: R_Γ ≈ 1.6 ± 0.2.**

**Pendiente espectral q_RP3 (ℓ = 5):**

| | q |
|---|---|
| Drenadas, densas | 0.70 |
| Drenadas, poco densas | 0.59 |
| U6 sin drenar | 1.03 ± 0.12 |
| U9 sin drenar | 0.74 ± 0.07 |

Sin el baño, el espectro de las cuerdas sale **más azul**: tiene más energía en las ondas cortas. Los errores bootstrap de las drenadas (≤ 0.01, con dos semillas) no recogen las sistemáticas. Lo único firme es que **q depende de cómo se trate el baño**, entre 0.6 y 1.0, y que el estimador de q a ℓ ≈ 5 es frágil.

## 4. Qué significa

**Para la medida de la red:**
- El estimador de emisión de GHV (radiación libre = 2 × energía cinética tangente fuera de los núcleos) **cuenta como emisión la conversión del baño masivo de la condición inicial**. En nuestras redes es el 25–60 % de la señal, según el protocolo.
- No es la contaminación que describen Correia et al. (el campo propio de las cuerdas, que es otro efecto). Es una sistemática nueva y medible. Afecta sobre todo a RP³, que tiene siete modos masivos frente a uno de U(1).
- Con el baño drenado, la dependencia del protocolo casi desaparece. **R_ξ sigue siendo robusta (≈ 1.4–1.6)** y **R_Γ limpia ≈ 1.6–1.8**.

**Para el sector oscuro:**
- **Si a ℓ ≈ 70 el espectro tiene q ≈ 1**, la R_Γ limpia (≈ 1.7) está donde U7c dejaba pasar casi todas las filas (7/8 en B2κ). La predicción tibia sobrevive, al borde del límite.
- **Si q < 1**, como apuntan las redes drenadas, la emisión deja de contar (U5f). La materia oscura sería fría, de desalineamiento: pasa Lyα sin problemas, pero **pierde la firma tibia**.
- **En ningún caso la limpieza mata al sector oscuro.** Lo que hace es mover la pregunta de vuelta a q.

**Para el artículo:** la sistemática del baño masivo, y cómo se corrige drenándolo, es un resultado metodológico publicable por sí solo. Las simulaciones de cuerdas globales con varios modos masivos necesitan quitar o medir ese baño antes de citar una eficiencia de emisión.

## 5. Siguiente

1. **Repetir U6 drenado a 576³** (6 + 4 semillas). Daría la R_Γ y la q limpias con la estadística de U6.
2. **Una X exacta en la red para RP³**, con energías de enlace como la identidad de U(1). Permitiría corregir también las redes drenadas.
3. **U10b** (aproximación desde arriba), pendiente.
