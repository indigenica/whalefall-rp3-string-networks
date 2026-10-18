# U11 · U6 drenado: la R_Γ y la q limpias de la red RP³

*11 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U11.md` (10 de octubre, ≈19:50, antes de producir).
- *Resultados:* `network/results/u4/U11_RESULT.json` (reproducido en la nube con el mismo script) y las redes `runs/*_N576_*_dr.json`.
- *Producción:* 6 RP³ y 4 U(1) a 576³, con ξ_i = 0.5 y las semillas de U6, con el baño masivo inicial drenado hasta ℓ = 3.8. Unas 14 h en WSL. Esta vez N = 576 se mantuvo fijo, como estaba preregistrado.

## 0. En una frase

**Veredictos preregistrados: LIMPIA_ZONA_PERMITIDA y Q_MENOR_QUE_1.**
- Sin el baño de la condición inicial, la red RP³ radia **R_Γ = 1.56 ± 0.09** veces lo que la U(1); con la corrección de U(1) por el intercambio medido, 1.54 ± 0.09.
- Su espectro tiene pendiente **q = 0.87 ± 0.06** en ℓ = 5, que **baja con ℓ** (q ≈ 0.66 al final de la ventana).
- R_ξ = 1.61 se confirma una vez más.

## 1. Los números

| | U6 (sin drenar) | **U11 (drenado)** |
|---|---|---|
| R_Γ | 1.389 ± 0.062 | **1.562 ± 0.094** |
| R_Γ con U(1) corregida por X | — | **1.544 ± 0.087** |
| R_ξ | 1.617 ± 0.053 | **1.607** |
| q_RP3 (ℓ = 5) | 1.03 ± 0.12 | **0.873 ± 0.059** |
| dq/dℓ | −0.12 ± 0.29 | −0.32 ± 0.18 |
| q al final de la ventana (ℓ = 5.66) | 0.95 | 0.66 |
| ρ_M/ρ_a al final, RP³ / U(1) | 0.90 / 0.60 | 0.14 / 0.19 |

**Comparación semilla a semilla** (cada red drenada frente a su gemela de U6): la emisión que queda tras el drenaje es

| | Por semilla | Media |
|---|---|---|
| RP³ | 0.75–0.86 | **0.81 ± 0.05** |
| U(1) | 0.67–0.76 | **0.72 ± 0.04** |

Es decir, el baño de la condición inicial era el **19 % de la emisión de RP³** y el **28 % de la de U(1)** en las redes de U6. Como inflaba más a U(1), la R_Γ de U6 estaba sesgada a la baja: 1.39 en lugar de 1.56.

**El intercambio X en las U(1) drenadas es casi nulo** (−0.3 % a −2.3 % de la emisión). Una vez drenado el baño, lo que mide el estimador es emisión de las cuerdas, al menos en U(1).

**Coherencia con U10a** (512³, dos semillas): la R_Γ drenada con redes poco densas fue 1.63 ± 0.02, compatible. La q drenada fue 0.59. Ahora, con seis semillas, sale 0.87 ± 0.06. **La q sigue siendo el estimador más inestable** entre tamaños y conjuntos de semillas.

## 2. Veredictos

- **R_Γ: LIMPIA_ZONA_PERMITIDA.** R_Γ + 1σ = 1.66 ≤ 1.69, el valor con el que en U7c pasaban 7 de 8 filas de B2κ. La versión corregida da lo mismo.
- **q: Q_MENOR_QUE_1.** q + 2σ = 0.99 < 1. **Por muy poco:** 0.01 por debajo del umbral.

## 3. Qué significa para el sector oscuro

La condición de vida de U5 sigue siendo la misma: el sector oscuro solo vive si, a ℓ ≈ 70, ξ deja de crecer y q ≤ 1.025. Lo que cambia es dónde caemos dentro de ella.

| Si a ℓ ≈ 70… | Materia oscura | Lyα | Firma |
|---|---|---|---|
| q < 1 (como en nuestras redes limpias a ℓ ≈ 5) | casi toda de desalineamiento; la emisión aporta < 0.3 % (U5f con q = 0.9) | pasa sin problemas (δA ≈ 0) | m_G ≈ 2–10×10⁻¹⁵ eV, **sin componente tibia** |
| q ≈ 1 (Correia et al.) | 15–55 % emitida, con R_Γ ≈ 1.56 en la zona permitida | pasa, cerca del límite | supresión de pequeña escala δA ≈ 0.1–0.2 |
| q > 1.025 | demasiado tibia | **excluida** | — |

- **Nuestros propios datos ya no apuntan al lado excluido** en q, al contrario de lo que parecía en U5 (q ≈ 1.1–1.2, medido con el baño dentro).
- **La pendiente que medimos baja con ℓ**, pero eso no se puede extrapolar a ℓ ≈ 70. En U(1), unos grupos ven q creciendo con ℓ y otros (Correia et al.) ven q → 1.
- **El riesgo que queda para la materia oscura** es el de siempre: que la densidad de cuerdas siga creciendo (U5a, U9), o que q cruce 1 a ℓ grande.
- **El riesgo que queda para la predicción** es nuevo: si q < 1 se mantiene, la ballena pierde su firma tibia. Seguiría siendo falsable por la masa (m_G en el hueco de la superradiancia) y por n_s y r, pero no por Lyα.

## 4. Para el artículo de la red (PRD)

Los resultados que se pueden defender ya son:
1. **R_ξ = 1.6**, robusta frente a la densidad inicial (U9) y al drenaje (U10a, U11).
2. **Una sistemática del estimador de emisión.** El baño masivo de la condición inicial se cuenta como emisión: un 19 % (RP³) y un 28 % (U(1)) con redes poco densas, y más de la mitad con redes densas. Se corrige drenándolo, y en U(1) lo confirma una medida directa del intercambio.
3. **R_Γ limpia = 1.56 ± 0.09**, independiente de la densidad inicial dentro de ≈ 7 % (U10a).
4. **q_RP3 ≈ 0.87 ± 0.06 en ℓ = 5**, que baja con ℓ. Es el resultado menos robusto, porque el estimador depende del conjunto de semillas y del tamaño.
5. **La densidad de RP³ crece unas tres veces más deprisa que la de U(1)** en ℓ ≤ 6, con independencia de cómo empiece la red (U9).

## 5. Nota práctica

En el comando de lanzamiento, el primer `tail` no encuentra el log porque Python todavía no lo ha creado (tarda un segundo en arrancar). Y `tail -f` no termina nunca por diseño: sigue el fichero hasta que se corta con Ctrl+C, sin afectar a la producción. En adelante, el lanzamiento será `sleep 5; tail -f …`.
