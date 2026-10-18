# U10b · El atractor, acotado desde arriba y desde abajo

*11 de octubre de 2026.*
- *Preregistro:* `docs/preregistros/ESPEC_U10b.md` (≈10:40, antes de producir).
- *Resultados:* `network/results/u4/U10B_RESULT.json` (reproducido en la nube con el mismo script) y las redes `runs/*_ab.json`. *Figura:* `paper/zenodo_whale/figures/fig_bracket.pdf`.
- *Producción:* tres RP³ y tres U(1) a 576³, lanzadas tarde (ℓ_i = 3.5) y muy densas (ξ_i = 3 con dos semillas, ξ_i = 6 con una), drenadas hasta ℓ = 3.8. Unas 5.5 h en WSL.

## 0. En una frase

**Veredicto preregistrado: ATRACTOR_ACOTADO, en RP³ y en U(1).**
- Las redes que llegan desde arriba bajan durante toda la ventana, y las que llegan desde abajo (U11) suben. Dentro de la ventana no se encuentran.
- Al final de la ventana (ℓ = 5.61), el atractor queda **entre 1.12 y 1.73 para RP³** y **entre 0.67 y 1.10 para U(1)**.
- El crecimiento que medíamos desde U4b es compatible con un acercamiento desde abajo a una densidad constante, como sostienen Correia et al.
- **El cociente de densidades RP³/U(1) vale lo mismo desde los dos lados: R_ξ ∈ [1.57, 1.68].**

## 1. Densidad de cuerdas en la ventana [4.36, 5.61]

| | Desde arriba (U10b) | Desde abajo (U11) | Δ_fin = ξ_arriba/ξ_abajo |
|---|---|---|---|
| **RP³**, ξ al final | 1.73 (ξ_i = 3: 1.69; ξ_i = 6: 1.82) | 1.12 | 1.55 |
| RP³, dξ/dℓ | −0.20 (94 % de los tramos bajando) | +0.29 | |
| **U(1)**, ξ al final | 1.10 (ξ_i = 3: 1.08; ξ_i = 6: 1.16) | 0.67 | 1.66 |
| U(1), dξ/dℓ | −0.26 (88 % de los tramos bajando) | +0.14 | |

- **La puerta G_A** (que la prueba de verdad se acerque desde arriba) **pasa en los dos modelos.** Es lo que faltaba en U9.
- **Las redes desde arriba se aplanan al final.** En RP³, las de ξ_i = 3 bajan solo −0.12 por unidad de ℓ, y una sube un poco en la última instantánea. Las de ξ_i = 6 siguen bajando más deprisa (−0.37).
- **Informativo (extrapolación, no predicción):** las rectas ajustadas en las últimas 0.7 unidades de ℓ se cruzarían hacia ℓ ≈ 7.3, en ξ ≈ 1.6 (RP³) y ξ ≈ 0.9 (U(1)).
  - Ese ξ ≈ 0.9 para U(1) es del orden del ξ* ≈ 1.2 de Hindmarsh et al. Pero las definiciones de longitud no son idénticas: la nuestra es la longitud aparente, sin corregir por la velocidad.

## 2. El cociente RP³/U(1)

| | R_ξ | R_Γ | R_Γ/R_ξ |
|---|---|---|---|
| Desde abajo (U11, drenadas) | 1.61 | 1.56 ± 0.09 | 0.97 |
| Desde arriba (U10b) | 1.50 | 1.30 ± 0.14 | 0.87 |
| Acotamiento de R_ξ al final de la ventana | [1.57, 1.68] | | |

- **R_ξ ≈ 1.6** sale lo mismo con redes poco densas, densas, drenadas, sin drenar, desde abajo y desde arriba. Es el resultado más sólido de la red.
- **La emisión por unidad de longitud, comparada con U(1), es casi la misma con los dos acercamientos** (0.87 y 0.97). La red RP³ radia por metro de cuerda como una U(1), y radia más en total porque tiene más cuerda.
- **La pendiente espectral de las redes desde arriba no es comparable:** q(ℓ = 5) = 0.19, que sube a 0.70 al final de la ventana. Es la emisión de una red que se está deshaciendo: lazos colapsando y núcleos que se reajustan. No se usa.

## 3. Qué significa

**Para la red** (artículo de PRD):
- en la ventana que podemos simular, ℓ ≤ 5.7, **el atractor de escalado está acotado por los dos lados**, y el crecimiento de las redes poco densas es compatible con un acercamiento desde abajo;
- no demostramos que el atractor sea constante hasta ℓ ≈ 70. La ventana no distingue entre «constante» y «crece despacio»;
- U9 decía que la red RP³ crece igual empiece como empiece, pero ninguna de aquellas redes bajaba de verdad. U10b muestra que, si se empieza por encima del atractor después del desplome, la red baja. Las dos cosas son compatibles: en U9 todas las redes acababan por debajo del atractor después del desplome.

**Para el sector oscuro:**
- ATRACTOR_ACOTADO es la lectura preregistrada que **apoya la rama B2** (densidad constante), la que pasa Lyα en U4a, U7 y U7c. La rama B4 (ξ creciendo con nuestra pendiente) quedaba excluida en U5a; ahora tiene menos apoyo en nuestros datos.
- Con U11 (R_Γ limpia en la zona permitida y q < 1), **las dos condiciones de vida del sector oscuro son ahora compatibles con nuestras redes.** Nuestras redes ya no se pueden citar contra la ballena. Pero la condición se refiere a ℓ ≈ 70, y eso solo lo pueden decidir simulaciones mucho mayores.

## 4. Advertencias

- Las redes desde arriba nacen de un campo aleatorio proyectado en ℓ = 3.5, no de una historia física. Es el procedimiento estándar para explorar condiciones iniciales, como hacen Correia et al., pero no es la red cosmológica.
- Son tres redes por modelo. La de ξ_i = 6 todavía baja deprisa al final.
- La ventana dura 1.25 unidades de ℓ. El cruce extrapolado en ℓ ≈ 7.3 no es un resultado.
