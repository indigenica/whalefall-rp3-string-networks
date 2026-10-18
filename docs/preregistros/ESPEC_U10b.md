# U10b · Redes que se acercan al atractor desde arriba · Preregistro

*11 de octubre de 2026, ≈10:40 (Madrid). Escrito antes de producir ninguna red de U10b.*

## Motivo

- **U9** mostró que la densidad de cuerdas de RP³ crece igual dentro de la ventana empiece la red densa o poco densa. Pero ninguna red se acercó desde arriba: las densas se desploman hacia ℓ ≈ 3.3–4.3, cuando el horizonte es solo 25–75 veces más grueso que la cuerda, y después suben.
- **El argumento de Correia et al.** es que el crecimiento que ven otros grupos es un acercamiento desde abajo a una densidad constante ξ*.
- **U10b lo pone a prueba directamente.** Lanza redes muy densas **tarde**, en ℓ_i = 3.5, cuando el desplome genérico ya ha pasado, para que lleguen al atractor bajando. Se comparan con las de U11, que llegan subiendo.

## Exploración previa (solo para elegir el diseño; no cuenta)

En la nube, con redes pequeñas y una semilla:

**U(1), 192³, ℓ_fin = 4.56:**
- arrancando en ℓ_i = 3.0 con ξ_i = 2 o 4, o en ℓ_i = 3.5 con ξ_i = 3, la densidad baja de forma monótona: en ℓ = 4.5 vale 1.00, 1.11 y 1.39 respectivamente;
- la red poco densa estándar sube de 0.38 a 0.46.

**RP³, 128³, ℓ_fin = 4.16:**
- arrancando en ℓ_i = 3.2 con ξ_i = 3, la densidad baja de 2.5 a 1.66;
- la red poco densa sube de 0.46 a 0.58.

En las dos, las redes «desde arriba» bajan despacio y la distancia con las poco densas no se cierra en esas ventanas cortas.

## Producción

`u4b_produce.py 8`, con **N = 576 fijo** y la misma parada por falta de RAM que en U11.
- Arranque en **ℓ_i = 3.5**.
- Drenaje de la radiación masiva hasta ℓ = 3.8, como en U11, para quitar el baño del arranque.
- Espectros en cada instantánea; X en las U(1).
- Etiqueta `_ab`.

| Red | ξ_i | Semillas |
|---|---|---|
| RP³ | 3 | 71, 72 |
| RP³ | 6 | 73 |
| U(1) | 3 | 71, 72 |
| U(1) | 6 | 73 |

**Orden:** RP³ ξ3 s71; U(1) ξ3 s71 y s72 en pareja; RP³ ξ3 s72; U(1) ξ6 s73; RP³ ξ6 s73.

**Tiempo estimado:** ≈ 1.5 h por RP³ (empiezan más tarde) y ≈ 0.4 h por U(1); unas 5–6 h en total.

**Referencia «desde abajo»:** las redes drenadas de U11 (RP³ 61–66 y U(1) 61–64, ξ_i = 0.5, 576³).

## Análisis (`network/u10b_analyze.py`, fijado antes de los datos)

La ventana es W = [4.36, ℓ_fin − 0.05], con ℓ_fin = 5.66, y se usan las medias por conjunto. Para cada modelo (RP³ y U(1)), con todas las redes `_ab` juntas:
- **s_ab y s_low:** las pendientes dξ/dℓ en W de las desde arriba y de las desde abajo (U11).
- **Δ_fin = ξ_ab / ξ_low** en ℓ_fin − 0.05.
- **Puerta G_A**, la de que la prueba sea válida: s_ab < 0 y más de la mitad de los tramos de W bajando en la media de las `_ab`.

**Veredicto por modelo:**
- **CONVERGENCIA:** Δ_fin < 1.15. Desde arriba y desde abajo llegan al mismo valor dentro de la ventana; se informa ese valor y su pendiente.
- **ATRACTOR_CRECIENTE:** no hay convergencia y s_ab ≥ 0. Las redes desde arriba también crecen: el crecimiento es del atractor.
- **ATRACTOR_ACOTADO:** no hay convergencia, G_A pasa y s_low > 0. En ℓ_fin, el atractor queda entre ξ_low y ξ_ab. El crecimiento de las desde abajo es compatible con un acercamiento desde abajo, como dicen Correia et al.
- **INDETERMINADO:** cualquier otro caso.

**Informativo:**
- dónde y en qué valor se cruzarían las dos rectas ajustadas en las últimas 0.7 unidades de ℓ. Es una extrapolación; no es una predicción;
- el acotamiento de la densidad relativa: el cociente RP³/U(1) del límite inferior y del superior;
- R_Γ, R_ξ y q_RP3 de las redes `_ab` frente a U11, y la emisión por unidad de longitud, R_Γ/R_ξ, en los dos casos;
- los resultados por densidad inicial (ξ_i = 3 y 6).

## Lectura para el sector oscuro

- **ATRACTOR_ACOTADO o CONVERGENCIA con ξ casi constante:** apoyan la rama B2 (ξ constante, U4a y U7).
- **ATRACTOR_CRECIENTE:** apunta a B4 (ξ creciente), que U5a excluye.

En cualquier caso, la ventana (Δℓ ≈ 1.3) no permite extrapolar a ℓ ≈ 70.
