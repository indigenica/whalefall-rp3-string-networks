# U24 · ¿Cambia la emisión si se drena el baño de otra manera? · Preregistro

*18 de octubre de 2026, 14:40 (Madrid). Escrito antes de producir ninguna red de U24 y antes de escribir su código.*

## Motivo

El artículo de la red (`paper/prd_network/`) llama «bath-free» a R_Γ = 1.56 tras drenar el baño masivo inicial (U10a, U11).
- **El drenaje es fijo:** la velocidad masiva se frena a ritmo físico Γ_d = m_r hasta ℓ_d = 3.8, solo donde el campo está cerca del vacío y cada 4 pasos (`kern.c`, `damp_massive`).
- **El único control de que no altera a las cuerdas es indirecto:** R_ξ cambia menos de un 1 %. Una densidad igual no garantiza una emisión igual.
- **En RP³, además, el intercambio X no está validado cuantitativamente** (V5 de U10a), así que no hay un segundo método que confirme lo que el drenaje quita.

La revisión de Codex del 18 de octubre (`paper/prd_network/PENDIENTE_ENVIO_PRD.md` §3) pide medir cuánto depende R_Γ del momento en que se apaga el drenaje y de su intensidad.

**Pregunta:** con todo lo demás igual, ¿cuánto cambian R_Γ, ε_RP³ y ε_U(1) al variar ℓ_d y Γ_d alrededor del protocolo de producción?

## Producción (`network/u24_produce.py`)

**Protocolo base, el de U11 a menor tamaño:**
- N = 384, τ_fin = N/2 = 192 y sfac = 1 (m_r·a·dx = 1 al final), con lo que ℓ_fin = ln 192 = 5.257;
- ξ_i = 0.5, ℓ_i = 1, dτ = 0.25;
- **45 instantáneas** (Δℓ ≈ 0.095) y espectros cada 3.
- Las mismas semillas dan las mismas condiciones iniciales en todas las variantes, así que la comparación se hace **semilla a semilla**.

**Variantes del drenaje** (RP³ y U(1), semillas 61, 62 y 63):

| Etiqueta | ℓ_d | Γ_d | Qué prueba |
|---|---|---|---|
| `D0` (referencia) | 3.8 | m_r | el protocolo de U11 |
| `Dearly` | 3.4 | m_r | apagar antes: queda más baño |
| `Dlate` | 4.1 | m_r | apagar después, más cerca de la ventana |
| `Dweak` | 3.8 | m_r/2 | drenar la mitad de fuerte |
| `Dstrong` | 3.8 | 2 m_r | drenar el doble de fuerte |

**Cambio de código.** `lattice.run` admite un argumento nuevo, `damp_rate` (por defecto `None`, que significa m_r como hasta ahora); se pasa a `set_damp(tau_d, damp_rate)`. Con `damp_rate=None` la evolución tiene que ser **idéntica bit a bit** a la actual. Antes de producir se comprueba con una red de N = 48 drenada: `None` frente a `m_r` explícito, todos los campos iguales. La prueba se añade a `network/test_network.py`.

**Orden y coste:**
- Se produce primero la semilla 61 de las cinco variantes en los dos modelos, luego la 62 y después la 63. El script salta las redes que ya existen, de modo que se puede parar y seguir.
- Medido en la nube con 4 hilos (U16, r = 0.5, τ_fin = 192): unos 32 min por red RP³ de 384³ y unos 8 min por red U(1). Son 15 + 15 redes, unas **10 h**.
- Salida: `network/results/u24/runs/<MODELO>_N384_s<semilla>_<variante>.json` y `network/results/u24/produce.log`.

## Análisis (`network/u24_analyze.py`, fijado ahora)

**Ventana fija W = [4.36, 5.21]:**
- empieza donde la de producción (576³: [4.36, 5.61]), de modo que la referencia `D0` deja el mismo hueco de 0.56 en ℓ entre el fin del drenaje y la ventana;
- 10 puntos equiespaciados para `analyze.window_ratio`, máscara rm2.5 y los mismos estimadores que U11 (local, ε_lat).

**Cantidades, por variante v:**
- R_Γ^v y R_ξ^v (cocientes de medias sobre semillas), y las medias ε_RP³^v, ε_U(1)^v, ξ_RP³^v, ξ_U(1)^v;
- **cambio relativo a la referencia:** δ_v = R_Γ^v / R_Γ^D0 − 1, con su error por **bootstrap emparejado**: se remuestrean semillas y la misma semilla entra en todas las variantes. Igual para R_ξ, ε_RP³ y ε_U(1).

**Veredicto (sobre R_Γ):**
- **DRENAJE_ROBUSTO:** para las cuatro variantes, |δ_v| ≤ 0.05 y σ(δ_v) ≤ 0.03. La expresión «bath-free» se mantiene y la dependencia del drenaje queda acotada en el 5 %.
- **DRENAJE_SENSIBLE:** en alguna variante, |δ_v| > 0.05 y |δ_v| > 2σ(δ_v). La mayor |δ_v| entra como sistemática del drenaje en el error de R_Γ, y el texto deja de decir «bath-free» sin matiz: diría «with the bath drained by our protocol».
- **INDETERMINADO:** cualquier otro caso, sobre todo si σ(δ_v) > 0.03. No falsa nada; se informa la potencia que faltó.

**Lecturas informativas, que no cambian el veredicto:**
- **Si el drenaje se lleva emisión de las cuerdas:** con `Dstrong` frente a `D0` bajarían ε_RP³ y ε_U(1). Si en cambio solo quitaba baño, deberían cambiar poco una vez el baño está agotado.
- **Si en la referencia queda baño:** con `Dlate` frente a `D0` bajaría ε. Así se lee cuánto baño queda en la ventana.
- **R_ξ^v:** se espera |cambio| < 1 %, como en U11.
- **D0 a 384³ frente a U11 a 576³:** distintos N y ℓ_fin; solo como comparación.

## Lo que U24 no resuelve

- **Es 384³ y no 576³.** La ventana dura 0.85 en ℓ en lugar de 1.25, con 3 semillas en lugar de 6/4. Mide la *sensibilidad* al drenaje en un régimen cercano al de producción; no repite U11.
- **No prueba drenajes de otra forma**, como otro umbral de «cerca del vacío» o drenar también cerca de los núcleos.
- **El muestreo de instantáneas (U18) es el mismo en todas las variantes**, así que se cancela en δ_v, pero no en los valores absolutos.

## Enmienda 1 · 18 de octubre de 2026, 20:53 (Madrid): sin veredicto con una semilla, errores informativos y notas de registro. Escrita tras ver la semilla 61 y antes de ver las semillas 62 y 63

Una auditoría independiente del código y de la semilla 61 (que reproduce los números de U24_RESULT.json hasta la última cifra) señala varios puntos. **Los criterios del veredicto no cambian.**

1. **Corrección de código.** Con una sola semilla, todas las réplicas del bootstrap son la misma semilla y σ sale 0. «σ ≤ 0.03» se cumplía entonces por defecto, y el análisis provisional decía «DRENAJE_ROBUSTO (provisional)».
   - Con menos de dos semillas, ahora no hay σ ni veredicto: **SIN_VEREDICTO (n = 1)**. Vale también para U25.
2. **Error informativo.** Con tres semillas, el bootstrap de semillas da un σ más pequeño que el error típico habitual: un factor ≈ √(2/3), con solo 10 remuestreos distintos.
   - Se mantiene como estaba preregistrado.
   - Se informa además, sin que cuente para el veredicto, el error típico entre semillas: la desviación de los δ por semilla, con ddof = 1, entre √n (`se_seeds`).
3. **La ventana de ε tiene 9 puntos, no 10.** El último punto (ℓ = 5.21) cae después de la última instantánea con derivada centrada (ℓ = 5.164).
   - R_Γ y las medias de ε promedian sobre [4.36, 5.12], y R_ξ sobre los 10 puntos.
   - Es lo mismo que en U11 y ya lo decía el código. Aquí se declara.
4. **Referencia a 576³, informativa.** Las redes de U11, evaluadas en esta misma ventana, dan R_Γ = 1.562. Con las tres semillas, R_Γ(D0) a 384³ se comparará con ese valor. Si difieren en más de 2σ combinados, se anota como **efecto de tamaño finito** de 384³ en el informe. No cambia el veredicto de U24, que mide diferencias entre variantes con la misma semilla.
   - Con la semilla 61, R_Γ(D0) = 1.20. Una sola semilla de 576³ dispersa ≈ 18 % en el mismo τ/N, así que ese valor está a ≈ −1.3σ.
5. **Lo que el drenaje toca, ya visible en la semilla 61 (informativo).** La expectativa «R_ξ cambia < 1 %» no se cumple:
   - con `Dstrong`, ξ_U(1) sube un 13 % y R_ξ cambia un −9.7 %;
   - con `Dweak`, R_ξ cambia un +4.1 %.
   - Las emisiones siguen a ξ, y R_Γ se mueve menos porque los dos modelos se compensan en parte. El umbral de «cerca del vacío» del drenaje incluye probablemente la parte externa de los núcleos de las cuerdas en movimiento.
   - El informe lo dirá, sea cual sea el veredicto, y el artículo no usará «R_ξ no cambia» como prueba de que el drenaje no toca a las cuerdas.
6. **Notas de registro.**
   - Este ESPEC se subió en el commit `55070bd`, a las 14:20 (Madrid). La cabecera dice 14:40 por error.
   - El analizador se fijó en `109ce51`, a las 14:30 (Madrid): dos minutos después de arrancar la producción y una hora antes de la primera red terminada.
   - Las horas de `produce.log` están en UTC: la producción empezó a las 14:28 de Madrid.
   - El coste real es ≈ 55–66 min por red RP³ y ≈ 11 min por red U(1), frente a 32 y 8. U24 tarda ≈ 17 h y U25 ≈ 7 h, no 10 y 4.
