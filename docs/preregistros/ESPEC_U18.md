# U18 · ¿Por qué difieren los dos estimadores de R_Γ al cambiar la resolución? · Preregistro

*14 de octubre de 2026, 00:45 (Madrid). Escrito antes de producir ninguna red de U18 y sin haber analizado de nuevo las de U16.*

## Lo que se sabe antes (U16, `docs/informes/INFORME_U16.md` §1 y §5)

Al resolver mejor el núcleo (m_r·a·dx de 1 a 0.5), con condiciones iniciales emparejadas:
- **el estimador local de R_Γ** (derivada centrada de a⁴ρ_a sobre tres instantáneas, promediada en la ventana) sube un +4.0 % ± 0.25 %, con un zigzag de ±2 % (bajada en r = 0.6) que se repite en las tres semillas;
- **el integral** (Δ(a⁴ρ_a) entre los bordes de la ventana) no cambia (−0.5 %).

Las instantáneas son 30 entre ℓ_i = 1 y ℓ_fin = 4.56: una cada Δℓ ≈ 0.12.

**Hipótesis, fijada ahora.** a⁴ρ_a lleva una oscilación coherente: la condición inicial con π = 0 deja un baño masivo en fase en toda la caja (auditoría U2, hallazgo 1), y el estimador cinético la ve. Su periodo en ℓ es Δℓ ≈ π e^(−ℓ)·m_r/ω, del orden de 0.04–0.07 en la ventana, más corto que el paso entre instantáneas. La derivada local la muestrea con aliasing. La dispersión del leapfrog, que depende de m_r·a·dτ y por tanto de r, desplaza la fase de esa oscilación, y con ella el aliasing. Si es eso:
- con instantáneas mucho más densas, el estimador local debe dejar de depender de r y acercarse al integral;
- a⁴ρ_a debe mostrar una oscilación resoluble.

## Producción (`network/u18_produce.py`)

- **Semilla 61**, RP³ y U(1), con las mismas condiciones emparejadas de U16 en dos resoluciones: (N, τ_fin, r) = (192, 96, 1) y (320, 160, 0.6).
- **Instantáneas:** las 30 de U16 (los mismos τ) más una rejilla densa de 240 τ geométricos en el mismo intervalo. Con la unión, las 30 originales siguen ahí.
- **Espectros** solo en las 30 originales. La evolución no depende de cómo se trocee ni de las medidas, que solo leen los campos.
- **Validación V1:** con solo las 30 originales, R_ξ y R_Γ (local e integral) de cada resolución deben reproducir los de la semilla 61 de U16 a ≤ 10⁻⁶ relativo. Si no, U18 se informa como **INVALIDO** y no se lee nada más.
- **Coste:** unas 2 h con 2 hilos.

## Análisis (`network/u18_analyze.py`)

Para cada resolución:
- R_Γ local con las 30 instantáneas (= U16);
- R_Γ local con todas;
- R_Γ local con todas, pero con la derivada suavizada sobre una ventana de ±0.15 en ℓ (regresión lineal local);
- R_Γ integral.

Además, la amplitud y el periodo (en ℓ) del residuo de a⁴ρ_a respecto a un ajuste suave (polinomio de grado 3 en ℓ en la ventana), en RP³ y en U(1).

**Cambios con la resolución:** Δ_X = X(0.6) − X(1), de la semilla 61. En U16 eran +0.023 (local, 30 instantáneas) y −0.012 (integral).

**Clasificación:**
- **MUESTREO:** con las instantáneas densas, |Δ_local,denso| ≤ 0.010 y |R_Γ local,denso − R_Γ integral| a r = 1 y a r = 0.6 baja a menos de la mitad que con 30 instantáneas. El efecto de U16 era aliasing y el estimador integral es el bueno.
- **FISICO:** con las instantáneas densas, Δ_local,denso ≥ +0.015 (al menos dos tercios de lo de U16). El efecto es real en la tasa local; el integral lo diluye. Se dirá si viene de alguna parte concreta de la ventana.
- **INDETERMINADO:** cualquier otro caso.

**Lo que se informa en cada caso:** si hay oscilación resoluble (amplitud del residuo > 3 veces su ruido entre instantáneas vecinas), su periodo en ℓ y si se desfasa entre r = 1 y r = 0.6.

## Lo que U18 no resuelve

- **Una sola semilla.** En U16 las tres semillas dieron el mismo zigzag, así que una basta para distinguir aliasing de efecto real. No sirve para reducir el error de R_Γ.
- **No se prueba r = 0.5** (384³), por coste: 240 medidas a 384³ son unas 4 h.


## Corrección de horas · 14 de octubre de 2026, 00:44 (Madrid)

La hora de la cabecera (00:45) está mal: el preregistro se fijó en el commit `c7ff41b` a las 00:32:50. La producción de U18 se lanzó después, hacia las 00:35. El orden no cambia.
