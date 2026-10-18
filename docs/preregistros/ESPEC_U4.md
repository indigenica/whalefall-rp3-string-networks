# U4 · Salir de dudas: CLASS para Lyman-α y redes grandes para el atractor · Preregistro

*6 de octubre de 2026, ≈21:40 (Madrid). Escrito antes de ejecutar nada. Las dos pruebas corren en la máquina WSL del equipo (`env/INSTALAR_WSL.md`), no en la nube.*

U3 dejó el sector oscuro con dos dudas:
- (i) la prueba de Lyα se hizo con un resolvedor propio, no validado según su preregistro;
- (ii) las redes de 256³ no llegan al atractor de escalado.

U4a resuelve (i) con CLASS. U4b ataca (ii) con redes más grandes, hasta donde dé la memoria.

---

## U4a · Lyman-α con CLASS (`lya/u4a_class.py`)

**Modelo.**
- La cosmología es la misma de U3a: h = 0.674, ω_b = 0.0224, ω_dm = 0.120, n_s = 0.9649, A_s = 2.1×10⁻⁹, N_ur = 3.044.
- La parte fría (desalineamiento) es CDM. La tibia (emisión) es una especie ncdm con su distribución de momentos en fichero (`use_ncdm_psd_files`), construida como en U3a: KDE en ln q de los paquetes, con la cola de número por encima del 99.9 % descartada.
- Por invariancia de escala se usa una masa sustituta de 1 keV con la misma distribución de q/m. CLASS renormaliza la amplitud a la Ω pedida, porque se le pasan m_ncdm y Ω_ncdm a la vez.
- P(k) en z = 3, con k_max = 1000 h/Mpc.

**Validaciones (sin ellas no hay veredicto):**
- **C0, lectura del fichero:** un fichero de distribución Fermi–Dirac tiene que dar la misma P(k) que la ncdm térmica analítica de CLASS (< 0.5 % en 0.5–50 h/Mpc).
- **C1, precisión:** con ℓ_max_ncdm = 50, el doble de puntos de momento y tolerancias 10 veces menores, δA de dos filas decisivas cambia < 2 %.
- **C2 (informativa):** CLASS frente a `lya/boltz.py` en las mismas filas, en T² y en δA. Es la validación que le faltó a U3a.

**Puertas por fila** (filas de S6 con K_sw ∈ {0, 1} que pasan P1–P4):
- **LYA-δA:** δA (Murgia et al. 2017) ≤ δA de la referencia, con umbral estricto (mín.) y amplio (máx.) como en U3a, pero calculados con CLASS.
- **LYA-forma (puerta principal, nueva).** Se ajusta T²(k) en 0.5–20 h/Mpc con la familia de mezclas térmicas (fracción f_w, masa m_w) calculada con CLASS, sobre una malla de 9 masas (0.5–10 keV) por 9 fracciones (0.05–1).
  - **Permitida (estricta)** si f_w ≤ 0.14 (m_w/keV)^1.1 (ajuste de Garcia-Gallego et al. 2025, 2σ) para m_w ≤ 4 keV, o si m_w ≥ 5.3 keV para cualquier f_w (Iršič et al. 2017).
  - **Permitida (amplia)** con la tabla de ese artículo: f_w ≤ 0.16, 0.35, 0.50 y 0.67 para 1, 2, 3 y 4 keV, interpolada en log m.
  - Si el mejor ajuste tiene un residuo rms en ln T² > 0.1, la fila se marca FORMA_NO_REPRESENTABLE y decide δA.
- **Veredicto U4a:** la lógica de ESPEC_U2 §6 con P5 sustituida por LYA-forma estricta. Se informa también con δA y con los umbrales amplios.

---

## U4b · Redes grandes (`network/u4b_*.py`)

**Hecho conocido (U2, U3b).** ξ(ℓ) a un ℓ dado depende de la historia en ℓ, no del tamaño de la caja. Por eso una caja mayor no cambia nada en ℓ = 4.5. Lo que sí hace es **alargar la ventana**: ℓ_fin = ln(N/2), que da 5.55 para N = 512, 5.95 para 768 y 6.24 para 1024.

**Plan, adaptado a la máquina** (`network/u4b_plan.py` elige el N mayor que cabe en el 75 % de la RAM, con el consumo medido en la nube):
- **Fase 1 (solo U(1), 5 veces más barata que RP³):** al N máximo, ξ_i ∈ {0.25, 2.0} (semilla 61) y ξ_i = 0.5 (semillas 61 y 62).
- **Fase 2 (solo si la fase 1 pasa sus puertas):** RP³ al N mayor que quepa, con ξ_i = 0.5 y semillas 61 y 62, y U(1) al mismo N si es distinto del de la fase 1.

**Puertas (al final de la ventana, ℓ_ref = ℓ_fin − 0.3):**
- **G1′:** ξ_U1(ℓ_ref) ∈ [0.5, 2.5].
- **G2′:** memoria de la condición inicial |ξ(2.0) − ξ(0.25)| / media < 20 % en ℓ_ref.
- **G3′ (informativa):** pendiente dξ/dℓ del U(1) en [ℓ_ref − 0.7, ℓ_ref], frente a 0.24–0.25 (GHV, Buschmann et al.) o 0 (Hindmarsh et al.). Es el primer dato propio sobre la disputa.

**Si G1′ y G2′ pasan y hay redes RP³:**
- R_Γ en la ventana [ℓ_ref − 1, ℓ_fin − 0.05], con el bootstrap de semillas y los sistemáticos de U2 (máscara, dτ).
- Pendiente espectral q de RP³ a ℓ_fin.
- S6 y U4a se recalculan con ese R_Γ. Si |ΔR_Γ| > 0.3, se repite U4a.

**Si no pasan:** INCONCLUSO, y se informa del ℓ al que llegó y de la tendencia de la memoria con N (256 frente a N_max). Si la memoria baja con N, se extrapola el N necesario, como dato y no como veredicto.

---

## Qué decide U4, sobre el sector oscuro de la ballena 5

| U4a (CLASS, forma) | U4b (red) | Lectura |
|---|---|---|
| B2/B3 permitidas | gates pasan y R_Γ ≈ 1.7 se mantiene | **CIERRE_SEGÚN_RAMA, firme.** Queda pendiente solo la disputa de la literatura (B1 frente a B2), y G3′ aporta un dato propio |
| B2/B3 excluidas | — | **FALSADA_EN_FORMA_ACTUAL** (el sector oscuro necesita otro mecanismo) |
| B2/B3 permitidas | gates fallan | sigue INCONCLUSO por la red; Lyα deja de ser una duda |

---

## Enmienda 1 · 7 de octubre de 2026, ≈15:00 (Madrid) · cuadratura en momentos de las filas de U4a

*Escrita antes de tener el δA de ninguna fila con la configuración definitiva, salvo la fila de prueba que se indica abajo.*

**Qué pasó.** La primera ejecución en WSL (`make u4a-quick`, 7 oct por la mañana) superó C0 (0.37 %) y calculó las referencias térmicas. Después, ninguna de las 10 filas terminó en más de dos horas. Se reprodujo en la nube con CLASS v3.4.1 (el mismo commit, 814ce59) y con los mismos ficheros de distribución, idénticos byte a byte:
1. **Muestreo.** El muestreo adaptativo de CLASS (estrategia 0) pone 230 nodos de momento en estas distribuciones; una térmica usa 8. La causa es el borde abrupto a q bajo (dlnN/dlnq ≈ 20) y la cola dN/dln q ∝ 1/q, que ocupa tres décadas. Con la jerarquía completa hasta k = 1000 h/Mpc, el coste crece como N^≈2.7, así que cada fila tardaría días.
2. **Aproximación de fluido.** La de CLASS (`ncdm_fluid_approximation = 2`) no sirve aquí. Toma la velocidad del sonido de ⟨v²⟩ de toda la especie, que con esta cola crece con el q más alto muestreado. En la fila de prueba, δA = 0.221, 0.287 y 0.359 con 20, 40 y 80 nodos, sin converger. Se mantiene la jerarquía completa (`= 3`, que en CLASS es `ncdmfa_none`), como ya estaba en el código.
3. **Reglas fijas.** Las reglas fijas de CLASS (Laguerre, trapecio en q, trapecio en t = 1/(1+q)) no pueden repartir los nodos de forma uniforme en ln q, y convergen de forma errática: ±40 % en δA con 20–40 nodos, y < 1 % solo a partir de unos 120 nodos (comprobado con `lya/boltz.py`).

**Cambio.**
- Las filas usan una regla del trapecio **uniforme en ln q sobre el rango del fichero**, con **45 nodos**: `ncdm_quadrature_strategy = 4`.
- Esa regla se añade a CLASS con `env/class_qm_trapz_log.patch`, 18 líneas que solo tocan `tools/quadrature.c` y el `enum` de `include/quadrature.h`. La física de CLASS no cambia.
- **Nueva validación C0b:** el fichero Fermi–Dirac de C0, con estrategia 4 y 20 nodos, tiene que reproducir la ejecución analítica en lo que usa U4a: δA con un error relativo < 0.5 % y T² a < 0.5 % en 0.5–20 h/Mpc, el rango de la puerta de forma. Un CLASS sin parche falla aquí (y en `env/check_env.py`).
  - Prueba en la nube: δA = 0.4773 frente a 0.4773; T² a ≤ 0.05 % hasta k = 20 y 1.3 % en k = 50, donde T² ≈ 0.14.
  - Con 30 nodos sobre [0.01, 30], 0.8 % en k = 50. Por eso C0b no exige el 0.5 % hasta 50 h/Mpc de C0.
- **C1 no cambia:** el doble de nodos (90), ℓ_max_ncdm = 50 y tolerancias diez veces menores.
- En `--quick` no se ejecutan C1 ni C2.

**Pruebas que informaron el cambio.** Son todas sobre una sola fila, B2 central, espectro log-plano, x₀ = 10, K_sw = 0. Quedan registradas y no cuentan como resultado.

| Configuración | δA (k_max = 1000) | T² en k = 5 / 10 / 20 / 50 h/Mpc |
|---|---|---|
| CLASS, fluido (2), t-trapecio, 20 / 40 / 80 nodos | 0.221 / 0.287 / 0.359 | — |
| CLASS, jerarquía completa, t-trapecio, 20 / 40 nodos | 0.134 / 0.156 | — |
| CLASS, jerarquía completa, t-trapecio en unidades de 3× la mediana, 40 nodos | 0.185 | 0.968 / 0.947 / 0.905 / 0.796 |
| **CLASS, jerarquía completa, ln q, 45 nodos (parche)** | **0.179** | **0.980 / 0.957 / 0.912 / 0.798** |
| `boltz.py` con los mismos 45 nodos en ln q | — | 0.981 / 0.959 / 0.917 / 0.809 |

- **Convergencia en `boltz.py`** (δA con k_max = 120, solo para comparar entre sí): 30 / 45 / 60 / 90 / 150 / 600 nodos en ln q dan 0.1543 / 0.1574 / 0.1573 / 0.1573 / 0.1573 / 0.1573.
- **Con los mismos nodos y sin aproximación de fluido, CLASS y `boltz.py` coinciden** a ≤ 1.4 % en T² hasta k = 50 h/Mpc. También coinciden en una componente térmica caliente (0.5 keV, f = 0.02), a ≤ 0.3 %.
- **Lectura.** La discrepancia de U3a con el ajuste de Viel (W1) es de `boltz.py` con las distribuciones térmicas a k alto (T² en k = 50: 0.18 frente a 0.135 de CLASS para 3.5 keV). No aparece en estas distribuciones anchas.
- **Fila de prueba frente al umbral.** Su δA con CLASS (0.179) queda por debajo del umbral estricto calculado con CLASS (0.234, de la ejecución de la mañana). La configuración se eligió por convergencia, frente a la referencia de 600 nodos, no por su posición respecto al umbral. La misma fila daba 0.197 con `boltz.py` en U3a.
