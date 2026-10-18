# U2 · ¿Cuánto emite de verdad la red RP³? Medida de ε en la red y cierre físico · Preregistro

*6 de octubre de 2026, 12:30 (Madrid) · claude_dual/unification/u2_red. Escrito antes del código. Las referencias se consultaron antes de escribir.*

## 0. Pregunta

U1 dejó el cierre cosmológico pendiente de un solo número: la eficiencia ε con la que la red de cuerdas ℤ₂ de RP³ emite Goldstones libres. Se mide respecto a la normalización GHN que usa el libro oscuro:

Γ_GHN = 8π H³ f² ℓ · ξ_ref(ℓ), con ξ_ref = 0.24 ℓ y ℓ = ln(m/H).

U1 exige ε ≲ 10⁻³–6×10⁻². U2 mide ε en una simulación 3D de la teoría de campos y vuelve a ejecutar el cierre con la ε medida.

## 1. Modelos (los mismos en el mismo código)

**RP³ (el de la acción común, `cierre_accion/MODELO.md`).** Q es una matriz real simétrica 4×4, con 10 componentes en una base de Frobenius ortonormal. El lagrangiano, dividido por Z_Q:

ℒ = ½Tr(Q̇²) − ½Tr(∇Q·∇Q) − U(Q)/Z_Q,  U = λ(TrQ − 1)²/4 + λ_P Tr[(Q² − Q)²]/4,

con λ = 320, λ_P = 192, Z_Q = 400.
- El vacío es Q = nnᵀ con n ∈ S³, es decir RP³.
- Espectro: 3 Goldstone; un quintete y un singlete con m² = 0.24; un modo pesado con m² = 1.84. Los comprueba V2.
- Constante de Goldstone efectiva: f² = ½ en unidades de código, igual a (F/2)², que es la f_a = F/2 del libro.

**Control U(1).** φ ∈ ℝ², V = (λ_u/4)(φ² − v²)², con v² = ½ (misma f) y λ_u = 0.24 (m_r² = 0.24, igual que el quintete). La masa de núcleo y la normalización de la tensión, μ = πf²ℓ, son idénticas.

**Dinámica.** Era de radiación, a = sτ en tiempo conforme, cuerdas físicas (sin el truco PRS). Leapfrog simplectico con π = a²φ'. Laplaciano de 7 puntos. dx = 1, dτ = 0.25. Se elige s = 1/(m τ_fin), de modo que m·a·dx ≤ 1 y ℓ = ln(τ²/τ_fin).

**Condición inicial.**
- En τ_i = 2: campo gaussiano suavizado (P(k) ∝ e^{−k²l_c²/2}, l_c = 2 celdas), proyectado al vacío.
- **RP³:** n es el autovector principal de una matriz simétrica gaussiana, con signo arbitrario, así que las cuerdas ℤ₂ aparecen en las líneas de degeneración. Q = nnᵀ, Q̇ = 0.
- **U(1):** φ = v·ĝ, φ̇ = 0.

## 2. Validaciones previas (V1 a V4; sin ellas no hay producción)

- **V1, energía en Minkowski** (a = 1, sin fricción) con la condición inicial no lineal, 400 pasos: deriva relativa de la energía total < 10⁻³ en los dos modelos.
- **V2, espectro lineal** alrededor del vacío: frecuencias de un modo k de prueba en cada sector (Goldstone, quintete, singlete, pesado) frente a ω² = k_lat² + m², con error < 10⁻³ (tras corregir la dispersión de la red y del leapfrog).
- **V3, detector de cuerdas:**
  - Una cuerda recta ℤ₂ (n gira π alrededor de z) da exactamente una plaqueta perforada por rebanada z.
  - El vórtice U(1) igual.
  - Un estado uniforme da 0.
- **V4, estimador de radiación sin cuerdas.** Un paquete de Goldstone libre en era de radiación debe cumplir:
  - a⁴ρ_a constante dentro de un 3 % una vez dentro del horizonte (es decir, Γ_falsa ≈ 0);
  - no se cuenta ninguna cuerda.

## 3. Medidas (en ~24 instantes equiespaciados en ln τ)

- **Densidad de cuerdas.** ξ = L_fís t²/V_fís = L_c τ²/(4V_c), con L_c = (2/3)·N_plaq·dx. Una plaqueta está perforada cuando la holonomía ℤ₂ Π sign(n_i·n_j) = −1 (RP³) o el devanado ≠ 0 (U(1)).
- **Radiación de Goldstones libres.** ρ_a = ⟨2K_G⟩ en las celdas a una distancia física mayor que r_m = 2.5/m de cualquier plaqueta perforada.
  - K_G = |Q'n|² − (n·Q'n)², entre a⁴ (parte tangente al vacío), para RP³.
  - K_G = (φ₁φ₂' − φ₂φ₁')²/(2|φ|²a⁴) para U(1).
  - Es el estimador de GHN: radiación libre = 2 × energía cinética fuera de los núcleos.
  - Validez: fracción no enmascarada ≥ 0.5.
- **Emisión.** Γ = H·a⁻⁴ d(a⁴ρ_a)/d ln τ (en RD, d ln τ/dt = H), por diferencias centradas suavizadas sobre tres instantes.
  - ε_lat(ℓ) = Γ/(8πH³f²ℓ·0.24ℓ).
  - η(ℓ) = Γ/(8πH³f²ℓ·ξ_med), la eficiencia por unidad de energía de cuerda. Debe ser O(1) por conservación de energía.
- **Radiación masiva.** ρ_M = ⟨2K_M⟩ fuera de los núcleos, como diagnóstico de la fracción que no va a Goldstones.
- **Espectro (secundario).** Pendiente q del espectro instantáneo F(k) ∝ k^{−q} en 10H < k < m/4 (física), a partir de la variación entre instantes del espectro de la velocidad Goldstone enmascarada.

## 4. Producción

| Conjunto | Modelo | N | τ_fin | ℓ_fin | Semillas | Uso |
|---|---|---|---|---|---|---|
| A | RP³ y U(1) | 256 | 128 | 4.85 | 4 + 4 | principal |
| B | RP³ y U(1) | 192 | 96 | 4.56 | 2 + 2 | volumen (caja/horizonte igual, más ruido) |
| C | RP³ y U(1) | 256 | 128, s × 0.6 | 4.34 | 1 + 1 | resolución del núcleo (m·a·dx ≤ 0.6) |

Ventana de análisis: ℓ ∈ [3.0, ℓ_fin], con la fracción no enmascarada ≥ 0.5. **Producto principal:** los cocientes R_ξ = ξ_RP³/ξ_U1 y R_Γ = Γ_RP³/Γ_U1 al mismo ℓ, que se calculan sobre las medias de cada conjunto. La incertidumbre combina:
- la dispersión entre semillas;
- la diferencia entre A y B;
- la diferencia entre A y C;
- la sensibilidad a r_m ∈ {1.5, 4}/m.

Se suman en cuadratura.

## 5. Del ℓ de la red al ℓ físico (~50–87): dos ramas de la literatura

Se evalúan las dos; ninguna se elige a posteriori.
- **B1, GHN** (Gorghetto–Hardy–Villadoro, arXiv:2007.04990): ξ_U1 ≈ 0.24ℓ, es decir ε_U1 = 1. Entonces **ε_RP³ = R_Γ**, constante.
- **B2, escalado estándar** (Hindmarsh et al., arXiv:2102.07723: ζ\* = 1.20 ± 0.09): ξ_U1 → 1.20. Entonces **ε_RP³(ℓ) = R_Γ · 1.20/(0.24ℓ)**, que depende de la época.
- Si R_Γ muestra tendencia significativa con ℓ dentro de la ventana (> 2σ), se añade una tercera rama: extrapolación lineal en 1/ℓ hasta ℓ_fís.

**Espectro.** Se evalúan los dos brazos de U1 (ir_delta y energy_logflat, x0 = 10 y 30). El q medido indica cuál se parece más, pero en la red (q ≲ 1) no se puede extrapolar con fiabilidad. Se informa.

**Cierre.** El libro de U1 se vuelve a ejecutar con la fuente multiplicada por ε(ℓ) de cada rama, en los extremos ±1σ de R_Γ, con K_sw ∈ {0, 1} y las mismas puertas P1–P5 (P6 informativa).

## 6. Veredictos preregistrados

- **CIERRE_FÍSICO:** las dos ramas, con R_Γ en su extremo +1σ, quedan permitidas en al menos un brazo espectral con K_sw = 0.
- **CIERRE_SEGÚN_RAMA:** solo una rama queda permitida (con R_Γ al +1σ). La ballena se cierra si la literatura resuelve en ese sentido.
- **FALSADA_EN_FORMA_ACTUAL:** ninguna rama queda permitida ni siquiera con R_Γ al −1σ. El sector oscuro de la ballena 5 exige un cambio: otro mecanismo de masa, otra f, o una red que no escale.
- **INCONCLUSO:** las bandas ±1σ cruzan el borde.

**Puerta de validez del control:** el U(1) tiene que dar η ∈ [0.3, 3] y ξ_U1(ℓ = 4.5) ∈ [0.5, 2.5] (rango de la literatura en ℓ ≈ 4–5). Si no, el estimador no es fiable y el veredicto pasa a INCONCLUSO.

## 7. Qué no hace (declarado)

- No resuelve el ℓ físico. La extrapolación es de la literatura y está en disputa, por eso hay dos ramas.
- No usa AMR ni la radiación covariante del protocolo de `cierre_cosmologia/SIGUIENTE_PASO.md`. Es el método estándar de la literatura (GHN, Hindmarsh et al.) con un control U(1) en el mismo código: la parte robusta es el cociente.
- Mallas pequeñas (256³, ℓ ≤ 4.9, una caja de unos 2 horizontes al final). Si el resultado está en el borde, se prepara un job más grande para WSL o GPU.
- No hay gravedad local, NPQT, χ ni acoplamiento ξ_R R S. Q es espectador en radiación pura, lo cual es válido mucho después del recalentamiento.

## Enmiendas antes de la producción (6 oct, 13:00; escritas tras V1–V4 y antes de cualquier red de producción)

1. **V1 falla tal como estaba preregistrado.** Deriva máxima: 2.7 % (RP³) y 2.5 % (U(1)).
   - Es un desplazamiento constante que se fija en el primer paso y escala exactamente como dτ²: 2.7e-2, 6.7e-3 y 1.7e-3 con dτ = 0.25, 0.125 y 0.0625.
   - La variación secular posterior es < 1e-3: 6.5e-4 (RP³) y 7.6e-4 (U(1)).
   - Es el desfase de la energía sombra del leapfrog con una condición inicial rugosa a escala de red; el integrador no deriva. Se mantiene dτ = 0.25 y se añade el **conjunto D**: N = 128, τ_fin = 64, una semilla por modelo, con dτ = 0.25 y con 0.125.
   - La diferencia en R_Γ entre los dos dτ entra en la incertidumbre de la §4.
2. **V3, prueba de devanado completo.** El devanado completo de n (no topológico en RP³) da holonomía ambigua solo en las plaquetas que contienen el punto singular, porque ahí los pasos son de 90°. Es inevitable a escala de red. La prueba pasa a ser «ningún flujo ℤ₂ a más de 2 celdas de los núcleos»: pasa.
   - La prueba de cuerda recta pasa a ser un par cuerda-anticuerda, porque una cuerda sola no es periódica: da exactamente 2N.
3. **Radios de máscara.** Se limitan a N/4 celdas. Los instantes con r_m > N/4 se marcan como no válidos, lo que solo afecta a ℓ < 1.
4. **Condición inicial (cambia la §1; se decide tras un piloto con N = 128, archivado en `_superseded/pilot_ic_v0/`).**
   - **Qué falló.** Con τ_i = 2 y l_c = 2, Q decae casi a cero mientras H ≫ m y luego se ordena de golpe en ℓ ≈ 2.5. La red nace con ξ ≈ 20 y colapsa hasta ξ ≈ 0.5 justo dentro de la ventana de análisis. Toda la ventana sería un transitorio, lejos del escalado.
   - **Nueva condición.**
     - Se empieza en ℓ_i = 1, es decir τ_i = √(e·τ_fin/s_fac), con el campo ya en el vacío.
     - La longitud de correlación se elige para que el ξ inicial sea ≈ 0.5 en los dos modelos. La calibración (5 valores de l_c, N = 128) da ξ·l_c² = 50 (RP³) y 27.2 (U(1)) en τ_i = 18.65.
   - **Conjunto E nuevo** (atractor): ξ_i ∈ {0.25, 1.0}, una semilla por modelo, N = 256. Comprueba si el cociente al mismo ℓ depende de la condición inicial; la diferencia con A entra en la incertidumbre.
   - Las semillas de A pasan a 4 (U(1)) y 2 + 2 (RP³), por orden de ejecución.
5. **Estimador integral complementario** (declarado antes de ver el conjunto A, después del piloto D de N = 128, cuyas derivadas locales fluctúan entre instantes).
   - Definición: ε_int = Δ(a⁴ρ_a) / ∫ 8πH³f²ℓ·0.24ℓ·a⁴ dt sobre la ventana ℓ ∈ [3.0, ℓ_fin − 0.05]. Es la integral temporal de Γ.
   - Para el mismo N, τ_fin y s_fac, la normalización es idéntica en los dos modelos. Así, R_Γ,int = Δ(a⁴ρ_a)_RP³ / Δ(a⁴ρ_a)_U1.
   - El estimador preregistrado (§3) sigue siendo el principal. Si los dos difieren en más de la σ total, la diferencia se suma como sistemático.
