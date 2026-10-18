# U3 · Las dos pruebas que deciden el sector oscuro · Preregistro

*6 de octubre de 2026, ≈18:17 (Madrid) · claude_dual/unification/u3. Escrito antes del código (el kernel de U3b es de las 18:18).*

*Nota: las horas de las enmiendas se corrigieron con las marcas de tiempo de los ficheros; la primera versión llevaba horas mal estimadas, de 18:40 a 20:05.*

Hay dos pruebas: U3a sustituye el indicador λ_fs de Lyα por una transferencia de Boltzmann lineal; U3b repite la medida de U2 con redes preparadas para llegar al atractor de escalado.

**Orden de ejecución.**
1. U3b (redes).
2. Si U3b pasa sus puertas, se recalcula S6 con su R_Γ.
3. Después U3a sobre las filas de S6 resultantes. Si U3b no pasa, U3a se aplica sobre las filas de S6 de U2.

---

## U3b · Redes en el atractor (preparación con núcleo de anchura comóvil fija)

**Problema (U2, auditoría).** Arrancando con cuerdas físicas en ℓ = 1, la red pasa un transitorio y la ventana de análisis [3, 4.85] empieza en el mínimo del colapso. El control U(1) queda en ξ = 0.45 < 0.5.

**Método (crecimiento de núcleo de tipo PRS, como en Hindmarsh et al. 2021).**
- Hasta τ_c, el término de potencial usa un factor de escala congelado, a_p = s·τ_c, en lugar de a. Así la anchura del núcleo es constante en coordenadas comóviles: w = 1/(m s τ_c) = τ_end/τ_c celdas. La fricción y el resto de la ecuación no cambian.
- A partir de τ_c, a_p = a: cuerdas físicas, exactamente como en U2. La anchura es continua en τ_c.
- Durante la fase gruesa, el logaritmo efectivo ℓ_eff = ln(m a_p τ) crece desde ℓ_eff(τ_i) = 1. Con eso la red tiene un rango en ln τ mucho mayor (τ_i = e·w) para relajarse antes de la ventana.
- **Parámetros:**
  - w = 3, N = 256, τ_end = 128;
  - τ_c = 42.67, con ℓ_c = ln(τ_c/w) = 2.65;
  - τ_i = 8.15, dτ = 0.25;
  - mismo kernel, mismas medidas, mismas máscaras.
- **Condición inicial:** la de U2 (campo en el vacío, l_c escalado con ξ·l_c² = K_ξ τ_i²/18.65²).

**Conjunto F** (5 redes por modelo como máximo, según el tiempo):
- F1: ξ_i = 0.5, semillas 41 y 42, RP³ y U(1).
- F2: ξ_i = 2.0, semilla 43, RP³ y U(1).

**Piloto** (N = 128, τ_end = 64, una semilla por modelo). Solo comprueba que el código funciona: la continuidad de la energía en τ_c y que las cuerdas se detectan en la fase gruesa. Si el piloto obliga a cambiar w o τ_i, se escribe como enmienda antes de lanzar F.

**Puertas de U3b:**
- **G1, control U(1):** ξ_U1(ℓ = 4.5) ∈ [0.5, 2.5] y η_U1 ∈ [0.3, 3]. Es la misma puerta de U2.
- **G2, atractor:** en ℓ = 4.5, |ξ(F2) − ξ(F1)| / media < 0.20 en los dos modelos. Es decir, la memoria de la condición inicial es menor del 20 %.

**Si pasan G1 y G2:**
- R_Γ se calcula sobre F1 (en la ventana [3.0, 4.85], con bootstrap de semillas).
- La diferencia F1–F2 se suma como sistemático, junto a los sistemáticos de U2 (caja, núcleo, dτ y máscara, con los valores medidos en U2).
- S6 se recalcula con ese R_Γ y se aplica la lógica de veredictos de ESPEC_U2 §6.

**Si no pasan:** el veredicto de la red sigue siendo INCONCLUSO y se informa de cuál falla.

---

## U3a · Lyman-α con una transferencia de Boltzmann lineal

**Problema.** La puerta P5 de U1 (λ_fs medio frente a un WDM térmico de 5.3 keV, con una franquicia del 10 %) es un indicador de orden de magnitud. La versión exploratoria con la cota de materia oscura mixta mostró que el veredicto depende de esa puerta.

**Herramienta.** CLASS no se puede instalar en esta sesión: PyPI y GitHub responden 403 a través del proxy, también desde el equipo. Se escribe un **resolvedor lineal de Einstein–Boltzmann** propio, en gauge síncrono (Ma & Bertschinger 1995):
- CDM, con los bariones tratados como fríos;
- radiación sin masa (fotones y neutrinos, N_eff = 3.044) con jerarquía ℓ_max = 8 y la aproximación de radiación libre (RSA, Blas et al. 2011) a partir de kτ = 30;
- especie no fría con su distribución de momentos completa f0(q): jerarquía Ψ_l(q), ℓ_max = 12 y 48 nodos log en q.

**Validaciones previas (sin ellas no se aplica la puerta):**
- **W1:** un WDM térmico de 2, 3.5 y 5.3 keV debe reproducir el ajuste de Viel et al. 2005 en k_{1/2} con un error < 10 %, y |T/T_Viel − 1| < 0.10 donde T_Viel > 0.3.
- **W2:** una especie no fría de masa enorme (fría) debe dar el cociente frente a CDM igual a 1 con error < 10⁻³.
- **W3:** duplicar ℓ_max y N_q debe cambiar δA en < 2 %.

**Distribución de las partículas emitidas.**
- La de los paquetes de `s4_dark.source` para el nodo, como densidad de núcleo gaussiano en ln q (anchura 0.15).
- La cola de q por encima del cuantil 99.9 % del número se descarta. Se informa de su fracción de masa.
- Por invariancia de escala, se usa una masa sustituta de 1 keV con la misma distribución de q/m y la misma densidad: es gravitatoriamente idéntica, y la presión cuántica es despreciable para m_G ≫ 10⁻²⁰ eV.
- Las perturbaciones iniciales son adiabáticas en a_i = 10⁻⁹. Los modos de interés están fuera del horizonte cuando termina la emisión, en T_stop ≈ MeV.

**Estimador de Lyα:** δA de Murgia et al. 2017.
- P_1D(k) = (1/2π) ∫_k^∞ k′P(k′) dk′.
- A = ∫ P_1D/P_1D^CDM dk en 0.5–20 h/Mpc, con P(k) en z = 3.
- δA = 1 − A/A_CDM.

**Referencias, con el mismo código:**
- WDM térmico puro de 5.3 keV (Iršič et al. 2017, 2σ);
- las mezclas en el borde de Garcia-Gallego et al. 2025 (arXiv:2504.06367): f_max = 0.14 (m/keV)^1.1 para m = 1, 2, 3 y 4 keV.

**Puerta LYA (sustituye a P5):**
- **LYA_estricta:** δA ≤ mín(δA de las referencias).
- **LYA_amplia:** δA ≤ máx(δA de las referencias).
- Se informan las dos. El veredicto usa la estricta, y la amplia se da como sensibilidad.

**Aplicación:** a todas las filas de S6 con K_sw ∈ {0, 1} que pasan P1–P4, para las tres ramas y los tres valores de R (−1σ, central, +1σ). El veredicto sigue ESPEC_U2 §6, con P5 sustituida por LYA_estricta.

**Fuera de alcance (declarado):**
- No es CLASS: no hay acoplamiento fotón–barión, ni recombinación, ni lentes. El cociente frente a CDM es lo que cuenta, y lo validan W1 y W2.
- δA es un criterio aproximado para formas no térmicas. Se dejará un job de WSL con CLASS (`classy`) para confirmarlo cuando la red lo permita.

## Enmienda 1 (6 oct, ≈18:22; tras el piloto, que terminó a las 18:20, y antes de lanzar F)

**Piloto** (N = 128, τ_end = 64, w = 3; en `net/results/runs/*N128*`).
- **El código funciona:** ξ es continuo en τ_c y la fase gruesa detecta cuerdas.
- **La preparación no elimina el colapso.** Durante la fase gruesa, ξ sube hasta ≈ 3.7 (campo desordenado: los núcleos se solapan y la fracción sin máscara es ≈ 0). Después cae a ≈ 0.4–0.55 en ℓ ≈ 3.1–3.3, igual que en U2.
- En U2 el colapso ocurría en el mismo ℓ con cualquier condición inicial (E25, A, E100). Eso sugiere que lo controla el cociente horizonte/núcleo (e^ℓ ≈ 20), no la historia previa.

**Decisión: se lanza F sin cambiar w ni τ_i.** El piloto no da motivo técnico para cambiarlos; cambiar la preparación a la vista del piloto sería un ajuste a posteriori.

**Expectativa declarada:** probablemente G1 vuelva a fallar. Si falla, la conclusión es que con N = 256 (ℓ_max = 4.85) no se puede llegar al atractor por ningún método de preparación de este tipo, y hace falta más rango dinámico (N ≥ 1024, ℓ ≳ 6.2).

## Enmienda 2 (6 oct, ≈18:37; tras las validaciones de U3a y antes de aplicar la puerta a las filas)

**W1 falla tal como estaba preregistrado.** El k_{1/2} propio es mayor que el del ajuste de Viel et al. 2005: un 12 % (2 keV), un 14 % (3.5 keV) y un 17 % (5.3 keV). W2 pasa (10⁻¹¹), y ℓ_max, N_q, el paso de tiempo y el umbral de RSA no cambian T, con efecto < 1 %.

**Ancla posterior W1′.** Se compara con la relación δA(m) obtenida con CLASS por DeRocco et al. (arXiv:2112.09105, ec. 3.2):

δA_WDM = [1 + 0.19 (m/keV)^1.63]^−0.75

| m | δA propio | CLASS |
|---|---|---|
| 2 keV | 0.677 | 0.707 |
| 3.5 keV | 0.483 | 0.508 |
| 5.3 keV | 0.355 | 0.362 |

La concordancia es del 2–5 %. El mismo artículo documenta que los ajustes tipo Viel tienen errores de este orden.

**Decisión.**
- Se informan los dos estados: «validado según el preregistro» = NO y «validado con W1′» = SÍ.
- La puerta LYA se calcula igualmente, con la etiqueta correspondiente.

**Advertencia adicional (del mismo artículo).** δA es poco fiable para distribuciones multimodales, y las nuestras lo son: frío + tibio. Ya se ve en las referencias propias: las mezclas en el borde de Garcia-Gallego dan δA = 0.23–0.32, menos que el WDM puro de 5.3 keV (0.355). Por eso el umbral estricto (0.232) y el amplio (0.355) difieren mucho. El veredicto se da con los dos.

## Enmienda 3 (6 oct, ≈18:54; decisión de cómputo, tomada tras ver solo las redes U(1) de F)

**Hechos.**
- Las tres redes U(1) de F (F1: semillas 41 y 42; F2: semilla 43) han terminado.
- ξ_U1(4.5) vale 0.522 y 0.352 en F1 (media 0.44) y 0.60 en F2. Así **G1 falla** (0.44 < 0.5), y **G2 falla para U(1)**: la memoria de la condición inicial es del 31 %, por encima del 20 %.
- Las dos puertas dependen solo del control U(1). Por tanto el resultado de U3b, GATES_FAILED, ya no puede cambiar.

**Decisión.** Las redes RP³ de F (unas 3 h con CPU compartida) se detienen para liberar CPU para U3a. La que estaba en curso (F2 RP³, en τ ≈ 27) se archiva incompleta. Quedan preparadas para `make network-u3` en una máquina con más núcleos. Su único papel habría sido informativo: la robustez de R_Γ frente a la preparación.

## Enmienda 4 (6 oct, ≈20:15; tras la auditoría independiente `audit/u3/AUDIT_U3.md`, antes de la versión final de U3a)

1. **Convergencia en k_max.** Con k_max = 120 h/Mpc, el δA de las filas tibias está infravalorado entre un 20 % y un 31 %, porque T² se congela más allá de la rejilla. Converge en ~2–3 % entre 500 y 1000 h/Mpc (auditoría, a5).
   - La producción pasa a k_max = 1000 h/Mpc, con 46 nodos.
   - Solo se recalculan las filas con δA(k_max = 120) < 0.5. Las demás ya están lejos de los dos umbrales, y su δA solo puede crecer con k_max; se conservan marcadas como cota inferior.
2. **W1′ se retira.** arXiv:2112.09105 es de Dienes et al. (no de DeRocco). Su ec. 3.2 se construye con el ajuste de Viel (CLASS solo aporta P_CDM), así que no es un ancla independiente. **U3a queda NO VALIDADA según el preregistro**: W1 falla, porque el δA propio está un 4–10 % por debajo del basado en Viel. Los umbrales se calculan con el mismo código; según la auditoría (a4), las simplificaciones del resolvedor mueven el cociente fila/umbral ≤ 0.3 %.
3. **Etiqueta del veredicto.** ESPEC_U2 §6 fuerza INCONCLUSO cuando falla el control U(1) de la red, y falla (U2, U3b). Por eso el **veredicto preregistrado de U3a es INCONCLUSO**. «CIERRE_SEGÚN_RAMA» se informa como lectura **condicional a la red**.
4. **W3 se rehace sobre las filas decisivas** (las de δA más bajo y más cercano al umbral con K_sw = 0): nodos, ℓ_max, paso y k_max = 1500. La tolerancia pasa del 2 % al 3 %, porque la W3 preregistrada se había evaluado sobre filas saturadas (δA ≈ 0.99) y no comprobaba nada.
