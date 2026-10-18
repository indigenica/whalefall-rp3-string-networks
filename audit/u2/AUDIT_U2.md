# Auditoría independiente U2 · R_Γ = Γ_RP³/Γ_U(1) y cierre s6

*6 oct 2026 · auditor independiente · solo lectura fuera de `audit_u2/`. No se ejecutó `produce.py` ni `analyze.py` ni `s6_measured.main()`, y no se tocó el job de producción.*

## Cómo se hizo

Todo el código es propio y está en `audit_u2/`. Se ejecutó con `OMP_NUM_THREADS=1` y `nice`, con N ≤ 48.

| Script | Qué hace | Salida |
|---|---|---|
| `audit_kernel.py` | Comprueba la física del kernel (checks 1–3). Usa una copia compilada de `kern.c` (`_build/kern_audit.so`, idéntica byte a byte al `_build/kern.so` de producción). El potencial, el gradiente, el Hessiano, el estimador Goldstone y el detector de plaquetas están reimplementados en numpy a partir de ESPEC. | `audit_kernel.json` |
| `audit_analysis.py` | Recalcula el análisis (check 4) desde `results/runs/*.json`, sin importar `analyze.py`. | `audit_analysis.json` |
| `audit_closure.py` | Comprueba el cierre (check 6). Rehace la lógica de veredicto desde las filas de `s6_measured_closure.json`, reproduce nodos y evalúa la rama B3. Importa funciones de `s6_measured` y no guarda nada fuera de `audit_u2`. | `audit_closure.json` |
| `audit_bath.py` | Prueba exploratoria propia: cuánto Goldstone produce un baño masivo sin cuerdas. | `audit_bath.json` |

**Conjunto C.** `RP3_N256_s21_f0.6_dt0.25_xi0.5.json` apareció durante la auditoría (14:5x). Está incluido en el recómputo. `U2_RESULT.json` y `s6_measured_closure.json` son anteriores y no lo incluyen.

## Tabla de comprobaciones

| # | Comprobación | Resultado | Números |
|---|---|---|---|
| 1a | Fuerza de `kern.c` = a²(∇²q − a²∂U/∂q), con U ya dividido por Z_Q (LAM_Z = 320/400 = 0.8, LAMP_Z = 192/400 = 0.48); base de Frobenius ortonormal (fuera de la diagonal ×1/√2) | **PASS** | Error relativo máximo frente al gradiente numérico propio, en 300 sitios genéricos fuera del vacío (N = 16, a = 0.5): 2.0e-5 (RP³) y 2.1e-5 (U(1)), en el nivel de float32 |
| 1b | EOM conforme con π = a²φ' en radiación (a = sτ, a'' = 0, sin término a''/a) | **PASS** | Modos lineales frente a la ODE φ'' + (2/τ)φ' + (k²_lat + a²m²)φ = 0 resuelta con solve_ivp, τ: 4 → 60. Error/amplitud: Goldstone 1.7e-4 (ambos modelos), quinteto 4.5e-4, radial U(1) 2.6e-3, pesado 4.0e-3 (error O(dτ²) del leapfrog) |
| 1c | Espectro en Q = diag(1,0,0,0) a partir del Hessiano propio | **PASS** | Autovalores {0, 0, 0, 0.24 ×6, 1.84}: 3 Goldstone, quinteto + singlete degenerados con m² = λ_P/(2Z_Q) = 0.24, y pesado con m² = 0.24 + 2λ/Z_Q = 1.84. ∇U(vacío) = 2e-11. U(1): {0, 0.24} |
| 2a | K_G (RP³) = \|Q'n\|² − (n·Q'n)² = ½Tr(Q'_G²) = \|n'\|², y K_G (U(1)) = (φ×φ')²/(2\|φ\|²), con la misma normalización (energía cinética física; entre a⁶ en variables π, es decir \|φ'\|²/a² con φ' conforme) | **PASS** | Onda Goldstone de amplitud idéntica en unidades de f (ángulo de periodo 2π; en RP³ θ₂ = 2·ángulo de n), A = 0.05, 1.0 y 2.5: K_G^RP³/K_G^U1 = 1.0000000 y K_G/(½f²θ̇²) = 1.0000000. No hay ningún factor 2 entre modelos |
| 2b | Lo mismo de forma dinámica: paquete idéntico evolucionado en radiación (N = 32) | **PASS** | a⁴ρ_a^RP³/a⁴ρ_a^U1 = 1 ± 3e-7 en 5 instantes. Advertencia: la onda en un plano de RP³ es un subsector U(1) exacto, así que esta prueba fija la normalización pero no la dinámica S³ |
| 3a | Holonomía ℤ₂ (Π sign(n_i·n_j) = −1) y devanado U(1) | **PASS** | Detector numpy propio sobre campos aleatorios suavizados (N = 24): RP³ 1580 = 1580 y U(1) 816 = 816. El autovector principal del kernel (Jacobi) coincide con `numpy.eigh` (\|n·n'\| ≥ 0.99999996) |
| 3b | L_c = (2/3)N_plaq | **PASS** | 24 anillos de radio R = 9 con orientación aleatoria: (2/3)N_plaq/(2πR) = 1.008 ± 0.031 (igual en RP³ y U(1)) |
| 3c | ξ = L_c τ²/(4V_c) | **PASS** | Coincide con L_fís t²/V_fís (t = 1/2H, H = 1/(sτ²), a = sτ) hasta 1e-15 |
| 4a | Recómputo de R_ξ, R_Γ (local e integral) del conjunto A | **PASS** (coincide con U2_RESULT hasta 1e-15) | R_ξ = 1.562 ± 0.087. R_Γ local = 1.6945 ± 0.038 (boot); el error clásico entre semillas es 0.044. R_Γ,int = 1.6652 ± 0.037. ε_lat: RP³ 0.415, U1 0.245. η: RP³ 0.60, U1 0.54 |
| 4b | Combinación de sistemáticos (ESPEC §4 y enmiendas 1, 4 y 5) | **PASS** con matices | Se suman en cuadratura: semillas, B, D (dτ), E25, E100, r_m ∈ {1.5, 4} y C (nuevo). La diferencia entre estimadores (0.029) es menor que σ, así que no se añade (enmienda 5). Total sin C = 0.2861, que coincide. **Con C = 0.2871** (Δ_C = 0.024). El sistemático dominante es IC E100: 0.267 |
| 4c | Prueba de tendencia de R_Γ en ℓ (ESPEC §5, > 2σ → rama B3) | **FAIL (bug)** | `analyze.py:118`: `np.polyfit` sobre `per_ell` con el último punto NaN da pendiente NaN, así que `significant` es siempre False. La pendiente real es −0.48/ℓ. Con el bootstrap del propio código σ = 0.056, z = 8.6, y B3 debería haberse activado. Con la dispersión OLS σ = 0.33, z = 1.5. Ver hallazgo 1 |
| 4d | Puerta de control U(1) | **FAIL de la puerta, bien aplicada** | ξ_U1(4.5) = 0.451 (semillas 0.430, 0.433, 0.500, 0.440), por debajo de 0.5. η_U1 = 0.54 ∈ [0.3, 3]. La puerta falla y el veredicto pasa a INCONCLUSO, como está preregistrado. Es un fallo marginal: B da 0.52, E25 0.46, E100 0.54 y C 0.43 (este último en ℓ = 4.2) |
| 4e | ¿Está la ventana [3, 4.85] en el atractor? | **NO del todo** | Ver la sección «Transitorio» |
| 5 | Literatura | Ver la tabla de literatura | — |
| 6a | Fórmulas de rama en s6: ε_B1 = R_Γ y ε_B2 = R_Γ·1.20/(0.24ℓ) | **PASS** | eps_range de 8 filas igual a la fórmula sobre ell_range (rtol 1e-6). Dos nodos reproducidos exactamente (m_G, f_em, allowed) |
| 6b | Lógica de veredicto (ESPEC §6) y anulación por la puerta U(1) | **PASS** | Antes de la puerta: ninguna fila de 48 está permitida, ni siquiera al −1σ → FALSADA_EN_FORMA_ACTUAL. Fallos: P5 en 48 filas y P4 en 10. Después de la puerta: «INCONCLUSO (control U(1) fuera de rango)», correcto según §6 |
| 6c | Etiquetado de las piezas exploratorias | **MATIZADO** | Ver hallazgo 3 |
| 6d | Rama B3 que faltaba (evaluada por el auditor) | **sin efecto en el veredicto** | Ajuste R(ℓ) = −0.022 + 6.57/ℓ, que da ε ≈ 0.04–0.28 en ℓ_fís ∈ [25, 100]. Excluida en los 12 nodos con K_sw = 0 (−1σ, central, +1σ; dos brazos espectrales; x0 = 10 y 30), siempre por P5 |

## Hallazgos (bugs y desviaciones) y su efecto sobre R_Γ

1. **Bug en la prueba de tendencia** (`analyze.py:118`). `per_ell[-1]` es NaN porque el último punto de la malla, 4.80, cae más allá del último instante con derivada centrada. `np.polyfit` devuelve NaN y `abs(nan) > 2σ` es False.
   - Efecto sobre el valor de R_Γ: **ninguno**.
   - Efecto sobre la lógica: la rama B3 preregistrada nunca se evaluó.
   - ¿Es real la tendencia? Probablemente no. Las oscilaciones de `per_ell` son coherentes entre semillas: la caída de ℓ ≈ 3.25 aparece en todas las curvas, lo que apunta a una oscilación global del modo masivo heredada de la condición inicial con π = 0. Por eso el bootstrap sobre semillas subestima el ruido en ℓ.
   - El estimador integral por subventanas no muestra tendencia monótona: 1.65 en [3, 3.6], 2.29 en [3.6, 4.2] y 1.33 en [4.2, 4.8]; por mitades, 1.64 frente a 1.69.
   - El auditor evaluó B3 y queda excluida, así que **el veredicto no cambia**.
   - Recomendación: arreglar con `mk = np.isfinite(...)` y declarar qué σ se usa.
2. **Ventana efectiva del estimador local:** [3.0, 4.64], no [3.0, 4.80], porque el último punto es NaN. Efecto despreciable: el estimador integral cubre [3, 4.80] y difiere en 0.029.
3. **s6, piezas exploratorias.**
   - Están marcadas como «exploratory, not preregistered».
   - Pero la puerta mixta Lyman-α y el factor de supresión se guardan dentro de la clave `exploratory_massless_variant`, cuya nota describe solo la variante sin masa.
   - La cita «Hervás Peters et al.» es **incorrecta**: arXiv:2504.06367 es de Garcia-Gallego, Iršič, Haehnelt, Viel y Bolton. La fórmula sí es la del artículo.
   - El código extrapola fuera de 1–4 keV y pone f_max = 1 para m ≥ 5.3 keV, con un salto desde 0.88. Sin efecto en el veredicto preregistrado.
4. **V4 implementado distinto del texto de ESPEC.**
   - El texto pide a⁴ρ_a constante al 3 %. El código pide que la energía comóvil total varíe menos del 3 % (da 2.94 %, justo dentro) y que ⟨ρ_a⟩/E esté a menos del 15 %. El cambio no está en las enmiendas.
   - La prueba 2b de este informe cubre la normalización entre modelos, que es lo que importa para R_Γ.
5. **La calibración de la condición inicial (enmienda 4) no se cumple.**
   - ξ medido en el primer instante (ℓ ≈ 1.16): ≈ 2.0–2.2 para el nominal 0.5, ≈ 1.2–1.7 para 0.25 y ≈ 3.2–3.7 para 1.0.
   - Las etiquetas son nominales. Aun así, E cubre un factor de ~3 en ξ inicial, así que la prueba de dependencia con la condición inicial sigue siendo informativa.
6. **Estadística.** El conjunto A de RP³ tiene 2 semillas: el lote `A_RP3b` no se ejecutó. Con n = 2 el bootstrap es pobre, pero el error clásico (0.044) lo confirma y es despreciable frente a los sistemáticos.
7. **Sistemático no presupuestado: conversión del baño masivo en Goldstones** (prueba propia, exploratoria, orden de magnitud).
   - En ℓ = 3, a⁴ρ_M ≈ 1.1 × Δ(a⁴ρ_a) de toda la ventana, en ambos modelos. Al final, ρ_M/ρ_a = 0.98 (RP³) y 0.72 (U1).
   - Un baño gaussiano masivo sin cuerdas, con la misma energía y con semilla Goldstone, añade en la ventana ≈ 10 % del ΔG de la red en U(1) y ≈ 1 % en RP³.
   - Si se restara, **R_Γ subiría ~9 %**, es decir +0.15, del mismo orden que la mitad de σ_total. Va en contra del cierre, no a favor.
   - No es un bug. Es un sistemático que ESPEC no contempla.
8. **El conjunto C no está en `U2_RESULT.json`.** Recomputado aquí: R_Γ(C) = 1.793 frente a 1.817 de A en la misma ventana (Δ = 0.024; integral: 0.018). σ_total pasa de 0.286 a 0.287.
   - La resolución del núcleo sí cambia R_ξ: 1.75 (C) frente a 1.52 (A) en [3, 4.29]. Viene de que ξ_RP³ sube un 16 % con m·a·dx ≤ 0.6.
   - R_Γ es robusto a la resolución; R_ξ no.

**No se encontró ningún error de factor en kernel, EOM, normalización Goldstone, detector ni ξ.** R_Γ = 1.69 ± 0.29 queda confirmado tal como se definió. El único bug cambia la lógica de ramas, no el número, y una vez corregido no cambia el veredicto.

## Transitorio y atractor

- **Transitorio.** La red parte con ξ ≈ 2 en ℓ = 1.2, con el campo proyectado al vacío cuando H ≈ m/e y núcleos de ~6 celdas. Sube a ≈ 2.2–2.3 (U1 en ℓ ≈ 2.1, RP³ en ℓ ≈ 1.3–2.3) y se desploma hasta el mínimo de ℓ ≈ 3.25: 0.36 (U1) y 0.50 (RP³).
- **Lectura del colapso.** Es la aniquilación de la red sobredensa que se forma cuando el término de masa todavía no domina, junto con la conversión de esa energía en radiación masiva: en ℓ = 3, a⁴ρ_M/a⁴ρ_a = 2.5 (RP³) y 1.2 (U1).
- **Con qué encaja.** GHV describen exactamente este comportamiento («initially overdense simulations show first a drop and then a universal increase»). Aquí la caída **sobrepasa** y deja a U(1) por debajo de la literatura.
- **La ventana [3, 4.85] empieza justo en el mínimo**, no sobre el atractor:
  - dξ/dℓ en la ventana: U1 0.065–0.10 y RP³ 0.16–0.20 (E100 RP³ 0.01, porque todavía está bajando: ξ = 2.03 en ℓ = 3). En la literatura c₁ ≈ 0.24–0.25.
  - Memoria de la condición inicial en ℓ = 4.5: U1 0.45/0.46/0.54 y RP³ 0.74/0.66/0.77 (A/E25/E100), un ~20 % sin converger.
  - R_ξ por ℓ sube de 1.46 a 1.66 dentro de la ventana.
- **Qué significa para R_Γ.** Es un cociente medido sobre una red en aproximación, no en escalado. Su sensibilidad a la condición inicial (Δ = 0.27, el término dominante) ya está en σ. Que el control U(1) quede por debajo del rango de la literatura es justo lo que detecta la puerta: la puerta funciona como debe.

## Literatura (verificada en la web)

| Afirmación | Estado | Nota |
|---|---|---|
| GHV arXiv:2007.04990: ξ ≈ c₁ℓ con c₁ ≈ 0.24 | **CONFIRMADO** | Ec. (3): ξ = c₁ log + c₀ + c₋₁/log + c₋₂/log², con c₁ = 0.24(2). Escalado a partir de log ≈ 4. Extrapolación ξ(70) = 15(2). B1 usa solo el término 0.24ℓ (ignora c₀): aceptable para ℓ_fís ≳ 25 |
| Hindmarsh et al. arXiv:2102.07723: ζ* = 1.20 ± 0.09, escalado constante | **CONFIRMADO** | Abstract: «ζ* = 1.20 ± 0.09, v* = 0.609 ± 0.014». Punto fijo consistente entre cuerdas físicas y de anchura comóvil. Atribuyen el crecimiento logarítmico a una aproximación incompleta desde ζ bajo. La longitud se mide en el sistema de reposo cósmico, igual que ξ aquí |
| ξ típico de cuerdas U(1) físicas en ℓ ≈ 4–5 | **MATIZADO** | La literatura da O(1). Buschmann et al. (arXiv:2108.05368) ajustan c₁ = 0.254 ± 0.002 en log ∈ (4, 9). No se ha verificado un número exacto en ℓ = 4.5; el rango [0.5, 2.5] de la puerta es amplio y razonable. El 0.45 medido está por debajo en las 4 semillas de A |
| arXiv:2504.06367: f_WDM < 0.14 (m/keV)^1.1 | **MATIZADO** | La fórmula, el 2σ y el rango de 1–4 keV son correctos (f_WDM < 0.16 / 0.35 / 0.50 / 0.67 para 1 / 2 / 3 / 4 keV). La autoría citada («Hervás Peters et al.») es **INCORRECTA**: son Garcia-Gallego, Iršič, Haehnelt, Viel y Bolton. Solo se usa en una pieza exploratoria |

Fuentes: [arXiv:2007.04990](https://arxiv.org/pdf/2007.04990), [arXiv:2102.07723 (ar5iv)](https://ar5iv.labs.arxiv.org/html/2102.07723), [arXiv:2108.05368 (ar5iv)](https://ar5iv.labs.arxiv.org/html/2108.05368), [arXiv:2504.06367](https://arxiv.org/abs/2504.06367).

## Distancia al borde (contexto)

Según el factor de supresión exploratorio de s6, para que una rama pase con K_sw = 0 haría falta:
- B1: R_Γ ≲ 0.008–0.03.
- B2: R_Γ ≲ 0.12–0.48.

El valor medido, 1.69 ± 0.29, está a más de 4σ de cualquiera de los dos. Ni un error de factor 2 en el estimador ni el sistemático del baño masivo (+9 %) acercarían el resultado al borde. Lo único que podría cambiarlo es la extrapolación de la literatura (ε_U1 en el ℓ físico), que es justo lo que las ramas parametrizan.

## Veredicto del auditor

**R_Γ = 1.69 ± 0.29 queda confirmado** (con C: ± 0.29, Δ_C = 0.024). No hay ningún error de normalización entre modelos.

**Veredicto preregistrado: INCONCLUSO, correctamente aplicado.** La puerta U(1) falla: ξ_U1(4.5) = 0.45 < 0.5. Sin esa puerta, el resultado sería FALSADA_EN_FORMA_ACTUAL, con B1, B2 y también B3 excluidas incluso al −1σ.

Un bug (la prueba de tendencia, que da NaN) suprimió la rama B3. Corregido, no cambia nada. Una vez etiquetadas como exploratorias, las piezas que no se preregistraron no afectan al veredicto.
