# Auditoría independiente U4–U11 · Plan y criterios

*11 de octubre de 2026, 20:08 (Madrid). Escrito antes de escribir ni ejecutar ningún código de la auditoría.*

## Reglas

- El auditor (Claude, sesión en la nube) **no lee ni importa** `network/*_analyze.py`, `network/u9_extra.py`, `network/u10a_extra.py` ni `lya/*.py`. Reimplementa los estimadores a partir de las definiciones de los preregistros (ESPEC_U2 §3 y §6.5, ESPEC_U6, ESPEC_U9, ESPEC_U10a, ESPEC_U10b, ESPEC_U11) y de los datos versionados en `network/results/u4/runs/`.
- Las cifras de referencia son las de los JSON de resultados (`U6_RESULT`, `U9_RESULT`, `U10A_RESULT`, `U10B_RESULT`, `U11_RESULT`) y las de los informes.
- Solo se escribe en `audit/u4_u11/`.

## Qué se audita

1. **R_Γ y R_ξ** en la ventana [ℓ_fin − 1.3, ℓ_fin − 0.05] para U6, U9 (densas), U10a (drenadas, densas y poco densas), U10b (redes lanzadas tarde) y U11, con bootstrap sobre semillas.
   - Γ = H a⁻⁴ d(a⁴ρ_a)/d ln τ con la máscara rm2.5, diferencias centradas suavizadas sobre tres instantes; también el estimador integral Δ(a⁴ρ_a).
2. **Drenaje:** emisión drenada sobre emisión de la gemela de U6, semilla a semilla (U11).
3. **Intercambio X en U(1):** f_X = ⟨X⟩/Γ en la ventana y R_Γ^corr (U11).
4. **q_RP3(ℓ):** pendiente de la emisión instantánea F(k) ∝ k^(−q) entre instantáneas consecutivas en k_fís ∈ [10H, m/4] con ≥ 8 bins; ajuste q0 + q1(ℓ − 5) para ℓ ≥ 4.6, bootstrap sobre semillas (U6, U9, U10a, U11).
5. **Pendientes dξ/dℓ** y la lectura del atractor de U10b, a título de comprobación.
6. **Veredictos**: se recalcula cada clasificación preregistrada con los números auditados.

**No se audita aquí** la puerta Lyα de U7c: necesita CLASS, que no se pudo instalar en esta sesión. Queda pendiente para una sesión con CLASS o para WSL.

## Criterios

Para cada número con σ de bootstrap publicada:
- **PASA:** |auditado − publicado| ≤ 0.5 σ.
- **MATIZADO:** entre 0.5 σ y 1 σ.
- **DISCREPANCIA:** más de 1 σ.

Para números sin σ (cocientes por semilla, ρ_M/ρ_a, f_X): PASA si la diferencia relativa es ≤ 3 %, MATIZADO hasta 10 %, DISCREPANCIA por encima.

**Veredictos:** FAIL si un veredicto preregistrado cambia con los números auditados.

**Sensibilidad a decisiones de implementación** que los preregistros no fijan del todo (esquema de derivada, interpolación en ℓ, cómo se asigna ℓ a un par de instantáneas): se ensayan al menos dos variantes razonables. Si su dispersión supera la σ publicada, se informa como sistemática no contabilizada.

## Addenda 1 · 11 de octubre de 2026, 23:16 (Madrid): la puerta δA de U7c sobre dos filas

Ahora que CLASS funciona en la sesión, se audita la puerta δA de U7c sobre dos filas decisivas de B2κ con x₀ = 30: R_Γ = 2.09, K_sw = 0 (publicada con δA = 0.253, fuera del umbral) y R_Γ = 1.69, K_sw = 0 (0.219, dentro).

**Qué se reimplementa.**
- El estimador δA (ESPEC_U3: P_1D = (1/2π)∫_k^∞ k′P dk′; A = ∫ P_1D/P_1D^CDM dk en 0.5–20 h/Mpc, a z = 3), con cuadratura propia: trapecios en ln k sobre P(k) muestreado por el auditor, y cola más allá de k_max extrapolada como ley de potencias.
- Las referencias (WDM de 5.3 keV y las mezclas de 1–4 keV de Garcia-Gallego) y el umbral estricto (el mínimo).
- κ(x₀ = 30) de U7a, recalculado con una integral propia del número de Goldstones emitidos.

**Qué no se reimplementa**, y se dice: la distribución en momento de la fila (`chain/s6_measured.source_eps` y `lya/u4a_class.psd_from_packets`, auditadas en U1–U3) y los parámetros de CLASS (`thermal_params`, `packet_params`). CLASS se ejecuta aquí, sin caché.

**Criterio.** Para cada fila y referencia: PASA si |δA_auditado − δA_publicado| ≤ 3 % relativo; MATIZADO hasta 10 %; DISCREPANCIA por encima. FAIL si alguna de las dos filas cambia de lado del umbral estricto.
