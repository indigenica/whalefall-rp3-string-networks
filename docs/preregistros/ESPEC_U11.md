# U11 · U6 drenado: R_Γ y q limpias con buena estadística · Preregistro

*10 de octubre de 2026, ≈19:50 (Madrid). Escrito antes de producir ninguna red de U11.*

## Motivo

U10a (`docs/informes/INFORME_U10a.md`) mostró que entre el 25 y el 60 % de la emisión que medíamos era radiación masiva de la condición inicial que se convertía en Goldstones. Al drenar ese baño, la diferencia entre protocolos casi desaparece y salen R_Γ ≈ 1.6–1.8 y q ≈ 0.6–0.7. Pero U10a se hizo a 512³ y con dos semillas por protocolo.

U11 repite exactamente U6 (las mismas semillas, el mismo N = 576 y ξ_i = 0.5), con el drenaje de U10a:
- velocidad masiva frenada a ritmo físico m_r hasta ℓ_d = 3.8, solo donde el campo está cerca del vacío;
- aplicado cada 4 pasos.

Al usar las mismas condiciones iniciales que U6, cada red drenada tiene su gemela sin drenar, y la comparación se puede hacer **semilla a semilla**.

## Producción

`u4b_produce.py 7`, con **N = 576 fijo**: si no cabe en la RAM, el script se para en lugar de bajar N. Espectros en cada instantánea.

| Red | Semillas | X |
|---|---|---|
| RP³ drenada, ξ_i = 0.5 | 61–66 | no (en RP³ es cualitativo, V5) |
| U(1) drenada, ξ_i = 0.5 | 61–64 | sí (cuantitativo, V5) |

- **Orden:** RP³ 61; U(1) 61 y 62 en pareja; RP³ 62 y 63; U(1) 63 y 64 en pareja; RP³ 64, 65 y 66.
- **Tiempo estimado:** ≈ 2.1 h por RP³ y ≈ 0.5 h por pareja de U(1), unas 14 h en total.
- **RAM:** cada RP³ necesita ≈ 40 GB. El script exige al menos 42 GB libres al empezar cada red.

## Análisis (`network/u11_analyze.py`, fijado antes de los datos)

**Ventana:** la de U6, [ℓ_fin − 1.3, ℓ_fin − 0.05], con ℓ_fin = 5.663. Los mismos estimadores que U6 (rm2.5), con bootstrap sobre semillas.

**Medidas:**
1. R_Γ^dren y R_ξ^dren: las 6 RP³ frente a las 4 U(1).
2. **R_Γ^corr:** el mismo cociente, pero con la emisión de U(1) corregida con X, Γ_U1 − ⟨X⟩. Es la corrección cuantitativa de U(1). RP³ queda sin corregir, y se dice.
3. q_RP3^dren(ℓ = 5) y su pendiente, con el ajuste de U6 (`q_fit`).
4. **Pares por semilla:** la emisión de cada red drenada dividida por la de su gemela de U6, para RP³ y para U(1). De ahí sale la fracción de baño de cada modelo con su dispersión.
5. ρ_M/ρ_a al final de la ventana, drenadas frente a U6.

**Clasificación de R_Γ^dren**, con 1σ de bootstrap, frente a los valores de U7c:
- **LIMPIA_ZONA_PERMITIDA:** R_Γ^dren + 1σ ≤ 1.69. Con 1.69 pasaban 7 de las 8 filas de B2κ.
- **LIMPIA_ZONA_EXCLUIDA:** R_Γ^dren − 1σ ≥ 2.09. Con 2.09 fallaba x₀ = 30.
- **LIMPIA_INTERMEDIA:** en otro caso.

Se informa también la misma clasificación con R_Γ^corr.

**Clasificación de q_RP3^dren(ℓ = 5)**, con 2σ de bootstrap:
- **Q_MENOR_QUE_1:** q + 2σ < 1.
- **Q_COMPATIBLE_CON_1:** |q − 1| ≤ 2σ.
- **Q_MAYOR_QUE_1:** q − 2σ > 1.

**Lectura** (U5f: q_crit = 1.025):
- con Q_MENOR_QUE_1, la emisión no deja huella en Lyα y la materia oscura sería fría, de desalineamiento, si eso se mantiene a ℓ ≈ 70;
- con Q_COMPATIBLE_CON_1, decide R_Γ.

## Lo que U11 no resuelve

- **El comportamiento a ℓ ≈ 70.**
- **La corrección de RP³ por X**, que necesita una X exacta en la red.
- **El sesgo propio del drenaje** (el rozamiento hasta ℓ = 3.8), que no se puede separar con estas redes. U10a mostró que el drenaje deja la densidad de cuerdas casi igual.
