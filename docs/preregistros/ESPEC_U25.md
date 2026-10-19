# U25 · ¿Cuánto depende R_Γ del paso de tiempo cerca del régimen de producción? · Preregistro

*18 de octubre de 2026, 14:45 (Madrid). Escrito antes de producir ninguna red de U25 y antes de escribir su código.*

## Motivo

La validación V1 de U2 (energía en Minkowski, deriva < 10⁻³) falló tal como estaba preregistrada.
- **Lo que se midió:** un desfase constante de la energía sombra del leapfrog, del 2.5–2.8 %, que escala exactamente como dτ², con una deriva secular por debajo de 10⁻³.
- **La única prueba del efecto sobre los observables** es el conjunto D de U2: 128³, una semilla por modelo, sin drenaje. Allí R_Γ cambia 0.016 al pasar de dτ = 0.25 a dτ = 0.125.
- **Toda la producción usa dτ = 0.25.** La revisión de Codex del 18 de octubre (`paper/prd_network/PENDIENTE_ENVIO_PRD.md` §4) pide una comparación más cercana al régimen de las conclusiones.

**Pregunta:** a 384³, con el protocolo drenado de U11 y en la ventana de U24, ¿cuánto cambian R_Γ y R_ξ al reducir el paso a la mitad?

## Producción (`network/u25_produce.py`)

- **Es el protocolo `D0` de U24** (384³, τ_fin = 192, sfac = 1, ξ_i = 0.5, ℓ_d = 3.8, Γ_d = m_r, 45 instantáneas), con **dτ = 0.125**. RP³ y U(1), semillas 61, 62 y 63.
- **La pareja con dτ = 0.25** son las redes `D0` de U24, las mismas, no otras. Esto queda dicho aquí, antes de verlas.
- **Las mismas condiciones iniciales y las mismas instantáneas.** τ₀ y los τ de las instantáneas se redondean a la rejilla de 0.25 en los dos casos; esa rejilla también lo es de 0.125.
  - `lattice.run` admite un argumento nuevo, `tau_grid` (por defecto `None`, que significa dτ, como hasta ahora), que solo afecta a esos redondeos.
  - Con `tau_grid=None` la evolución tiene que ser idéntica bit a bit a la actual; se comprueba como en U24.
  - El drenaje se aplica cada 4 pasos con el factor acumulado, así que su ritmo físico no depende de dτ (`kern.c`).
- **Coste:** unos 64 min por red RP³ y 16 min por red U(1); unas **4 h** en total.
  - Se lanza **después** de las semillas 61 de U24, para que la pareja exista.
  - Salida: `network/results/u25/runs/<MODELO>_N384_s<semilla>_dt0125.json`.

## Análisis (`network/u25_analyze.py`, fijado ahora)

- **La misma ventana W = [4.36, 5.21] y los mismos estimadores que U24.**
- **Cambios:** Δ_X = X(dτ = 0.125) / X(dτ = 0.25) − 1, para X = R_Γ, R_ξ, ε_RP³, ε_U(1), ξ_RP³ y ξ_U(1), con bootstrap emparejado por semilla.
- **Extrapolación de Richardson, informativa:** si el error es O(dτ²), el sesgo de dτ = 0.25 respecto al continuo es ≈ (4/3)·|Δ_X|.

**Veredicto:**
- **PASO_ROBUSTO:** |Δ_RΓ| ≤ 0.02 y |Δ_Rξ| ≤ 0.01, las dos con σ ≤ 0.015. El paso de tiempo queda acotado en el error de R_Γ por debajo del 2 % (≈ 0.03), y así se dice en el texto.
- **PASO_SENSIBLE:** |Δ_RΓ| > 0.03 con |Δ_RΓ| > 2σ, o bien |Δ_Rξ| > 0.015 con |Δ_Rξ| > 2σ. (4/3)·|Δ| entra como sistemática del paso de tiempo en el error de R_Γ (o de R_ξ) y se declara en el texto.
- **INDETERMINADO:** cualquier otro caso. Incluye 0.02 < |Δ_RΓ| ≤ 0.03, que se informa con su valor sin añadirlo como sistemática.

**Informativo:**
- si el cambio está en ε_RP³, en ε_U(1) o en los dos;
- si se concentra al principio de la ventana, cerca del fin del drenaje.

## Lo que U25 no resuelve

- **No es 576³.** Si a 384³ y 576³ el paso pesa distinto, U25 no lo ve; con m_r·a·dx = 1 al final en los dos, se espera que pese igual.
- **No repite V1 en Minkowski.** Mide el efecto sobre los observables, que es lo que importa para el artículo.
- **No prueba dτ = 0.0625**, por coste: unas 8 h solo las redes RP³.

## Enmienda 1 · 18 de octubre de 2026, ≈16:50 (Madrid): el drenaje es de primer orden en dτ; extrapolación y umbrales corregidos, antes de cualquier red de U25

**La revisión de Codex en la PR #10** (comentario sobre las líneas 30–33) tiene razón en los dos puntos. Se comprueba en `network/kern.c`:

1. **El drenaje es de primer orden.** `damp_massive` se aplica cada 4 pasos con el factor `exp(−Γ_d·a(t₁)·4dτ)`, es decir, con la a del extremo final del intervalo. Es una cuadratura de primer orden, y el drenaje entra como un operador separado en instantes discretos, cuyo último instante antes de τ_d también depende de dτ. En las redes drenadas, el error de la evolución es O(dτ) en esa parte, no O(dτ²) como el del leapfrog. **El sesgo de dτ = 0.25 hay que estimarlo con Richardson de primer orden:** b ≈ 2|Δ|, no (4/3)|Δ|.
2. **La frase de PASO_ROBUSTO era incoherente.** Con |Δ_RΓ| ≤ 0.02 y (4/3)|Δ|, el sesgo podía llegar al 2.7 %, no «por debajo del 2 %» como decía el texto.

**No se cambia el kernel.** Hacer exacta la integral del drenaje cambiaría el protocolo de producción (U11), que la referencia D0 de U24 reproduce. U25 mide justo el efecto del paso de tiempo sobre ese protocolo tal como es.

**Cambios** (antes de producir ninguna red de U25):
- **Sesgo del paso grueso:** b_X = 2·|Δ_X|, de primer orden y conservador. (4/3)·|Δ_X| se informa también, como el valor si todo fuera de segundo orden.
- **Veredicto:**
  - **PASO_ROBUSTO:** b_RΓ + 2σ_b ≤ 0.02 y b_Rξ + 2σ_b ≤ 0.01, con σ_b = 2σ(Δ). El texto dice que el paso de tiempo está acotado por debajo del 2 % en R_Γ y del 1 % en R_ξ.
  - **PASO_SENSIBLE:** |Δ_RΓ| > 0.015 con |Δ_RΓ| > 2σ, o |Δ_Rξ| > 0.0075 con |Δ_Rξ| > 2σ. b entra como sistemática del paso de tiempo en el error de R_Γ (o de R_ξ).
  - **INDETERMINADO:** cualquier otro caso. b **también** entra como sistemática, con su valor y su error: es más conservador que el texto original, que no lo añadía.
- `network/u25_analyze.py` se ajusta a esto.

## Nota · 18 de octubre de 2026, 20:53 (Madrid), antes de ver ninguna red de U25

- **Con una sola semilla, U25 no da veredicto,** igual que U24 (enmienda 1 de ESPEC_U24, punto 1).
- **Hora real del preregistro:** este ESPEC se subió en `55070bd`, a las 14:20 (Madrid). La cabecera dice 14:45 por error.
- **Coste real:** ≈ 2 h por red RP³ con dτ = 0.125 y ≈ 22 min por red U(1), no ~4 h en total.
