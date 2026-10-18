# Revisión pre-lanzamiento de U4 (U4a CLASS + U4b redes grandes)

*6 oct 2026, antes de lanzar en la máquina WSL. Solo lectura del código: no se ha tocado ningún fichero fuera de `audit/u4_prelaunch/`.*

**Límites de esta revisión.** CLASS/classy no se puede instalar en este sandbox (PyPI y GitHub bloqueados), así que la API de CLASS se ha comprobado con su documentación:
- `explanatory.ini` (copia del repositorio en GitHub, la parte visible);
- el artículo CLASS IV, Lesgourgues y Tram 2011, arXiv:1104.2935;
- lo que recuerdo del código fuente (`input.c`, `parser.c`, `fourier.c`, `python/setup.py`).

Lo que sale solo de memoria va marcado como *(memoria)*.

**Pruebas independientes:**
- `test_slab_edt.py` (log en `test_slab_edt.log`). Compara la transformada de distancia por placas con la distancia periódica exacta (cKDTree con `boxsize = N`). Pasa: ALL OK.
- Inspección de los ficheros psd y de la caché del mock de U4a.
- Tiempo real de la red U(1) con N = 384 que se está ejecutando aquí.

---

## Resumen

| # | Ítem | Veredicto | Arreglo concreto |
|---|---|---|---|
| 1 | `env/check_env.py`: la prueba ncdm-fichero llama a `c.pk(k, 3.0)` sin `z_pk`. El `z_max_pk` por defecto de CLASS es 0, así que `pk` lanza `CosmoSevereError` ("increase input parameter z_max_pk"). `check_env` termina en FALLA y **J091 sale con código 3 sin ejecutar U4a**. | **BUG** | Línea 58, en el diccionario `base`: añadir `"z_pk": 3.0,` (p. ej. `..."output": "mPk", "z_pk": 3.0,`). |
| 2 | `lya/u4a_class.py:168`: `wide = max(ref_dA.values())` incluye las referencias `table_*` (de Garcia-Gallego). El preregistro dice "estricto (mín.) y amplio (máx.) como en U3a", y U3a usa el máximo sobre 5.3 keV + mezclas 0.14 m^1.1. | **BUG** (desviación del preregistro; afecta solo a los informes "amplio") | `wide = max(v for k, v in ref_dA.items() if not k.startswith("table"))`. Añadir `out["deltaA_table_max"] = max(v for k, v in ref_dA.items() if k.startswith("table"))` como dato informativo. |
| 3 | `Makefile:4` `export OMP_NUM_THREADS ?= 2`. Con los comandos de la guía (`make u4b-phase1/2` lanzados a mano) el kernel corre con **2 hilos**, de 3 a 4 veces más lento que el plan, porque el plan usa `os.cpu_count()`. En `make u4a`, cada worker de CLASS hereda 2 hilos: se sobresuscribe la CPU ×2, y `setdefault` no lo corrige. | **BUG** (rendimiento) | (a) En las recetas `u4b-phase1` y `u4b-phase2`: `cd network && OMP_NUM_THREADS=$${U4_THREADS:-$$(nproc)} $(PY) -B u4b_produce.py ...`. No se cambia el valor global, que usa `make test`. (b) `lya/u4a_class.py:127`: sustituir `os.environ.setdefault("OMP_NUM_THREADS", "1")` por `os.environ["OMP_NUM_THREADS"] = "1"`. Funciona porque en el hijo (spawn) classy se importa después, dentro de `class_pk`. |
| 4 | `network/u4b_plan.py:32`: con `N_force` (`make ... N=512`, que la guía recomienda si falta memoria) **no se comprueba la memoria**. La fase 2 lanzaría RP³ a N = 512 (28 GB) aunque no quepa: en WSL eso es OOM o swap durante horas. Además contradice "RP³ al N mayor que quepa". | **BUG** | `pick = lambda model: next((N for N in CANDIDATES if mem_GB(model, N) <= budget and (N_force is None or N <= N_force)), None)`. Así `N_force` es un tope, no una orden. Si se quiere forzar de verdad, añadir una opción `--N-exact`. |
| 5 | CLASS: `N_ncdm=1`, `m_ncdm`, `T_ncdm`, `omega_ncdm` a la vez | OK | `explanatory.ini`: "if both (omega_ncdm, m_ncdm) are passed, CLASS will renormalise the psd in order to fulfill both conditions"; CLASS IV dice lo mismo. La normalización de f0 en el fichero da igual; m y T fijan la forma. |
| 6 | `T_ncdm` en unidades de T_γ; `q` del fichero en unidades de T_ncdm,0 | OK | `explanatory.ini` ("in units of photon temperature") y CLASS IV (q ≡ p/T_ncdm,0, ε = (q² + a² m²/T²)^½). |
| 7 | Formato del fichero psd | OK | 2 columnas, q y f0(q). CLASS interpola con un spline cúbico, toma f(q<q_min) = f(q_min) y usa la cola α e^(−βq) por encima de q_max (CLASS IV). `np.savetxt` no escribe cabecera, lo que importa porque la lectura es con `fscanf("%lf %lf")` y una cabecera la cortaría *(memoria)*. Los ficheros del mock: q monótona, f0 > 0 (mínimo ~1e-17), sin ceros. |
| 8 | `ncdm_fluid_approximation = 3` | OK (significado) / **RIESGO** (tiempo) | 3 = `ncdmfa_none` (enum mb = 0, hu = 1, CLASS = 2, none = 3) *(memoria)*: jerarquía completa hasta z = 3 para k ≤ 674 Mpc⁻¹. Es preciso pero caro, sobre todo porque las distribuciones de las filas son anchas (q de 0.09 a 1700, con el 99 % por debajo de 80): CLASS necesitará muchos nodos de q, y `tol_ncdm_synchronous = 1e-4` es 10 veces más estricto que el valor por defecto. Las "2–4 h" no están medidas. Hacer **primero `make u4a-quick`** y mirar `meta.seconds` en `cache/*.npz`. Si una fila tarda más de 20–30 min, plantear una enmienda (fluid approximation 2 con C1 como control) **antes** de la corrida completa. |
| 9 | `l_max_ncdm`, `tol_ncdm_synchronous`, `tol_ncdm_bg` | OK | Son nombres de parámetros de precisión válidos de CLASS v3. |
| 10 | `'P_k_max_h/Mpc': 1000.0`, `'z_pk': 3.0` (float) | OK | classy convierte los valores con `str()`. La lista de k de las perturbaciones termina en un valor ≥ k_max, así que `pk(1000·h)` está dentro del rango *(memoria)*. Margen barato: `P_k_max_h/Mpc = 1.05*KMAX`. Hacerlo ahora, porque si no la caché queda partida. |
| 11 | `c.pk(k·h, z)`: k en Mpc⁻¹ y P en Mpc³, lineal (no hay `non linear`) | OK | `deltaA` usa KH en h/Mpc y P en Mpc³. El factor h³ es constante y se cancela en P1D/P1D_CDM; la pendiente de la cola sale de la CDM. Correcto. |
| 12 | `get_background()['(.)rho_ncdm[0]']`, `['(.)rho_crit']` | OK | Las claves existen en v3 y la llamada está dentro de `try`. **Mejora:** comprobar `abs(Omega_ncdm_0*h² − om)/om < 1e-3` en cada ejecución, para cazar una renormalización ignorada. Hoy se guarda y nadie lo mira. |
| 13 | `struct_cleanup()` y `empty()` | OK | Existen en classy 3.x. |
| 14 | `omega_cdm = 0` (f_w = 1, y C0 de check_env) | OK *(memoria)* | CLASS v3 pone `Omega0_cdm_min_synchronous` (1e-10) en gauge síncrono. check_env lo prueba una vez arreglado el BUG 1. |
| 15 | Fila con error de CLASS (p. ej. muestreo en q que no converge) | **RIESGO alto** | Un `CosmoComputationError` en un worker tumba todo `pool.map` y la corrida entera, que al reanudar vuelve a fallar. En `class_pk`, envolver `c.compute()` en `try/except Exception as e:` y, si falla, devolver `np.full_like(KH, np.nan)` y escribir `meta["error"]=str(e)` (sin cachear). En `main`, marcar la fila con `error` y excluirla del veredicto. |
| 16 | Escritura de la caché no atómica (`np.savez`) | **RIESGO** | Si se mata el proceso a mitad, queda un `.npz` corrupto y la reanudación falla en `np.load`. Arreglo: `tmp = f.with_name(f.stem + ".tmp.npz"); np.savez(tmp, ...); os.replace(tmp, f)`. |
| 17 | Clave de caché: parámetros + `classy.__version__` + `classy.__file__` | OK, con un matiz | Separa bien el mock (`__version__ = "MOCK"`, otra ruta). No incluye el **contenido** del psd, solo su nombre. En la práctica T_ncdm = xmed·M/k_BT lo cambia casi siempre, pero si se rehace U4a con otro R_Γ conviene meter el sha1 del fichero en la clave: `tag["_psd_sha"] = hashlib.sha1(open(p["ncdm_psd_filenames"],"rb").read()).hexdigest()` si existe. |
| 18 | multiprocessing con spawn | OK | `_job` es de nivel de módulo; el `sys.path.insert` de arriba se ejecuta en los hijos (`__mp_main__`); no hay claves de caché duplicadas en un mismo `pmap`. Mejora: `pool.imap(_job, ..., chunksize=1)`, porque los tiempos por ejecución son muy desiguales. |
| 19 | Mapeo físico de la masa sustituta | OK | k_B T_ncdm = xmed·M_surr y q_file = x/xmed, así que p = q T = x·M_surr y M/T = 1/xmed. Entonces v = q/√(q² + a²/xmed²) = x/√(x² + a²) exactamente. La densidad relativista temprana n·⟨p⟩ es invariante, porque n ∝ 1/M_surr y ⟨p⟩ ∝ M_surr. deg_ncdm queda absorbido por la renormalización. Los T_ncdm reales (de 0.014 a 0.29) dan M = m/T entre 1.5e7 y 3e8, el mismo orden que un WDM térmico de 1–5 keV, que CLASS integra sin problema. |
| 20 | Convención térmica `(T_x/T_ν)³ = ω·94 eV/m` y `T_ncdm = (…)^(1/3)·(4/11)^(1/3)` | OK | Es la de Viel et al. 2005, con T_ncdm en unidades de T_γ. Con m y ω pasados, CLASS solo usa T para la forma (m/T), que es lo que fija la convención. C0: T = 0.1054 para 3.5 keV, coherente con check_env. |
| 21 | `fit_shape` (RegularGridInterpolator con dimensión de salida, Nelder-Mead con recorte) | OK / RIESGO bajo | Funciona; la malla en f entre 0.75 y 1 es gruesa para ln T². Las filas decisivas (energy, f_em de 0.1 a 0.5) están bien muestreadas. |
| 22 | `f_max_strict` entre 4 y 5.3 keV, y `f_max_wide` por debajo de 1 keV | **RIESGO** (no está en el preregistro) | El preregistro no define el tramo de 4 a 5.3 keV; el código interpola en log. Declararlo en una enmienda antes de ver resultados, o informar aparte de las filas con m_fit ∈ (4, 5.3). |
| 23 | C1: "doble de puntos de momento" | RIESGO bajo | `refine=2` duplica la malla del fichero, no los nodos de CLASS; esos los fija `tol_ncdm_*`, que ya baja 10 veces. Aceptable. Si se quiere literal: `ncdm_quadrature_strategy=3` con `ncdm_N_momentum_bins` el doble. |
| 24 | C1 relativo con δA pequeño | OK | El δA mínimo de U3a es ≈ 0.09, así que el 2 % son 0.0018 en valor absoluto, alcanzable. |
| 25 | C0 que falla no detiene la corrida | RIESGO (tiempo) | Si `not out["C0"]["passed"]`, escribir el JSON parcial y salir (código 4) salvo `--force`. Ahorra horas. |
| 26 | Lógica del veredicto (ESPEC_U2 §6 con P5 → LYA_primary) | OK | `CIERRE_FISICO`, `FALSADA`, `SEGUN_RAMA`, `INCONCLUSO`, sobre filas con K_sw = 0 y +1σ/−1σ, idéntico a U3a. La puerta de red sale de S6 o de U4B (`blocks[-1]`). |
| 27 | U4b: conjuntos de producción frente a ESPEC_U4 | OK | Fase 1: ξ_i = 0.25 y 2.0 (s61) y 0.5 (s61, s62). Fase 2: RP³ a ξ = 0.5 (s61, s62) más U(1) al mismo N (se salta si ya existe). |
| 28 | Puertas G1′, G2′, G3′ y la ventana de R_Γ | OK | ℓ_ref = ℓ_fin − 0.3; G1′ ∈ [0.5, 2.5]; G2′ = \|Δξ\|/media < 20 %; G3′ en [ℓ_ref − 0.7, ℓ_ref]; R_Γ en [ℓ_ref − 1, ℓ_fin − 0.05] con bootstrap y los sistemáticos de máscara y dτ de U2 (las claves existen en U2_RESULT.json). |
| 29 | Fase 2 a N₂ ≠ N₁ | **RIESGO** | El bloque N₂ no tiene corridas de memoria inicial, así que G2′ da NaN, el estado queda GATES_FAILED y R_Γ se publica en un bloque "fallido". `blocks[-1]` sigue siendo N₁ (está bien para J093 y U4a). Arreglo: en `block()`, si no hay U[0.25] ni U[2.0], `status = "R_GAMMA_ONLY"` y `gates_from = N` del bloque más grande con puertas. |
| 30 | El N cambia entre reinicios | **RIESGO** | `plan()` usa el `MemAvailable` del momento: al relanzar con menos memoria libre puede elegir otro N y empezar otra serie. El bloque mayor queda incompleto (G2′ NaN). Arreglo: guardar `results/u4/plan_phase{1,2}.json` la primera vez y reutilizar su N si existe. |
| 31 | Modelo de tiempo (`P_TIME = 3.2`) | **RIESGO** | Los pasos van ∝ τ_end = N/2, así que la evolución escala como **N⁴**. Medido aquí: U(1) a N = 384 con 2 núcleos tarda 1331 s frente a los 910 s del modelo, con algo de contención por el mock de U4a. Extrapolando con N⁴ y 8 núcleos: N = 768 ≈ 2 h por red (el plan dice 0.8) y N = 1024 ≈ 6 h por red (el plan dice 1.9). La fase 1 a 1024 llevaría **≈ 1 día**. Usar `T(N) = 1331·(N/384)^4 s` a 2 núcleos, o al menos `P_TIME = 4.0`. |
| 32 | Memoria (`B_CELL` 82 y 206 B/celda) | OK (conservador) | A N grande la placa de la EDT (slab 64 + 2·pad ≤ 132) pesa relativamente menos que a N = 192, donde se midió. Picos estimados a N = 1024 para U(1): ~47 GB en las condiciones iniciales y ~75 GB en las medidas con espectro, por debajo de los 88 GB del plan. |
| 33 | `distance_to_strings` por placas, "exacta para d ≤ pad" | **OK, verificado** | `test_slab_edt.py`: N ∈ {40, 48, 64}, pad desde 3 hasta más de N, slabs de 3 a 64. \|D − T\| ≤ 1e-6 (redondeo de float32) donde T ≤ pad; D > pad donde T > pad; las máscaras `dist > rc` coinciden exactamente para rc ≤ pad − 2. |
| 34 | Cambio de RCAP = min(N/4, 64) | OK | Para N ≤ 256 es igual que U2. Para N > 256 solo invalida máscaras con r_c > 64 celdas, es decir ℓ < ℓ_fin − 2 ln(64/r) (≤ ℓ_fin − 5.5 para r = 4). No toca ni las puertas ni la ventana de R_Γ. |
| 35 | Desbordamiento de enteros (kern.c y numpy) a N = 1024 | OK | N³ = 1.07e9 cabe en int32, pero 10·N³ (RP³) no, y el código usa `long` en todo índice lineal: `V`, `i`, `c*V`, `n = C*V`, `((long)z*N+y)*N+x`, `cnt`. `plaquettes` devuelve `c_long`. En numpy, `digitize` va en int64 y luego en int32 (< N/2); el resto son ints de Python. No hay desbordamientos. |
| 36 | Sin checkpoint en las redes de muchas horas | RIESGO | Un reinicio de Windows Update hace perder la red en curso. Desactivar los reinicios automáticos o fijar "horas activas" durante U4b. |
| 37 | La cadena no lee U4B | **RIESGO** (expectativa) | Nada en `chain/` lee `U4B_RESULT.json`: `make chain` en J093 **no** incorpora el R_Γ nuevo, y "S6 y U4a se recalculan con ese R_Γ" no está automatizado. Es un paso manual pendiente: decirlo en la guía y en J093. |
| 38 | environment.yml: nombres en conda-forge | OK | python, numpy, scipy, matplotlib, sympy, cython, setuptools, wheel, pip, make, git, c-compiler, cxx-compiler y psutil existen. `psutil` no se usa en ningún sitio. |
| 39 | Compiladores de conda frente a gcc del sistema al compilar classy | **RIESGO alto** | Con `c-compiler` activo, setuptools compila la extensión con `x86_64-conda-linux-gnu-cc` y su sysroot. El `setup.py` de CLASS añade `-lmvec` y la ruta de libgomp del `gcc` del sistema *(memoria)*, de donde viene el error típico "cannot find -lmvec". Arreglo en la guía (paso 5): `CC=gcc LDSHARED="gcc -pthread -shared" pip install --no-build-isolation -v .`, o quitar `c-compiler`/`cxx-compiler` de environment.yml (build-essential ya da gcc). |
| 40 | `pip install .` con aislamiento de compilación | RIESGO | Si falta `pyproject.toml`, el aislamiento no ve numpy ni Cython; si existe, compila contra otro numpy. Usar `--no-build-isolation` (el entorno ya trae numpy, cython, setuptools y wheel). Si falla con numpy 2 o Cython 3 (`np.float_t`…), añadir a "Problemas frecuentes": `conda install "numpy<2" "cython<3"` y recompilar. |
| 41 | El `make` de CLASS instala classy con `--user` en algunas versiones | RIESGO bajo | Deja una copia en `~/.local` que tapa la del entorno. Después de instalar, comprobar `python -c "import classy; print(classy.__file__)"`, que debe estar en `~/miniforge3/envs/whalefall/...`. Si sale `~/.local`, borrarla. |
| 42 | Jobs: `set -u` + `conda activate` | **RIESGO** | Los scripts de activación de los compiladores de conda-forge fallan con "unbound variable" bajo `set -u` en varias versiones, y entonces el job sale con código 2. Arreglo en J091/J092/J093: `set +u; source "$CONDA_SH" && conda activate whalefall \|\| {…}; set -u`. |
| 43 | Jobs: `cd claude_dual/unification` (en /mnt/f, ruta con espacios) | RIESGO bajo | No sigue la recomendación de la guía (~/whalefall/unification). La ruta con espacios llega a `ncdm_psd_filenames`. El lector de listas de cadenas de CLASS solo corta por comas *(memoria)*, así que debería funcionar, y si falla, CLASS falla pronto, en C0. Recomendado: `cd "${UNIF_DIR:-$HOME/whalefall/unification}"`. Además, J091 no mira el código de salida de U4a: usar `python -B lya/u4a_class.py > … 2>&1 \|\| { echo "U4a falló"; exit 5; }`. |
| 44 | Guía, paso 1 | RIESGO bajo | `rsync` y `tmux` se usan pero no se instalan: `sudo apt install -y build-essential git curl rsync tmux`. |
| 45 | Guía: límite de memoria de WSL2 | **RIESGO** (decide el N) | WSL2 limita por defecto la VM al 50 % de la RAM de Windows, y el plan usa el `MemAvailable` de dentro de WSL. Añadir un paso: `%UserProfile%\.wslconfig` con `[wsl2]` / `memory=<≈90 % de la RAM>` / `swap=8GB`, seguido de `wsl --shutdown`. Sin eso, una máquina de 64 GB no pasa de N ≈ 640 en U(1). |
| 46 | Resultados del mock en `lya/results/u4a/` | RIESGO bajo | El rsync del paso 4 copia el `U4A_RESULT.json`, `psd/` y `cache/` del mock (`classy_version: "MOCK"`). La clave no reutiliza esa caché, pero el JSON se puede confundir. `rm -rf lya/results/u4a` antes de `make u4a-quick`. |
| 47 | U4a y U4b a la vez | RIESGO | Si la cola ejecuta J091 y J092 en paralelo, los N workers de CLASS más una red dimensionada al 75 % de la RAM pueden llevar a OOM. Ejecutarlos en serie: J091 → J092 → J093. |

---

## Detalle de los BUG (para aplicar tal cual)

1. **env/check_env.py, línea 58** (diccionario `base`):
   ```python
   base = {"h": 0.674, "omega_b": 0.0224, "omega_cdm": 0.0, "A_s": 2.1e-9, "n_s": 0.9649, "N_ur": 3.044, "output": "mPk", "z_pk": 3.0,
   ```
2. **lya/u4a_class.py, línea 168**:
   ```python
   strict = min(v for k, v in ref_dA.items() if not k.startswith("table")); wide = max(v for k, v in ref_dA.items() if not k.startswith("table"))
   ```
   Y añadir `deltaA_table_max=max(v for k, v in ref_dA.items() if k.startswith("table"))` al `out.update(...)` de la línea 169.
3. **Hilos.** En `lya/u4a_class.py:127`:
   ```python
   os.environ["OMP_NUM_THREADS"] = "1"
   ```
   En `Makefile`, recetas de `u4b-phase1` y `u4b-phase2`:
   ```make
   cd network && OMP_NUM_THREADS=$${U4_THREADS:-$$(nproc)} $(PY) -B u4b_produce.py 1 $(if $(N),--N $(N),)
   ```
   (y lo mismo con `2`).
4. **network/u4b_plan.py, línea 32**:
   ```python
   pick = lambda model: next((N for N in CANDIDATES if mem_GB(model, N) <= budget and (N_force is None or N <= N_force)), None)
   ```

## Lo más urgente fuera de los BUG
- El RIESGO 15 (que una fila con error no tumbe la corrida).
- El 16 (caché atómica).
- El 42 (`set +u` en los jobs).
- El 39 y el 40 (compilar classy con `CC=gcc … --no-build-isolation`).
- El 45 (`.wslconfig`).
- El 31 (horas reales de U4b, ∝ N⁴).
- Hacer `make u4a-quick` antes que nada, para medir el tiempo por ejecución de CLASS sin fluid approximation (ítem 8).

## Fuentes consultadas
- CLASS `explanatory.ini` (copia en GitHub): https://github.com/standardgalactic/class_public/blob/master/explanatory.ini
- Lesgourgues y Tram, *CLASS IV: efficient implementation of non-cold relics*, arXiv:1104.2935: https://arxiv.org/pdf/1104.2935
- Repositorio oficial (página de `input.c`, vista truncada): https://github.com/lesgourg/class_public/blob/master/source/input.c
