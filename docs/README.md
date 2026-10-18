# Study codes and where they appear in the paper

Each study has a dated specification written before its networks were produced (`preregistros/ESPEC_U*.md`) and a
report with the verdict (`informes/INFORME_U*.md`). Both are the original records, in Spanish; amendments are dated
and appended at the end. Results are in `network/results/`.

| Code | What it measured | Results | Paper |
|---|---|---|---|
| U2 | First RP³ and U(1) networks (256³); kernel validations V1–V4 | `network/results/u2/` | §2 Validations; §3.1 (256³ memory) |
| U3 (U3b) | Networks with core preparation (the Lyα part of U3 belongs to the larger project) | `network/results/u3/` | — (method history) |
| U4 (U4b) | Large networks: U(1) attractor gates at 768³, RP³ vs U(1) at 576³ (the CLASS part of U4 belongs to the larger project) | `network/results/u4/` | §3.1 |
| U6 | Six RP³ and four U(1) seeds at 576³: R_ξ, undrained R_Γ, q(ℓ) | `network/results/u4/U6_RESULT.json` | §3.2, §3.3 |
| U9 | Dense starts (ξ_i = 2) | `network/results/u4/U9_RESULT.json` | §3.2, §3.3 |
| U10a | Drain of the initial massive bath and direct radial-to-Goldstone exchange; 512³ protocol test | `network/results/u4/U10A_*.json` | §3.3, §3.4 |
| U10b | Late dense starts approaching the attractor from above | `network/results/u4/U10B_RESULT.json` | §3.2, Fig. 1 |
| U11 | The U6 networks with the bath drained: clean R_Γ, R_ξ and q | `network/results/u4/U11_RESULT.json` | §3.3, §3.4, Fig. 2 |
| T7 | Fading of the bath contribution within the window | `network/results/u4/T7_DRAIN_FADE.json` | §3.3 |
| U16 | Resolution (192–384³ at fixed LH) and volume (×4.6 Hubble volumes) | `network/results/u16/` | §2 Validations |
| U18 | Snapshot sampling of the local and integral emission estimators | `network/results/u18/` | §2 Validations, R_Γ sampling error |
| Audit U4–U11 | Independent re-implementation of the estimators | `audit/u4_u11/` | §5 |

`informes/NOTA_DRENAJE.md` explains why draining the initial bath removes a lattice artefact and not part of the
cosmological production (§3.3, last paragraph).
