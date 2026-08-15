# poisselect — Poisson-Selektionsmodell für Zähldaten

## Überblick

`poisselect` schätzt ein Poisson-Selektionsmodell für Zähldaten per Maximum
Likelihood. Eine Zählvariable `y` wird nur für die selektierte Teilstichprobe
(`s = 1`) beobachtet; hängt die Selektionswahrscheinlichkeit über
unbeobachtete Größen mit `y` zusammen, ist ein klassisches Poisson-GLM auf
den beobachteten Daten verzerrt. Das Paket modelliert den Auswahlprozess mit
und korrigiert die Verzerrung durch simultane ML-Schätzung beider
Gleichungen; das Likelihood-Integral wird per Gauß-Hermite-Quadratur
(Golub-Welsch, selbst implementiert) approximiert und vollständig auf der
Log-Skala (Log-Sum-Exp-Trick) ausgewertet.

Das Modell in Kurzform:

```
Outcome:   ln mu_i = x_i'beta + eps_i,   y_i ~ Poisson(mu_i)
Selektion: s*_i    = z_i'gamma + u_i,    s_i = 1{s*_i > 0}
Fehler:    (eps_i, u_i) ~ N2(0, [sigma^2, rho*sigma; rho*sigma, 1])
```

`rho = Corr(eps, u)` misst die Stärke der Selektionsverzerrung; `rho = 0`
bedeutet ignorierbare Selektion.

## Installation

Getestet unter R 4.4.0 (macOS). Abhängigkeiten: `checkmate`, `stats`,
`graphics`, `utils` (Suggests: `testthat >= 3.0.0`); `dependencies = TRUE`
installiert `checkmate` automatisch von CRAN nach.

```r
remotes::install_local("poisselect.tar.gz",
                       dependencies = TRUE,
                       type = "source")
```

## Paketdateien

| Datei | Inhalt |
|---|---|
| `DESCRIPTION` | Paket-Metadaten, Version 0.1.0, Imports/Suggests |
| `NAMESPACE` | Von roxygen2 generiert — nie von Hand editiert |
| `LICENSE` | MIT-Lizenz (Jahr, Rechteinhaber) |
| `README.md` | Diese Datei |
| `R/poisselect-package.R` | Paket-Dokumentation (`_PACKAGE`-Sentinel) |
| `R/gauss_hermite.R` | `compute_gh_rule()`: Gauß-Hermite-Knoten/-Gewichte via Golub-Welsch-Eigenwertmethode |
| `R/check_inputs.R` | Alle Input-Checks: Argument-Ebene (checkmate) und datenabhängige Checks mit informativen Fehlermeldungen; `build_model_data()` baut y, s, X, Z |
| `R/loglik.R` | `compute_poisselect_loglik()`: approximierte Log-Likelihood, vektorisiert über n x K, Log-Sum-Exp; gemeinsamer Helfer `compute_node_matrices()` (DRY mit plot) |
| `R/starting_values.R` | Startwerte: Poisson-GLM (beta), Probit-GLM (gamma), sigma = 1, rho = 0, mit Fallbacks |
| `R/poisselect.R` | Hauptfunktion `poisselect()`: Checks, Quadratur, Optimierung (BFGS auf log(sigma)/atanh(rho)-Skala), Hesse/Kovarianz, S3-Objektbau, Warnungen |
| `R/methods_print.R` | `print.poisselect` |
| `R/methods_summary.R` | `summary.poisselect` + `print.summary.poisselect` (Koeffiziententabellen mit z- und p-Werten, rho, Log-Likelihood, AIC) |
| `R/methods_plot.R` | `plot.poisselect`: Plot 1 beobachtete vs. modell-implizierte Count-Verteilung, Plot 2 Log-Likelihood entlang rho |
| `R/methods_predict.R` | `predict.poisselect` mit `type = "response"/"link"/"pselect"` |
| `R/methods_extractors.R` | `coef`, `vcov`, `logLik` |
| `R/simulate.R` | `simulate_poisselect()`: Datengenerator exakt nach Modell |
| `man/*.Rd` | Von roxygen2 generierte Hilfeseiten (alle exportierten Funktionen und Methoden) |
| `tests/testthat.R` | Standard-testthat-Einstieg |
| `tests/testthat/helper.R` | Gemeinsame Test-Fixtures (seeded Datensatz + Fit) |
| `tests/testthat/test-gauss_hermite.R` | Bekannte Regeln K = 1/2, exakte GH-Identitäten, Fehlinputs |
| `tests/testthat/test-input_checks.R` | Ein Test pro Fehlinput aus der Spezifikation (Formeln, k, start, control, s, y, NA, Kollinearität, Fallzahlen) + Warnungen |
| `tests/testthat/test-loglik.R` | Endlichkeit, Log-Sum-Exp-Robustheit bei großen Counts, analytischer Grenzfall rho = 0, sigma -> 0 |
| `tests/testthat/test-fit_recovery.R` | Parameter-Recovery auf simulierten Daten, Vergleich mit GLM bei rho = 0, Stabilität in K, Startwert-Robustheit |
| `tests/testthat/test-methods.R` | Alle S3-Methoden inkl. predict-Handrechnung und Plot-Aufrufe |
| `tests/testthat/test-simulate.R` | Struktur des Datensatzes, NA-Muster, Fehlinputs |
| `tests/testthat/test-reproducibility.R` | Gleicher Seed -> identische Daten; gleiche Daten -> gleiche Schätzung |

## Verwendung

```r
library(poisselect)

set.seed(1)
d <- simulate_poisselect(800, beta = c(1, 0.5), gamma = c(0.5, 1, -0.8),
                         sigma = 0.5, rho = 0.6)
fit <- poisselect(y ~ x1, s ~ x1 + w1, data = d)   # K = 20 (Default)
summary(fit)
plot(fit)
head(predict(fit))                  # E[Y | x] = exp(x'beta + sigma^2/2)
predict(fit, type = "pselect")[1:5] # Selektionswahrscheinlichkeiten
```

## Tests ausführen

```r
devtools::test()                        # im Paketordner
# oder nach Installation:
testthat::test_package("poisselect")
```

## Zu beachten

* **Quadraturknoten `k`:** Argument von `poisselect()`, Default `K = 20`
  (zulässig 2–500). Mehr Knoten = genauere Approximation; die Schätzung ist
  ab ca. K = 10 sehr stabil (getestet K = 10 vs. 30).
* **`y` darf/soll `NA` sein, wo `s = 0`:** unbeobachtete Outcomes werden
  ignoriert. `NA` bei `s = 1` ist dagegen ein Fehler. Der
  Selektionsindikator selbst darf nie `NA` sein.
* **Reparametrisierung:** Optimiert wird unrestringiert über
  `log(sigma)` und `atanh(rho)`; die Rücktransformationen `exp`/`tanh`
  erzwingen `sigma > 0` und `-1 < rho < 1` automatisch. Die
  Standardfehler von `beta`/`gamma` sind davon unberührt: die Jacobi-Matrix
  der Rücktransformation ist auf dem Koeffizientenblock die Einheitsmatrix,
  der (beta, gamma)-Block der inversen Hesse-Matrix ist also invariant.
* **Warnungen und ihre Bedeutung:**
  * *Nicht-Konvergenz* (optim-Code != 0): mehr Iterationen über
    `control = list(maxit = ...)` versuchen.
  * *Randlösung* (`|rho| > 0.99` oder `sigma < 1e-4`): Schätzer am Rand
    des Parameterraums, Inferenz unzuverlässig.
  * *Fehlende Exclusion Restriction*: alle Selektionskovariablen kommen
    auch in der Outcome-Gleichung vor — das Modell ist dann nur über die
    funktionale Form identifiziert.
  * *Hesse-Matrix nicht invertierbar/positiv definit*: Standardfehler und
    `vcov` werden `NA`, die Punktschätzer bleiben nutzbar.
* **Bekannte Grenzen:** bei `rho` nahe +-1 oder sehr kleinen Stichproben
  kann die Likelihood flach sein (Randlösungs-Warnungen ernst nehmen);
  Startwerte lassen sich über `start = c(beta, gamma, sigma, rho)`
  vorgeben.
* **lintr/Checks:** `lintr::lint_package()` meldet 0 Lints (keine
  `.lintr`-Sonderkonfiguration nötig). Hinweis: der
  `object_usage_linter` von lintr kann paketinterne Funktionen erst
  auflösen, wenn das Paket installiert ist — vorher zeigt er
  Scheinbefunde. `R CMD check` läuft mit 0 Errors, 0 Warnings, 0 Notes;
  ohne Internetzugang erscheint eine umgebungsbedingte NOTE
  „unable to verify current time" (R kann die Systemzeit nicht gegen
  einen Zeitserver prüfen — kein Paketproblem).
