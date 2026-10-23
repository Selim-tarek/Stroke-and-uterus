# Part 6: verify P08 and test the subgroup claim formally (difference in C between subgroups with a bootstrap CI; heterogeneity of the
# linear predictor's effect by fibroid status). The 'migraine explains it' interpretation gets one exploratory check, labelled as such.
source("final_verification/R/V00_setup.R")
md_open("V04_fibroid.R", "V04 Fibroid subgroup verification")
b <- load_base(); b[, fib := as.integer(fibroids == 1)]
dv <- b[year(anchor) <= 2022]; vv <- b[year(anchor) >= 2023]
c_direct <- function(lp, dd, h) { yy <- dd$y * (dd$t <= h); tt <- pmin(dd$t, h); k <- concordance(Surv(tt, yy) ~ lp, reverse = TRUE); c(k$concordance, sqrt(k$var)) }
M <- list(C4 = CAND$C4_reduced_baseline, C2 = CAND$C2_current_nine_as_coded, C4_plus_migraine_as_coded = c(CAND$C4_reduced_baseline, "migraine_any"))
fits <- lapply(M, fitc2, dd = dv); for (n in names(M)) vv[, (paste0("lp_", n)) := lpof(fits[[n]], vv)]
tab <- rbindlist(lapply(names(M), function(n) rbindlist(lapply(0:1, function(g) { x <- vv[fib == g]
  data.table(model = n, subgroup = c("no fibroids", "fibroids")[g + 1], women = nrow(x), events_3y = x[y == 1 & t <= H, .N], c3 = ci_c(c_direct(x[[paste0("lp_", n)]], x, H), 2)) }))))
aud_tab(tab, "V04a_subgroup_c", "Temporal C at 3 years by fibroid status (develop 2018–2022 on all women, evaluate 2023+). Exact definition: fibroids = Codebook flag fibroids == 1 (any coded leiomyoma in the 2018+ extract); 'no fibroids' = endometriosis and/or adenomyosis only.")
set.seed(SEED)
bt <- rbindlist(lapply(1:1000, function(i) { r <- list()
  for (g in 0:1) { x <- vv[fib == g]; x <- x[sample(.N, replace = TRUE)]; for (n in names(M)) r[[paste0(n, "_", g)]] <- c_direct(x[[paste0("lp_", n)]], x, H)[1] }
  as.data.table(r) }))
q <- function(z) sprintf("%.3f (%.3f to %.3f)", mean(z), quantile(z, .025), quantile(z, .975))
aud_tab(data.table(contrast = c("C4: fibroids − no fibroids", "C2: fibroids − no fibroids", "no fibroids: C2 − C4", "fibroids: C2 − C4", "no fibroids: (C4 + migraine as coded) − C4"),
                   delta_c_95ci = c(q(bt$C4_1 - bt$C4_0), q(bt$C2_1 - bt$C2_0), q(bt$C2_0 - bt$C4_0), q(bt$C2_1 - bt$C4_1), q(bt$C4_plus_migraine_as_coded_0 - bt$C4_0))),
        "V04b_subgroup_differences", "1,000 bootstrap resamples drawn within each subgroup of the 2023+ set (stratified; within-subgroup contrasts are paired). Exploratory, not pre-specified.")
x3 <- trunc_h(vv, H)
g1 <- coxph(Surv(t, y) ~ lp_C4 * fib, data = x3); g2 <- coxph(Surv(t, y) ~ lp_C2 * fib, data = x3)
it <- function(g, k) sprintf("ratio of slopes %.2f (%.2f–%.2f), P %s", exp(coef(g)[k]), exp(confint(g)[k, 1]), exp(confint(g)[k, 2]), fp(summary(g)$coefficients[k, 5]))
aud_tab(data.table(model = c("C4", "C2"), interaction_lp_by_fibroids = c(it(g1, "lp_C4:fib"), it(g2, "lp_C2:fib"))),
        "V04c_interaction", "Heterogeneity test in the 2023+ set: Cox model of the outcome on the development linear predictor × fibroid status (3 years). A ratio ≠ 1 means the predictors' combined effect differs by subgroup.")
res_put("V04", list(c = tab))
