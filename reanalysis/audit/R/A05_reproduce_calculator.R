# Audit Q7: independent reproduction of the nine-item model and of the calculator (R vs stored spec vs JavaScript vs rendered web page).
source("audit/R/00_audit_setup.R")
md_open("A05_reproduce_calculator.R", "A05 Reproduction of the nine-item model and the calculator")
d <- load_dev()

# ---- 1. refit and re-run the bootstrap exactly as R/16 (same seed, same order of random draws) -----------------
set.seed(SEED)
f <- coxph(fml_of(PRED9), data = d, x = TRUE)
cidx <- function(lp, dd) harrell(lp, dd)$concordance
boot <- rbindlist(lapply(1:200, function(b) { i <- sample(nrow(d), replace = TRUE); fb <- coxph(fml_of(PRED9), data = d[i])
  lpo <- predict(fb, newdata = d, type = "lp")
  data.table(opt = cidx(predict(fb, newdata = d[i], type = "lp"), d[i]) - cidx(lpo, d), slope = coef(coxph(Surv(t, y) ~ lpo, data = d))) }))
shrink <- mean(boot$slope); S0 <- s0_at(f, 5)
rep_tab <- data.table(quantity = c("S0(5)", "shrinkage", "mean age (centre)", "mean Hb deficit (centre)", paste("beta", PRED9), "apparent C", "optimism-corrected C"),
  reproduced = c(S0, shrink, f$means["age"], f$means["hb_def"], coef(f)[PRED9], cidx(predict(f, type = "lp"), d), cidx(predict(f, type = "lp"), d) - mean(boot$opt)),
  stored = c(SPEC9$s0, SPEC9$shrink, SPEC9$means$age, SPEC9$means$hb_def, unlist(SPEC9$coef[PRED9]), as.numeric(SPEC9$validation[["Harrell C, apparent"]]), as.numeric(SPEC9$validation[["Harrell C, optimism-corrected"]])))
rep_tab[, abs_difference := signif(abs(reproduced - stored), 3)]
rep_tab[, `:=`(reproduced = signif(reproduced, 8), stored = signif(stored, 8))]
aud_tab(rep_tab, "A05a_model_reproduction", "Refit of the nine-item model and re-run of the 200-resample bootstrap (seed 20261022) vs output/calculator_model_spec.json. Stored C values are rounded to 3 decimals.")
md_note("Coefficients in the spec and in R8 are the ORIGINAL (unshrunk) Cox coefficients; the calculator multiplies the centred linear predictor by the shrinkage factor (", f2(shrink, 4),
        "). Hazard ratios printed in the methods summary are therefore unshrunk. The baseline survival is basehaz(centered = TRUE) at the coxph means, which are 0 for binary predictors — the same convention the calculator uses.",
        " Bootstrap optimism SD = ", f2(sd(boot$opt), 4), "; Monte-Carlo SE of the shrinkage factor = ", f2(sd(boot$slope) / sqrt(200), 4), ".")

# ---- 2. test profiles -----------------------------------------------------------------------------------------
P <- fread(text = "id,label,age,hb,migraine_any,estrogen,dm,htn,haemoglobinopathy,cad,ckd_pre
1,low risk,30,13.5,0,0,0,0,0,0,0
2,typical (cohort means),41.1,12.5,0,0,0,0,0,0,0
3,calculator example,48,9.4,1,0,0,1,0,0,0
4,young with migraine and estrogen,25,11,1,1,0,0,0,0,0
5,high risk,55,7.5,0,0,1,1,0,1,1
6,all items positive,18,3,1,1,1,1,1,1,1
7,age boundary below 45,44.99,12,0,0,0,0,0,0,0
8,age boundary at 45,45,12,0,0,0,0,0,0,0
9,Hb boundary 11.99,40,11.99,0,0,0,0,0,0,0
10,Hb boundary 8.0,40,8,0,0,0,0,0,0,0
11,Hb boundary 7.99,40,7.99,0,0,0,0,0,0,0
12,hemoglobinopathy,35,10,0,0,0,0,1,0,0
13,oldest no risk factors,60,20,0,0,0,0,0,0,0
14,CAD only,50,13,0,0,0,0,0,1,0")
P[, hb_def := pmax(13 - hb, 0)]
# independent R computation from the refit (not from the spec)
P[, r_independent := 100 * (1 - S0^exp(shrink * drop(sweep(as.matrix(.SD[, PRED9, with = FALSE]), 2, f$means[PRED9]) %*% coef(f)[PRED9])))]
P[, r_spec := 100 * risk_spec(P)]
P[, pts_r := as.integer(age < 45) + fifelse(hb >= 12, 0L, fifelse(hb >= 8, 1L, 2L)) + 2L * migraine_any + 2L * estrogen + dm + 2L * htn + 3L * haemoglobinopathy + 4L * cad + 2L * ckd_pre]
P[, pts_risk_r := 100 * (1 - SPEC14$score_baseline_survival_5y^exp(SPEC14$score_coef * (pts_r - 1)))]
fwrite(P, file.path(AUD, "tables", "A05_profiles_input.csv"))

# ---- 3. JavaScript (MODEL block evaluated in node) and the rendered page (Playwright; local file and the public URL) ----------
# the public page is fetched with curl (the headless browser does not trust the session's egress proxy CA, and TLS checks are not disabled);
# it is compared byte-for-byte with the repository copy and the downloaded copy is driven in the browser
dir.create(file.path(AUD, "public_snapshot"), FALSE)
snap <- file.path(AUD, "public_snapshot", "index.html")
system2("curl", c("-sS", "https://selim-tarek.github.io/Med-Hack-score/", "-o", snap))
same <- identical(tools::md5sum(snap)[[1]], tools::md5sum("calculator/standalone/index.html")[[1]])
md_note("Public page https://selim-tarek.github.io/Med-Hack-score/ fetched ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; byte-identical to calculator/standalone/index.html: ", same, ".")
js <- normalizePath("audit/js/calc_test.js")
st <- system2("node", c(js, normalizePath(file.path(AUD, "tables", "A05_profiles_input.csv")), normalizePath("calculator/standalone/index.html"),
                        normalizePath(file.path(AUD, "tables")), paste0("file://", normalizePath(snap))),
              env = paste0("NODE_PATH=", system("npm root -g", intern = TRUE)), stdout = TRUE, stderr = TRUE)
aud_log(paste(st, collapse = " / "))
J <- fread(file.path(AUD, "tables", "A05_js_web_output.csv"))
P <- merge(P, J, by = "id", all.x = TRUE)
out <- P[, .(id, label, age, hb, points_r = pts_r, points_js = js_points, points_web_local = web_local_points, points_web_public = web_public_points,
             risk_model_r_independent = round(r_independent, 4), risk_model_spec = round(r_spec, 4), risk_model_js = round(js_model, 4),
             risk_model_web_public = web_public_model, risk_points_r = round(pts_risk_r, 4), risk_points_js = round(js_points_risk, 4), risk_points_web_public = web_public_points_risk)]
aud_tab(out, "A05b_profile_predictions", "Predicted 5-year risk (%) and points for test profiles. Web values are as displayed (rounded to 1 decimal).")
disc <- P[, .(comparison = c("model: independent R vs spec", "model: spec vs JavaScript", "model: JavaScript vs public page (display rounding 0.05)",
                             "points: R vs JavaScript", "points: JavaScript vs local page", "points: JavaScript vs public page", "points-table risk: R vs JavaScript"),
              max_abs_difference = c(max(abs(r_independent - r_spec)), max(abs(r_spec - js_model)), max(abs(round(js_model, 1) - web_public_model), na.rm = TRUE),
                                     max(abs(pts_r - js_points)), max(abs(js_points - web_local_points)), max(abs(js_points - web_public_points), na.rm = TRUE), max(abs(pts_risk_r - js_points_risk))))]
disc[, max_abs_difference := signif(max_abs_difference, 3)]
aud_tab(disc, "A05c_discrepancies", "Maximum absolute discrepancy across the 14 profiles (risk in percentage points)")
md_note("Unit conversion on the public page: example profile entered as 94 g/L gave ", P[id == 3, web_public_gl], "; as 5.835 mmol/L gave ", P[id == 3, web_public_mmol],
        " (g/dL entry: ", P[id == 3, web_public_points], " points, ", P[id == 3, web_public_model], "% model risk).")
md_note("Two different 5-year risks are shown for every patient: the headline 'points' risk comes from a Cox model on the 14-predictor points total (risk_model_spec.json), the 'exact model' risk from the nine-item model. For the calculator example they are ",
        f2(P[id == 3, pts_risk_r], 2), "% and ", f2(P[id == 3, r_spec], 2), "%.")
res_put("A05", list(reproduction = rep_tab, discrepancies = disc))
