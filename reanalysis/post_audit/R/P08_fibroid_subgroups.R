# Part 9: is the model fibroid-specific? Performance of the reduced baseline model (C4) and the current nine (C2) by diagnosis, with CIs.
source("post_audit/R/P00_setup.R")
md_open("P08_fibroid_subgroups.R", "P08 Performance by uterine diagnosis")
b <- load_base(); b[, grp := fifelse(fibroids == 1, "any fibroids", "no fibroids (endometriosis and/or adenomyosis only)")]
dv <- b[year(anchor) <= 2022]; vv <- b[year(anchor) >= 2023]
out <- rbindlist(lapply(c("C4_reduced_baseline", "C2_current_nine_as_coded"), function(n) {
  fa <- coxph(fml_of(CAND[[n]]), data = b, x = TRUE, model = TRUE); fv <- coxph(fml_of(CAND[[n]]), data = dv, x = TRUE, model = TRUE)
  rbindlist(lapply(unique(b$grp), function(g) { x <- b[grp == g]; z <- vv[grp == g]; ca <- cstat(lpof(fa, x), x, H); cv <- cstat(lpof(fv, z), z, H); oe <- oe_counts(fv, z, H)
    data.table(model = n, group = g, women = nrow(x), events_3y = x[y == 1 & t <= H, .N], apparent_c3 = sprintf("%.2f (%.2f–%.2f)", ca[1], ca[1] - 1.96 * ca[2], ca[1] + 1.96 * ca[2]),
               temporal_women = nrow(z), temporal_events_3y = z[y == 1 & t <= H, .N], temporal_c3 = sprintf("%.2f (%.2f–%.2f)", cv[1], cv[1] - 1.96 * cv[2], cv[1] + 1.96 * cv[2]),
               temporal_oe3 = fmt_ci(oe[c("oe", "lo", "hi")], 2)) })) }))
aud_tab(out, "P08a_by_fibroid_status", "Apparent (all years) and temporal (develop 2018–2022, evaluate 2023+) performance by fibroid status, baseline cohort, primary outcome. Overlapping CIs are not evidence of equivalence.")
rt <- b[, .(women = .N, events = sum(y), rate_per_1000py = pois_ci(sum(y), sum(t))), by = grp]
aud_tab(rt, "P08b_rates", "Crude rate of acute-coded stroke by fibroid status (no comparison group without uterine disease exists)")
res_put("P08", list(by_fibroid = out))
