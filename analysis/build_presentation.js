// Anemia and stroke in benign uterine disease - presentation deck (numbers from outputs/_state.pkl, same run)
// Build: npm install pptxgenjs, then node analysis/build_presentation.js (run the pipeline first for figures)
const pptxgen = require("pptxgenjs");
const FIG = require("path").join(__dirname, "..", "outputs", "figures") + "/";
const OUT = require("path").join(__dirname, "..", "outputs", "presentation", "anemia_stroke_presentation.pptx");

const C = {
  crimson: "8E1B2C", dark: "2B1115", ink: "1F1F1F", muted: "5E5A5B", light: "FFFFFF",
  tint: "F7ECEE", tint2: "EFE3E5", navy: "2F3C7E", good: "2E7D5B", line: "D9CFD1",
};
const HF = "Cambria", BF = "Calibri";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5
pres.title = "Anemia and stroke in women with benign uterine disease";

const W = 13.333, H = 7.5, M = 0.6;

function title(slide, text, sub) {
  slide.addText(text, { x: M, y: 0.4, w: W - 2 * M, h: 0.8, fontFace: HF, fontSize: 34, bold: true,
    color: C.dark, margin: 0, isTextBox: true });
  if (sub) slide.addText(sub, { x: M, y: 1.18, w: W - 2 * M, h: 0.45, fontFace: BF, fontSize: 16,
    color: C.muted, margin: 0, isTextBox: true });
}
function shadow() { return { type: "outer", color: "000000", blur: 6, offset: 2, angle: 90, opacity: 0.12 }; }
function card(slide, x, y, w, h, fill) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: fill || C.tint },
    line: { color: fill || C.tint }, shadow: shadow() });
}
function badge(slide, x, y, n, fill) {
  slide.addShape(pres.shapes.OVAL, { x, y, w: 0.55, h: 0.55, fill: { color: fill || C.crimson }, line: { color: fill || C.crimson } });
  slide.addText(String(n), { x, y, w: 0.55, h: 0.55, fontFace: HF, fontSize: 18, bold: true, color: C.light,
    align: "center", valign: "middle", margin: 0, isTextBox: true });
}
function footer(slide, text) {
  slide.addText(text, { x: M, y: H - 0.45, w: W - 2 * M, h: 0.3, fontFace: BF, fontSize: 10, color: C.muted,
    margin: 0, isTextBox: true });
}

// 1. Title -----------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.dark };
  s.addShape(pres.shapes.OVAL, { x: 9.3, y: -1.2, w: 5.6, h: 5.6, fill: { color: C.crimson, transparency: 35 }, line: { color: C.crimson, transparency: 35 } });
  s.addShape(pres.shapes.OVAL, { x: 10.9, y: 3.9, w: 3.2, h: 3.2, fill: { color: C.crimson, transparency: 60 }, line: { color: C.crimson, transparency: 60 } });
  s.addText("Beyond heavy periods", { x: M, y: 1.7, w: 9, h: 0.6, fontFace: BF, fontSize: 20, color: "E8B9C0",
    italic: true, margin: 0, isTextBox: true });
  s.addText("Anemia and stroke risk in women with fibroids, adenomyosis and endometriosis",
    { x: M, y: 2.35, w: 9.4, h: 2.3, fontFace: HF, fontSize: 34, bold: true, color: C.light, margin: 0,
      valign: "top", isTextBox: true });
  s.addText("39,807 women · Mayo Clinic electronic health records · 2018–2026",
    { x: M, y: 4.85, w: 9, h: 0.4, fontFace: BF, fontSize: 16, color: "D9C4C7", margin: 0, isTextBox: true });
  s.addText("[Presenter name] · [Department] · [Date]", { x: M, y: 6.4, w: 9, h: 0.4, fontFace: BF, fontSize: 14,
    color: "B8A3A6", margin: 0, isTextBox: true });
  s.addNotes("Title. This study asks whether anemia, which is very common in women with benign uterine disease, is linked to stroke. It uses electronic health records from about 40,000 women.");
}

// 2. Why ask -------------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Why ask this question?");
  const items = [
    ["Anemia is common here", "Heavy bleeding from fibroids and adenomyosis makes anemia one of the most frequent problems in these women: 22% of our cohort were anemic."],
    ["Stroke in young women", "Stroke under 60 is uncommon but devastating, and conventional risk factors explain only part of it [ref]."],
    ["Hemoglobin is free", "Hemoglobin is measured routinely at gynecological diagnosis. If it flags vascular risk, the information is already in the chart."],
  ];
  items.forEach(([h, t], i) => {
    const x = M + i * 4.1;
    card(s, x, 1.9, 3.8, 4.3);
    badge(s, x + 0.35, 2.2, i + 1);
    s.addText(h, { x: x + 0.35, y: 3.0, w: 3.1, h: 0.9, fontFace: HF, fontSize: 21, bold: true, color: C.dark,
      margin: 0, valign: "top", isTextBox: true });
    s.addText(t, { x: x + 0.35, y: 3.9, w: 3.1, h: 2.1, fontFace: BF, fontSize: 15, color: C.ink, margin: 0,
      valign: "top", isTextBox: true });
  });
  s.addText("Question: is anemia at gynecological diagnosis associated with stroke?", { x: M, y: 6.5, w: W - 2 * M,
    h: 0.5, fontFace: BF, fontSize: 17, bold: true, color: C.crimson, margin: 0, isTextBox: true });
  s.addNotes("Three reasons. Anemia is very common in this group. Stroke in younger women is poorly explained by standard risk factors. And hemoglobin costs nothing extra because it is already measured. [ref] marks where a literature citation is needed.");
}

// 3. Study design ---------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Who we studied", "Retrospective study of electronic health records");
  const steps = [["52,827", "records reviewed"], ["39,807", "eligible women aged 18–60 with fibroids, adenomyosis or endometriosis"],
    ["32,683", "had a hemoglobin (Hb) test"], ["7,306", "were anemic (Hb <12 g/dL)"]];
  steps.forEach(([n, t], i) => {
    const y = 1.95 + i * 1.18;
    card(s, M, y, 5.6, 0.98, i === 3 ? C.crimson : C.tint);
    s.addText(n, { x: M + 0.3, y, w: 1.9, h: 0.98, fontFace: HF, fontSize: 28, bold: true,
      color: i === 3 ? C.light : C.crimson, valign: "middle", margin: 0, isTextBox: true });
    s.addText(t, { x: M + 2.2, y, w: 3.25, h: 0.98, fontFace: BF, fontSize: 14, color: i === 3 ? C.light : C.ink,
      valign: "middle", margin: 0, isTextBox: true });
  });
  const x2 = 7.0;
  s.addText("Anemia grade (WHO)", { x: x2, y: 1.95, w: 5.7, h: 0.45, fontFace: HF, fontSize: 20, bold: true,
    color: C.dark, margin: 0, isTextBox: true });
  const rows = [["Grade", "Hb (g/dL)", "Women"], ["None", "≥12", "25,377"], ["Mild", "10–11.9", "5,293"],
    ["Moderate", "8–9.9", "1,476"], ["Severe", "<8", "537"]];
  s.addTable(rows.map((r, i) => r.map(v => ({ text: v, options: { bold: i === 0, color: i === 0 ? C.light : C.ink,
    fill: { color: i === 0 ? C.crimson : (i % 2 ? "FFFFFF" : C.tint) } } }))),
    { x: x2, y: 2.5, w: 5.7, colW: [2.1, 1.8, 1.8], rowH: 0.42, fontFace: BF, fontSize: 14,
      border: { type: "solid", color: C.line, pt: 0.5 } });
  s.addText([
    { text: "Outcome: ", options: { bold: true } }, { text: "stroke or TIA (1,538 women; 1,189 ischemic)", options: { breakLine: true } },
    { text: "Adjusted for: ", options: { bold: true } }, { text: "age, race, BMI, blood pressure, diabetes, cholesterol, smoking, migraine, thrombophilia, hormones, bleeding, diagnosis group" },
  ], { x: x2, y: 4.85, w: 5.7, h: 1.6, fontFace: BF, fontSize: 14, color: C.ink, margin: 0, valign: "top",
    paraSpaceAfter: 6, isTextBox: true });
  s.addNotes("52,827 records, 39,807 eligible. 32,683 had a hemoglobin test and 7,306 of them were anemic. Grades follow WHO cut-offs. 1,538 women had a stroke or TIA. Every model adjusts for the standard vascular risk factors plus migraine, hormones and bleeding.");
}

// 4. Headline numbers ------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.tint };
  title(s, "The headline");
  const stats = [["22%", "of women were anemic"], ["1.8×", "odds of stroke with moderate anemia"], ["2.0×", "odds of stroke with severe anemia"]];
  stats.forEach(([n, t], i) => {
    const x = M + i * 4.1;
    card(s, x, 1.8, 3.8, 3.3, C.light);
    s.addText(n, { x, y: 2.05, w: 3.8, h: 1.6, fontFace: HF, fontSize: 72, bold: true,
      color: i === 0 ? C.navy : C.crimson, align: "center", valign: "middle", margin: 0, isTextBox: true });
    s.addText(t, { x: x + 0.3, y: 3.7, w: 3.2, h: 1.1, fontFace: BF, fontSize: 17, color: C.ink, align: "center",
      valign: "top", margin: 0, isTextBox: true });
  });
  s.addText("Compared with Hb ≥12 g/dL, after adjusting for vascular risk factors. Moderate OR 1.77 (95% CI 1.41–2.21); severe OR 2.04 (1.40–2.97); 31,472 women, 1,460 strokes.",
    { x: M, y: 5.5, w: W - 2 * M, h: 0.7, fontFace: BF, fontSize: 14, color: C.muted, margin: 0, isTextBox: true });
  s.addNotes("The headline: one in five women was anemic. Moderate anemia went with about 1.8 times the odds of stroke and severe anemia about 2 times, after adjusting for the usual risk factors.");
}

// 5. Dose response (native chart) ------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Lower hemoglobin, higher risk", "Adjusted odds ratio for stroke, by anemia grade");
  s.addChart(pres.charts.BAR, [{ name: "Adjusted OR", labels: ["None (≥12)", "Mild (10–11.9)", "Moderate (8–9.9)", "Severe (<8)"],
    values: [1.0, 1.19, 1.77, 2.04] }], {
    x: M, y: 1.8, w: 7.4, h: 4.9, barDir: "col", chartColors: ["C9A3A9", "B5646F", C.crimson, "5E1220"],
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.00", dataLabelFontSize: 14,
    dataLabelColor: C.ink, catAxisLabelColor: C.muted, valAxisLabelColor: C.muted, catAxisLabelFontSize: 13,
    valAxisMinVal: 0, valAxisMaxVal: 2.5, valAxisMajorUnit: 0.5, valGridLine: { color: "E5E5E5", size: 0.5 },
    catGridLine: { style: "none" }, showLegend: false, varyColors: true, barGapWidthPct: 60,
    valAxisTitle: "Odds ratio (reference = no anemia)", showValAxisTitle: true, valAxisTitleFontSize: 12,
    valAxisTitleColor: C.muted,
  });
  const x2 = 8.4;
  card(s, x2, 1.9, 4.3, 4.7);
  s.addText([
    { text: "95% confidence intervals", options: { bold: true, breakLine: true } },
    { text: "Mild 1.19 (1.03–1.38)", options: { breakLine: true } },
    { text: "Moderate 1.77 (1.41–2.21)", options: { breakLine: true } },
    { text: "Severe 2.04 (1.40–2.97)", options: { breakLine: true } },
    { text: " ", options: { breakLine: true } },
    { text: "Each step down in grade: ", options: { bold: true } }, { text: "+28% (OR 1.28, 1.18–1.38; P<0.001)", options: { breakLine: true } },
    { text: " ", options: { breakLine: true } },
    { text: "Ischemic stroke only: ", options: { bold: true } }, { text: "+33% per grade" },
  ], { x: x2 + 0.35, y: 2.15, w: 3.7, h: 4.3, fontFace: BF, fontSize: 15, color: C.ink, margin: 0, valign: "top",
    paraSpaceAfter: 4, isTextBox: true });
  s.addNotes("A clear staircase: risk rises with each anemia grade. Mild 1.19, moderate 1.77, severe 2.04. Each step adds about 28 percent, and the pattern is slightly stronger for ischemic stroke. A dose-response like this is one of the strongest signs that an association is real.");
}

// 6. Hb curve ---------------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Risk starts to rise below 13 g/dL", "Stroke risk across the whole hemoglobin range (reference 13 g/dL)");
  const iw = 8.3, ih = iw / 1.92;
  s.addImage({ path: FIG + "fig14_p2_hb_curves.png", x: M, y: 1.9, w: iw, h: ih });
  const x2 = M + iw + 0.35, w2 = W - M - x2;
  const pts = [["13 g/dL", "1.00"], ["11 g/dL", "1.19"], ["9 g/dL", "1.50"], ["8 g/dL", "1.68"]];
  s.addText("Odds of stroke vs 13 g/dL", { x: x2, y: 1.95, w: w2, h: 0.4, fontFace: HF, fontSize: 15, bold: true,
    color: C.dark, margin: 0, isTextBox: true });
  pts.forEach(([a, b], i) => {
    const y = 2.45 + i * 0.72;
    card(s, x2, y, w2, 0.6, i === 3 ? C.crimson : C.tint);
    s.addText(a, { x: x2 + 0.2, y, w: 1.6, h: 0.6, fontFace: BF, fontSize: 15, color: i === 3 ? C.light : C.ink,
      valign: "middle", margin: 0, isTextBox: true });
    s.addText(b + "×", { x: x2 + 1.8, y, w: w2 - 2.0, h: 0.6, fontFace: HF, fontSize: 20, bold: true,
      color: i === 3 ? C.light : C.crimson, align: "right", valign: "middle", margin: 0, isTextBox: true });
  });
  s.addText("Below 13 g/dL, each 1 g/dL lower Hb: about +11% odds, +10% stroke rate. Above 13: flat.",
    { x: x2, y: 5.45, w: w2, h: 1.2, fontFace: BF, fontSize: 14, color: C.ink, margin: 0, valign: "top", isTextBox: true });
  footer(s, "Restricted cubic spline, 4 knots; adjusted for covariates and conditions documented before the Hb. (A) cross-sectional odds; (B) rate after the Hb test.");
  s.addNotes("Treating hemoglobin as continuous: above 13 the risk is flat. Below 13 it rises steadily, about 10 percent for every gram lower. At 8 grams the odds are about 1.7 times those at 13. Panel B shows the same shape when we follow women forward from the blood test.");
}

// 7. Absolute risk (native chart) -------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "What this means in absolute terms", "Strokes per 1,000 women per year after the Hb test (31,400 women followed)");
  s.addChart(pres.charts.BAR, [{ name: "Rate", labels: ["None", "Mild", "Moderate", "Severe"], values: [4.32, 4.82, 7.43, 5.95] }], {
    x: M, y: 1.9, w: 7.2, h: 4.8, barDir: "col", chartColors: ["C9A3A9", "B5646F", C.crimson, "5E1220"], varyColors: true,
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.0", dataLabelFontSize: 14, dataLabelColor: C.ink,
    catAxisLabelColor: C.muted, valAxisLabelColor: C.muted, catAxisLabelFontSize: 13, valAxisMinVal: 0, valAxisMaxVal: 9,
    valAxisMajorUnit: 3, valGridLine: { color: "E5E5E5", size: 0.5 }, catGridLine: { style: "none" }, showLegend: false,
    barGapWidthPct: 60, valAxisTitle: "Strokes per 1,000 person-years", showValAxisTitle: true, valAxisTitleFontSize: 12,
    valAxisTitleColor: C.muted,
  });
  const x2 = 8.3, w2 = W - M - x2;
  card(s, x2, 1.95, w2, 2.1, C.light);
  s.addText("5-year risk of stroke", { x: x2 + 0.3, y: 2.1, w: w2 - 0.6, h: 0.4, fontFace: BF, fontSize: 14, color: C.muted, margin: 0, isTextBox: true });
  s.addText([{ text: "2.1%", options: { color: C.navy } }, { text: "  →  ", options: { color: C.muted } }, { text: "3.6%", options: { color: C.crimson } }],
    { x: x2 + 0.3, y: 2.5, w: w2 - 0.4, h: 0.9, fontFace: HF, fontSize: 32, bold: true, margin: 0, isTextBox: true });
  s.addText("no anemia vs moderate anemia", { x: x2 + 0.3, y: 3.4, w: w2 - 0.6, h: 0.4, fontFace: BF, fontSize: 14, color: C.ink, margin: 0, isTextBox: true });
  s.addText([
    { text: "Adjusted rate ratio, moderate anemia: ", options: { bold: true } }, { text: "1.60 (1.10–2.33)", options: { breakLine: true } },
    { text: " ", options: { breakLine: true } },
    { text: "Severe anemia had only 9 strokes, so its rate is imprecise, and these women more often died of other causes first." },
  ], { x: x2, y: 4.35, w: w2, h: 2.3, fontFace: BF, fontSize: 14, color: C.ink, margin: 0, valign: "top", isTextBox: true });
  s.addNotes("In absolute terms: about 4 strokes per 1,000 women per year without anemia, about 7 with moderate anemia. Over 5 years that is roughly 2 percent versus 4 percent. The severe bar looks lower but rests on 9 strokes, and these women died of other causes 5 times as often, which removes them before a stroke can happen.");
}

// 8. Stress tests ----------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Does it survive the obvious objections?", "Moderate anemia vs no anemia, odds ratio (95% CI)");
  const tests = [
    ["“The stroke caused the anemia”", "Hb measured ≥30 days before the stroke", "1.83 (1.24–2.71)"],
    ["", "Followed forward from the Hb test (rate ratio)", "1.60 (1.10–2.33)"],
    ["“It's kidney disease, sickle cell…”", "Adjusted for anemia-causing conditions", "1.49 (1.18–1.89)"],
    ["", "Those women excluded entirely", "1.46 (1.06–2.00)"],
    ["“They were already very ill”", "Cancer, heart failure, kidney, liver disease excluded", "1.64 (1.23–2.19)"],
    ["“Sicker women get tested more”", "Adjusted for number of blood tests", "1.38 (1.09–1.74)"],
  ];
  const y0 = 1.9, rh = 0.66;
  tests.forEach(([q, t, v], i) => {
    const y = y0 + i * rh;
    if (i % 2 === 0) s.addShape(pres.shapes.RECTANGLE, { x: M, y, w: W - 2 * M, h: rh * 2 - 0.08, fill: { color: C.tint }, line: { color: C.tint } });
    if (q) s.addText(q, { x: M + 0.25, y, w: 3.8, h: rh * 2 - 0.08, fontFace: HF, fontSize: 16, bold: true, italic: true,
      color: C.dark, valign: "middle", margin: 0, isTextBox: true });
    s.addText(t, { x: 4.75, y, w: 5.1, h: rh, fontFace: BF, fontSize: 15, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
    s.addText(v, { x: 9.9, y, w: 2.15, h: rh, fontFace: BF, fontSize: 15, bold: true, color: C.crimson, align: "right", valign: "middle", margin: 0, isTextBox: true });
    s.addShape(pres.shapes.OVAL, { x: 12.3, y: y + rh / 2 - 0.13, w: 0.26, h: 0.26, fill: { color: C.good }, line: { color: C.good } });
  });
  s.addText("Main result: 1.77 (1.41–2.21). Every estimate stays above 1. An unmeasured confounder would need a ~2.9-fold link with both anemia and stroke to explain it away (E-value).",
    { x: M, y: 6.0, w: W - 2 * M, h: 0.8, fontFace: BF, fontSize: 14, color: C.muted, margin: 0, isTextBox: true });
  s.addNotes("We tested the obvious objections. Reverse causation: the link holds when hemoglobin was measured well before the stroke and when we follow women forward. Confounding by diseases that cause anemia: it holds after adjusting for them and after excluding those women. Very ill women: holds. More testing in sicker women: it weakens a little but holds. An unmeasured confounder would need to nearly triple both anemia and stroke risk to explain it away.");
}

// 9. Repeated blood tests --------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Following hemoglobin over time", "420,102 blood tests; median 6 per woman over 5 years");
  const iw = 8.4, ih = iw / 2.071;
  s.addImage({ path: FIG + "fig18_p2_longitudinal.png", x: M, y: 1.85, w: iw, h: ih });
  const x2 = M + iw + 0.35, w2 = W - M - x2;
  const pts = [
    ["Updated at every test", "Moderate: 1.85× the stroke rate (1.28–2.67)"],
    ["Anemia often persists", "Of those retested, 65% (moderate) and 71% (severe) still anemic a year later"],
    ["Persistent vs resolved", "Persistent 1.45× (1.02–2.06); resolved 1.09× (no clear excess)"],
  ];
  pts.forEach(([h, t], i) => {
    const y = 1.85 + i * 1.35;
    badge(s, x2, y, i + 1);
    s.addText(h, { x: x2 + 0.7, y: y - 0.02, w: w2 - 0.7, h: 0.45, fontFace: HF, fontSize: 14, bold: true, color: C.dark, margin: 0, isTextBox: true });
    s.addText(t, { x: x2 + 0.7, y: y + 0.42, w: w2 - 0.7, h: 0.9, fontFace: BF, fontSize: 13.5, color: C.ink, margin: 0, valign: "top", isTextBox: true });
  });
  card(s, M, 6.1, W - 2 * M, 0.85, C.tint);
  s.addText([{ text: "Honest caveat: ", options: { bold: true, color: C.crimson } },
    { text: "when the hemoglobin from the 30 days before a stroke is excluded, the estimate falls to 1.28 (0.84–1.96). Part of the signal comes from illness just before the stroke." }],
    { x: M + 0.3, y: 6.1, w: W - 2 * M - 0.6, h: 0.85, fontFace: BF, fontSize: 14, color: C.ink, valign: "middle", margin: 0, isTextBox: true });
  s.addNotes("Most women had several blood tests, so we could track anemia over time. When anemia is updated at every test, moderate anemia still carries about 1.85 times the stroke rate. Most anemia persists. Anemia that persisted at one year carried more risk than anemia that resolved. The caveat: if we ignore hemoglobin from the month before a stroke, the estimate drops to 1.28 and is no longer significant, so part of the signal reflects illness just before the event. The true long-term effect is probably smaller than the headline figure.");
}

// 10. Type of anemia --------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "Which anemia matters?", "Red-cell pattern from the same blood count (MCV and RDW), odds ratio vs no anemia");
  const iw = 7.6, ih = iw / 1.863;
  s.addImage({ path: FIG + "fig15_p2_anaemia_pattern.png", x: M, y: 1.85, w: iw, h: ih });
  const x2 = M + iw + 0.4, w2 = W - M - x2;
  const rows = [["Pattern", "OR (95% CI)"], ["Iron-deficiency pattern", "1.17 (0.93–1.48)"], ["Normal cells, normal RDW", "1.30 (1.09–1.55)"],
    ["Normal cells, high RDW", "1.73 (1.40–2.14)"], ["Large cells (19 strokes)", "3.01 (1.74–5.23)"]];
  s.addTable(rows.map((r, i) => r.map((v, j) => ({ text: v, options: { bold: i === 0 || (i >= 3 && j === 1), color: i === 0 ? C.light : (i >= 3 && j === 1 ? C.crimson : C.ink),
    fill: { color: i === 0 ? C.crimson : (i % 2 ? "FFFFFF" : C.tint) } } }))),
    { x: x2, y: 1.95, w: w2, colW: [w2 * 0.56, w2 * 0.44], rowH: 0.5, fontFace: BF, fontSize: 13, border: { type: "solid", color: C.line, pt: 0.5 } });
  s.addText("The classic iron-deficiency pattern from menstrual loss carried the least risk. Anemia with other features may point to underlying illness.",
    { x: x2, y: 4.8, w: w2, h: 1.6, fontFace: BF, fontSize: 14, color: C.ink, margin: 0, valign: "top", isTextBox: true });
  footer(s, "Exploratory. Adjusted for main covariates; RDW-CV >14.5% = high; thalassemia-trait pattern too few strokes to estimate.");
  s.addNotes("Using red cell size and RDW from the same blood count: the classic iron-deficiency pattern carried the least risk. The strongest links were with normal-sized cells and high RDW, and with large cells, though that group had only 19 strokes. This suggests anemia that is not simply from menstrual blood loss may signal other illness. This part is exploratory.");
}

// 11. Limitations -----------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  title(s, "What we can and cannot say");
  const left = ["Large cohort: 39,807 women, 1,538 strokes/TIAs", "Clear dose-response with hemoglobin", "Consistent across multiple sensitivity analyses",
    "Repeated blood tests, not a single value"];
  const right = ["Association, not causation", "Part of the signal is from Hb shortly before the stroke", "Single center; strokes from records, not adjudicated",
    "Smoking often unrecorded; outpatient iron not captured"];
  [[left, "Strengths", C.good], [right, "Limitations", C.crimson]].forEach(([list, hdr, col], k) => {
    const x = M + k * 6.2;
    card(s, x, 1.7, 5.9, 4.9, k === 0 ? "EDF5F0" : C.tint);
    s.addText(hdr, { x: x + 0.4, y: 1.95, w: 5.1, h: 0.55, fontFace: HF, fontSize: 22, bold: true, color: col, margin: 0, isTextBox: true });
    s.addText(list.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < list.length - 1 } })),
      { x: x + 0.4, y: 2.65, w: 5.1, h: 3.7, fontFace: BF, fontSize: 16, color: C.ink, margin: 0, valign: "top", paraSpaceAfter: 12, isTextBox: true });
  });
  s.addNotes("Strengths: size, a clear dose-response, consistency across many analyses, and repeated blood tests. Limitations: this is an association, not proof of cause. Part of the signal reflects hemoglobin measured just before the stroke. It is one center, strokes come from records, smoking is often missing, and we only see iron given in hospital.");
}

// 12. Take-home -------------------------------------------------------------------------------------------
{
  const s = pres.addSlide();
  s.background = { color: C.dark };
  s.addText("Take-home messages", { x: M, y: 0.6, w: W - 2 * M, h: 0.8, fontFace: HF, fontSize: 36, bold: true, color: C.light, margin: 0, isTextBox: true });
  const msgs = [["Moderate or severe anemia", "is linked to roughly 1.5–2 times the risk of stroke in women with benign uterine disease."],
    ["The risk rises steadily", "below a hemoglobin of 13 g/dL, independent of standard vascular risk factors."],
    ["Hemoglobin is already measured", "at gynecological diagnosis: a free flag to prompt vascular risk assessment."]];
  msgs.forEach(([h, t], i) => {
    const y = 1.8 + i * 1.5;
    s.addShape(pres.shapes.OVAL, { x: M, y, w: 0.75, h: 0.75, fill: { color: C.crimson }, line: { color: C.crimson } });
    s.addText(String(i + 1), { x: M, y, w: 0.75, h: 0.75, fontFace: HF, fontSize: 24, bold: true, color: C.light, align: "center", valign: "middle", margin: 0, isTextBox: true });
    s.addText([{ text: h + " ", options: { bold: true, color: C.light } }, { text: t, options: { color: "E3D3D6" } }],
      { x: M + 1.1, y: y - 0.1, w: 10.8, h: 1.0, fontFace: BF, fontSize: 20, valign: "middle", margin: 0, isTextBox: true });
  });
  s.addText("Next: does correcting anemia change stroke risk? This needs a prospective study.", { x: M, y: 6.35, w: W - 2 * M, h: 0.5,
    fontFace: BF, fontSize: 16, italic: true, color: "D9A9B0", margin: 0, isTextBox: true });
  s.addNotes("Three messages. Moderate or severe anemia goes with roughly 1.5 to 2 times the stroke risk. The risk climbs below 13 grams. And hemoglobin is already measured, so it is a free prompt to check vascular risk. Whether treating anemia lowers stroke risk is the next question and needs a prospective study.");
}

pres.writeFile({ fileName: OUT }).then(f => console.log("wrote", f));
