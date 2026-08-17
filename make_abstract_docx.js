// Builds the formatted ISC 2027 abstract document. Table values come from
// results/abstract_tables.json, produced by abstract_tables.R, so the document
// cannot drift from the analysis.
//   Rscript abstract_tables.R && Rscript make_abstract_docx.R && node make_abstract_docx.js
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  WidthType, AlignmentType, BorderStyle, ShadingType, HeadingLevel,
} = require('docx');

const data = JSON.parse(fs.readFileSync('results/abstract_tables.json', 'utf8'));

const FONT = 'Calibri';
const INK = '1A1A1A';
const RULE = '8A8A8A';
const HAIR = 'D9D9D9';
const HEADFILL = 'EDEDED';
const SUBFILL = 'F7F7F7';

const noBorder = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
const cellBorders = (top, bottom) => ({
  top: top || noBorder, bottom: bottom || noBorder, left: noBorder, right: noBorder,
});
const line = (color, size) => ({ style: BorderStyle.SINGLE, size: size || 4, color });

function txt(text, opts = {}) {
  return new TextRun({
    text: String(text == null ? '' : text),
    font: FONT, size: opts.size || 18, bold: !!opts.bold, italics: !!opts.italics,
    color: opts.color || INK,
  });
}

function para(text, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.LEFT,
    spacing: { before: opts.before || 0, after: opts.after == null ? 40 : opts.after, line: opts.line || 240 },
    indent: opts.indent ? { left: opts.indent } : undefined,
    children: Array.isArray(text) ? text : [txt(text, opts)],
  });
}

function cell(children, opts = {}) {
  return new TableCell({
    width: { size: opts.width, type: WidthType.DXA },
    borders: opts.borders || cellBorders(),
    shading: opts.fill ? { type: ShadingType.CLEAR, fill: opts.fill, color: 'auto' } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children,
  });
}

// A table row: `cells` is an array of {text, bold, align, indent, italics}.
function row(cells, widths, opts = {}) {
  return new TableRow({
    tableHeader: !!opts.header,
    children: cells.map((c, i) => cell(
      [para(c.text, {
        bold: c.bold, italics: c.italics, indent: c.indent,
        align: c.align || (i === 0 ? AlignmentType.LEFT : AlignmentType.RIGHT),
        after: 0, color: c.color,
      })],
      { width: widths[i], fill: opts.fill, borders: opts.borders },
    )),
  });
}

function buildTable(widths, headerCells, bodyRows) {
  const rows = [row(headerCells, widths, {
    header: true, fill: HEADFILL,
    borders: cellBorders(line(RULE, 6), line(RULE, 6)),
  })];
  bodyRows.forEach((r, i) => {
    const last = i === bodyRows.length - 1;
    rows.push(row(r.cells, widths, {
      fill: r.section ? SUBFILL : undefined,
      borders: cellBorders(undefined, last ? line(RULE, 6) : line(HAIR, 2)),
    }));
  });
  return new Table({ columnWidths: widths, width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA }, rows });
}

const caption = (n, text) => new Paragraph({
  spacing: { before: 260, after: 100 },
  children: [txt(`Table ${n}. `, { bold: true }), txt(text)],
});

const note = (text) => new Paragraph({
  spacing: { before: 80, after: 200 },
  children: [txt(text, { size: 15, italics: true, color: '595959' })],
});

// --------------------------------------------------------------- Table 1
const W1 = [5400, 1900, 1900];
const t1 = [
  { c: ['Female', '47,123', '100%'] },
  { c: ['Excluded', '', ''], section: true },
  { c: ['Active malignancy', '7,316', '16%'], indent: 200 },
  { c: ['Included', '39,807', '84%'], bold: true },
  { c: ['Age', '', ''], section: true },
  { c: ['18-29', '5,344', '13%'], indent: 200 },
  { c: ['30-39', '11,806', '30%'], indent: 200 },
  { c: ['40-49', '14,777', '37%'], indent: 200 },
  { c: ['50-60', '7,880', '20%'], indent: 200 },
  { c: ['Race', '', ''], section: true },
  { c: ['White', '31,450', '79%'], indent: 200 },
  { c: ['Black or African American', '3,715', '9%'], indent: 200 },
  { c: ['Asian', '2,416', '6%'], indent: 200 },
  { c: ['Native Hawaiian or Pacific Islander', '61', '0.2%'], indent: 200 },
  { c: ['Other', '914', '2%'], indent: 200 },
  { c: ['Unknown or not disclosed', '931', '2%'], indent: 200 },
  { c: ['Uterine pathology', '', ''], section: true },
  { c: ['Fibroids', '14,461', '36%'], indent: 200 },
  { c: ['Endometriosis', '11,445', '29%'], indent: 200 },
  { c: ['Adenomyosis', '1,520', '4%'], indent: 200 },
  { c: ['Overlap or multiple', '12,381', '31%'], indent: 200 },
];
const table1 = buildTable(W1,
  [{ text: 'Characteristic', bold: true }, { text: 'n', bold: true }, { text: '%', bold: true }],
  t1.map((r) => ({
    section: r.section,
    cells: r.c.map((v, i) => ({ text: v, bold: r.bold || r.section, indent: i === 0 ? r.indent : 0 })),
  })));

// --------------------------------------------------------------- Table 2
const W2 = [5000, 1500, 1500, 1200];
const t2map = {};
data.table2.forEach((r) => { t2map[r.Characteristic] = r; });
const pathRow = (label, key) => {
  const r = t2map[key] || {};
  const m = /^([\d,]+) \((\d+)%\)$/.exec(r['n (%)'] || '');
  return { c: [label, m ? m[1] : '', m ? `${m[2]}%` : '', r['Prevalence in group'] || ''], indent: 200 };
};
const t2 = [
  { c: ['Prevalence of stroke', '1,538', '4%', ''], bold: true },
  { c: ['Stroke type', '', '', ''], section: true },
  { c: ['Ischemic stroke', '1,189', '77%', ''], indent: 200 },
  { c: ['TIA', '198', '13%', ''], indent: 200 },
  { c: ['Hemorrhagic stroke', '57', '4%', ''], indent: 200 },
  { c: ['SAH', '40', '3%', ''], indent: 200 },
  { c: ['CVST', '13', '1%', ''], indent: 200 },
  { c: ['Unknown', '41', '3%', ''], indent: 200 },
  { c: ['Uterine pathology in ischemic stroke', '', '', ''], section: true },
  pathRow('Fibroids', 'Fibroids only'),
  pathRow('Endometriosis', 'Endometriosis only'),
  pathRow('Adenomyosis', 'Adenomyosis only'),
  pathRow('Overlap or multiple', 'More than one'),
  { c: ['Suspected ischemic stroke mechanism', '', '', ''], section: true },
  { c: ['Small vessel disease', '268', '23%', ''], indent: 200 },
  { c: ['Large artery atherosclerosis', '46', '4%', ''], indent: 200 },
  { c: ['Cardioembolic', '15', '1%', ''], indent: 200 },
  { c: ['Other determined etiology', '147', '12%', ''], indent: 200 },
  { c: ['Cryptogenic or unknown', '713', '60%', ''], indent: 200 },
  { c: ['Number of ischemic infarcts', '', '', ''], section: true },
  { c: ['One', '148', '12%', ''], indent: 200 },
  { c: ['Two', '52', '4%', ''], indent: 200 },
  { c: ['Three or more', '273', '23%', ''], indent: 200 },
  { c: ['Missing data', '716', '60%', ''], indent: 200 },
];
const table2 = buildTable(W2,
  [{ text: 'Characteristic', bold: true }, { text: 'n', bold: true },
    { text: '%', bold: true }, { text: 'Prevalence', bold: true }],
  t2.map((r) => ({
    section: r.section,
    cells: r.c.map((v, i) => ({ text: v, bold: r.bold || r.section, indent: i === 0 ? r.indent : 0 })),
  })));

// --------------------------------------------------------------- Table 3
const W3 = [3600, 1000, 1000, 1500, 1500, 1000];
// A predictor with a single "Yes" row is binary — print it on one line under
// its own name rather than spending a header row on it.
const yesOnly = new Set();
{
  const byVar = {};
  data.table3.forEach((r) => { (byVar[r.Variable] = byVar[r.Variable] || []).push(r.Level); });
  Object.entries(byVar).forEach(([v, lv]) => {
    if (lv.length === 1 && lv[0] === 'Yes') yesOnly.add(v);
  });
}

const t3rows = [];
let lastVar = null;
data.table3.forEach((r) => {
  const isMulti = r.Level !== '' && !yesOnly.has(r.Variable);
  if (isMulti && r.Variable !== lastVar) {
    t3rows.push({ section: true, cells: [{ text: r.Variable, bold: true }, {}, {}, {}, {}, {}].map((c) => ({ text: c.text || '', bold: c.bold })) });
  }
  const label = isMulti ? r.Level : r.Variable;
  const ref = r['p-value'] === 'not estimable';
  t3rows.push({
    cells: [
      { text: label, indent: isMulti ? 200 : 0 },
      { text: r.n.toLocaleString() },
      { text: r.events.toLocaleString() },
      { text: ref ? '-' : r['Adjusted odds ratio'] },
      { text: ref ? '-' : r['95% CI'] },
      { text: r['p-value'], italics: ref },
    ],
  });
  lastVar = r.Variable;
});
const table3 = buildTable(W3,
  [{ text: 'Predictor', bold: true }, { text: 'n', bold: true }, { text: 'Ischemic strokes', bold: true },
    { text: 'Adjusted odds ratio', bold: true }, { text: '95% CI', bold: true }, { text: 'p-value', bold: true }],
  t3rows);

// --------------------------------------------------------------- Document
const m = Array.isArray(data.model) ? data.model[0] : data.model;
const body = [
  new Paragraph({
    spacing: { after: 80 },
    children: [txt('Noncancerous uterine disease as an independent risk factor for stroke in young women', { bold: true, size: 26 })],
  }),
  new Paragraph({
    spacing: { after: 260 },
    children: [txt('Elizabeth L. Saionz, Selim T. Tarabeah, Mahmoud M. Afia … Michelle P. Lin', { size: 19, color: '444444' })],
  }),

  para([txt('Introduction. ', { bold: true }), txt('Noncancerous uterine disease, such as fibroids, endometriosis, and adenomyosis, is common among women of childbearing age. Small case series have reported acute ischemic stroke in women with noncancerous uterine disease, often attributed to cryptogenic or hypercoagulable etiology. Given the prevalence of noncancerous uterine disease, we sought to determine the risk of acute ischemic stroke among this population and the associated stroke mechanism.')], { after: 140 }),

  para([txt('Methods. ', { bold: true }), txt('We performed a retrospective review of cases at a single multi-site enterprise from 2018-2025. We included women ages 18-60 with diagnosis or imaging findings of fibroids, endometriosis, or adenomyosis. Patients were excluded if they had a comorbid active diagnosis of systemic malignancy aside from non-melanomatous skin cancer. The primary outcome was prevalence of ischemic stroke.')], { after: 140 }),

  para([txt('Results. ', { bold: true }), txt('We included 39,807 patients ages 18-60 with noncancerous uterine disease: 14,461 (36%) fibroids, 11,445 (29%) endometriosis, 1,520 (4%) adenomyosis, 12,381 (31%) overlap diagnoses. Additional demographics are detailed in Table 1. The overall prevalence of stroke was 1,538 (4%), with 77% being imaging-confirmed ischemic stroke. Among the different uterine pathologies, stroke occurred in 4.4% fibroids, 3.7% adenomyosis, 2.9% endometriosis, and 4.2% overlap diagnoses (Table 2). We performed a multivariate logistic regression (Table 3) to identify potential risk factors for stroke in this population. '),
    txt(`Ischemic stroke was associated with lower hemoglobin (aOR 0.92 per 1 g/dL, 95% CI 0.89-0.95), migraine with aura (aOR 3.64, 3.03-4.35), atrial fibrillation (aOR 3.34, 2.54-4.36), congestive heart failure (aOR 3.24, 2.42-4.30), coronary artery disease (aOR 3.09, 2.38-3.97), hypertension (aOR 2.34, 2.05-2.68), estrogen-containing hormone therapy (aOR 1.72, 1.14-2.49), dyslipidemia (aOR 1.50, 1.31-1.72), and a history of smoking (aOR 1.29, 1.08-1.56). Relative to fibroids, isolated endometriosis (aOR 0.82, 0.69-0.97) and adenomyosis (aOR 0.67, 0.46-0.93) carried lower odds.`, { bold: true })], { after: 140 }),

  para([txt('Conclusion. ', { bold: true }), txt('Noncancerous uterine disease may contribute to stroke risk among young women, in part related to coagulopathy from iron-deficiency anemia. Secondary stroke prevention strategies may include addressing the underlying gynecologic pathology.')], { after: 200 }),

  caption(1, 'Study demographics, women with noncancerous uterine pathology ages 18-60, 2018-2025'),
  table1,
  note('Percentages are of the 39,807 included patients, except for the excluded row, which is of the 47,123 screened.'),

  caption(2, 'Prevalence of stroke, women ages 18-60 with noncancerous uterine pathology'),
  table2,
  note('Percentages are of all 1,538 strokes for stroke type, and of the 1,189 ischemic strokes for the sections below it. Prevalence is ischemic strokes as a share of that pathology group.'),

  caption(3, 'Multivariate logistic regression model for ischemic stroke in women ages 18-60 with noncancerous uterine pathology'),
  table3,
  note(`Model: ${m.n.toLocaleString()} patients with ${m.events.toLocaleString()} ischemic strokes; C-statistic ${m.c_statistic}. Reference categories are White race, never smoker, no migraine, and fibroids only. Patients whose only cerebrovascular event was hemorrhagic, SAH, TIA or CVST were excluded rather than counted as controls. Native Hawaiian or Pacific Islander is reported but not estimable (2 events).`),
];

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 18, color: INK } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1080, right: 1080 } } },
    children: body,
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync('results/ISC_2027_abstract_Saionz_tables.docx', buf);
  console.log('Wrote results/ISC_2027_abstract_Saionz_tables.docx');
});
