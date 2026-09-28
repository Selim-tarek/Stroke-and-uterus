// Render outputs/anaemia_manuscript.json (built by analysis/manuscript.py) to .docx
// usage: node analysis/build_docx.js in.json out.docx
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell, WidthType, ShadingType,
  ImageRun, AlignmentType, PageBreak, Footer, PageNumber, BorderStyle, TableLayoutType,
} = require("docx");

const [, , inPath, outPath] = process.argv;
const blocks = JSON.parse(fs.readFileSync(inPath, "utf8"));
const PAGE_W = 12240, MARGIN = 1440, CONTENT_W = PAGE_W - 2 * MARGIN; // US Letter, 1" margins
const FONT = "Arial";

function runs(text, opts = {}) {
  // **bold** markup -> bold runs
  const parts = String(text).split(/(\*\*[^*]+\*\*)/g).filter((s) => s.length);
  return parts.map((p) =>
    p.startsWith("**") ? new TextRun({ text: p.slice(2, -2), bold: true, font: FONT, ...opts })
      : new TextRun({ text: p, font: FONT, ...opts }));
}

function pngSize(path) {
  const b = fs.readFileSync(path);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20), data: b };
}

const border = { style: BorderStyle.SINGLE, size: 4, color: "BFBFBF" };
const borders = { top: border, bottom: border, left: border, right: border };

function table(b) {
  const n = b.columns.length;
  let widths = b.widths && b.widths.length === n ? b.widths.slice() : Array(n).fill(1);
  const sum = widths.reduce((a, c) => a + c, 0);
  widths = widths.map((w) => Math.floor((w / sum) * CONTENT_W));
  widths[n - 1] += CONTENT_W - widths.reduce((a, c) => a + c, 0);
  const size = (b.font || 8) * 2;
  const cell = (txt, i, header) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    borders,
    shading: header ? { type: ShadingType.CLEAR, color: "auto", fill: "D9E2F3" } : undefined,
    margins: { top: 40, bottom: 40, left: 60, right: 60 },
    children: [new Paragraph({ children: runs(txt, { size, bold: header }) })],
  });
  const rows = [new TableRow({ tableHeader: true, children: b.columns.map((c, i) => cell(c, i, true)) })];
  for (const r of b.rows) rows.push(new TableRow({ children: r.map((v, i) => cell(v, i, false)) }));
  const out = [new Paragraph({ spacing: { before: 200, after: 100 }, children: runs(b.title, { bold: true, size: 20 }) }),
    new Table({ width: { size: CONTENT_W, type: WidthType.DXA }, columnWidths: widths, layout: TableLayoutType.FIXED, rows })];
  if (b.note) out.push(new Paragraph({ spacing: { before: 80, after: 200 }, children: runs(b.note, { italics: true, size: 16 }) }));
  return out;
}

const children = [];
for (const b of blocks) {
  if (b.t === "title") {
    children.push(new Paragraph({ heading: HeadingLevel.TITLE, spacing: { after: 300 }, children: runs(b.text, { size: 32, bold: true }) }));
  } else if (b.t === "h1") {
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 300, after: 120 }, children: runs(b.text, { size: 28, bold: true }) }));
  } else if (b.t === "h2") {
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 80 }, children: runs(b.text, { size: 24, bold: true }) }));
  } else if (b.t === "p") {
    children.push(new Paragraph({ spacing: { after: 140, line: 360 }, children: runs(b.text, { size: 22 }) }));
  } else if (b.t === "note") {
    children.push(new Paragraph({ spacing: { after: 140 }, children: runs(b.text, { size: 18, italics: true, color: "7F7F7F" }) }));
  } else if (b.t === "table") {
    children.push(...table(b));
  } else if (b.t === "figure") {
    const s = pngSize(b.path);
    const wpx = Math.round((b.width || 6.3) * 96);
    const hpx = Math.round((wpx * s.h) / s.w);
    children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 200 },
      children: [new ImageRun({ type: "png", data: s.data, transformation: { width: wpx, height: hpx } })] }));
    children.push(new Paragraph({ spacing: { after: 240 }, children: runs(b.caption, { size: 18 }) }));
  } else if (b.t === "pagebreak") {
    children.push(new Paragraph({ children: [new PageBreak()] }));
  }
}

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 22 } } } },
  sections: [{
    properties: { page: { size: { width: PAGE_W, height: 15840 }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18 })] })] }) },
    children,
  }],
});
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(outPath, buf); console.log("wrote", outPath); });
