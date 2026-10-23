// Audit Q7: evaluate the calculator's own MODEL block in node, then drive the rendered page (local file and public URL).
// usage: node calc_test.js <profiles.csv> <standalone/index.html> <out_dir> <public_url>
const fs = require("fs"), vm = require("vm"), path = require("path");
const [csvPath, htmlPath, outDir, url] = process.argv.slice(2);
const lines = fs.readFileSync(csvPath, "utf8").trim().split(/\r?\n/);
const head = lines[0].split(",");
const P = lines.slice(1).map(l => { const v = l.split(","); const o = {}; head.forEach((h, i) => o[h] = isNaN(+v[i]) ? v[i] : +v[i]); return o; });
const html = fs.readFileSync(htmlPath, "utf8");
const block = html.slice(html.indexOf("// MODEL-START"), html.indexOf("// MODEL-END"));
const ctx = {}; vm.createContext(ctx); vm.runInContext(block + "\nthis.out = { score, scoreRisk, modelRisk };", ctx);
const { score, scoreRisk, modelRisk } = ctx.out;
const BIN = ["migraine_any", "estrogen", "dm", "htn", "haemoglobinopathy", "cad", "ckd_pre"];
const res = P.map(p => { const s = score(p).total; return { id: p.id, js_points: s, js_points_risk: scoreRisk(s), js_model: modelRisk(p) }; });

async function drive(page, target, tag) {
  await page.goto(target, { waitUntil: "load" });
  for (const [i, p] of P.entries()) {
    for (const k of BIN) await page.click(`button[data-k="${k}"][data-v="${p[k]}"]`);
    await page.selectOption("#hbUnit", "gdl");
    await page.fill("#age", String(p.age)); await page.fill("#hb", String(p.hb));
    const pts = await page.textContent("#pts"), risk = await page.textContent("#risk"), interp = await page.textContent("#interp");
    const m = interp.match(/Exact model risk\s*([\d.]+)%/);
    res[i][`web_${tag}_points`] = +pts; res[i][`web_${tag}_points_risk`] = parseFloat(risk); res[i][`web_${tag}_model`] = m ? +m[1] : NaN;
    if (p.id === 3) for (const [u, v] of [["gl", "94"], ["mmol", "5.835"]]) {
      await page.selectOption("#hbUnit", u); await page.fill("#hb", v);
      res[i][`web_${tag}_${u}`] = `${await page.textContent("#pts")} points, ${(await page.textContent("#interp")).match(/Exact model risk\s*([\d.]+)%/)?.[1]}% model risk`;
      await page.selectOption("#hbUnit", "gdl");
    }
  }
}
(async () => {
  const { chromium } = require("playwright");
  const browser = await chromium.launch({ executablePath: fs.existsSync("/opt/pw-browsers/chromium") ? undefined : undefined });
  const page = await browser.newPage();
  await drive(page, "file://" + path.resolve(htmlPath), "local");
  let pub = "public page reached";
  try { await drive(page, url, "public"); } catch (e) { pub = "public page NOT tested: " + e.message.split("\n")[0]; }
  await browser.close();
  const cols = Object.keys(res.reduce((a, r) => Object.assign(a, r), {}));
  fs.writeFileSync(path.join(outDir, "A05_js_web_output.csv"), [cols.join(","), ...res.map(r => cols.map(c => r[c] ?? "").map(String).map(s => s.includes(",") ? `"${s}"` : s).join(","))].join("\n"));
  console.log(`JS evaluated for ${res.length} profiles; local page driven; ${pub}`);
})();
