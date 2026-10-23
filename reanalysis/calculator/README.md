# MED-HACK Score calculator

- `bud_stroke_calculator.html`: the page published as a Claude artifact.
- `standalone/index.html`: the same page as a complete, self-contained HTML file. It needs no server code or database and runs entirely in the browser, so it can be hosted anywhere static (GitHub Pages, a university web server, Netlify) or opened directly from disk.
- The model values come from `../output/risk_model_spec.json` (produced by `R/15_risk_model.R`). The calculator's arithmetic was checked against R for test profiles: exact match.
- Status: research use only; not externally validated. Update the citation line when the paper is published.
