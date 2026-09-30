"""Run the complete analysis end to end:

    python -m analysis.run_all

Reads data/stroke_analysis_master.csv and data/stroke.xlsx (read-only) and
writes everything to outputs/.
"""
import time

from . import (data_prep, excel_out, figures, manuscript, manuscript_p3, master_export, paper1, paper2, paper2_ext,
               paper2_long, paper3, paper4, report, surgery)
from .utils import OUT_DIR, RESULTS, log


def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df, elig = data_prep.run()
    paper1.run(df, elig)
    paper2.run(elig)
    d2 = paper2_ext.run(elig)
    paper2_long.run(d2)
    paper3.run(elig)
    paper4.run(elig)
    surgery.run(df)
    figures.run()
    xlsx = excel_out.write()
    log("Run", f"Workbook written: {xlsx.name}; figures: {sum(len(v) for v in RESULTS['figures'].values())} files; "
               f"elapsed {time.time() - t0:.0f} s")
    report.run()
    manuscript.run()
    manuscript_p3.run()
    master_export.write(df, d2, RESULTS["_p2x_t"], data_prep.raw_columns())
    print("Done.")


if __name__ == "__main__":
    main()
