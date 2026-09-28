"""Run the complete analysis end to end:

    python -m analysis.run_all

Reads data/stroke_analysis_master.csv and data/stroke.xlsx (read-only) and
writes everything to outputs/.
"""
import time

from . import data_prep, excel_out, figures, manuscript, paper1, paper2, paper2_ext, paper3, report, surgery
from .utils import OUT_DIR, RESULTS, log


def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df, elig = data_prep.run()
    paper1.run(df, elig)
    paper2.run(elig)
    paper2_ext.run(elig)
    paper3.run(elig)
    surgery.run(df)
    figures.run()
    xlsx = excel_out.write()
    log("Run", f"Workbook written: {xlsx.name}; figures: {sum(len(v) for v in RESULTS['figures'].values())} files; "
               f"elapsed {time.time() - t0:.0f} s")
    report.run()
    manuscript.run()
    print("Done.")


if __name__ == "__main__":
    main()
