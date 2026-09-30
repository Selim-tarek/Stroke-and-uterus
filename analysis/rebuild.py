"""Rebuild figures, report and manuscripts from the cached state of the last full run (no model fitting):

    python -m analysis.rebuild
"""
import pickle

from . import figures, manuscript, manuscript_p3, report, utils
from .utils import OUT_DIR


def main():
    with open(OUT_DIR / "_state.pkl", "rb") as fh:
        state = pickle.load(fh)
    for k, v in state.items():
        getattr(utils, k).clear()
        getattr(utils, k).update(v) if isinstance(v, dict) else getattr(utils, k).extend(v)
    figures.run()
    report.run()
    manuscript.run()
    manuscript_p3.run()
    print("Rebuilt.")


if __name__ == "__main__":
    main()
