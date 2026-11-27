# lunar_setup.py -- packages and output helpers for the lunar regolith-TES sCO2 notebook
# Loaded by Step 1 of the notebook (keep this file next to the notebook).
# Defines: np, CP (CoolProp), plt, PdfPages, io, textwrap,
#          SHOW_IN_CELLS, OUT_PDF, RUN_LOG, REPORT, say(), rep()
import io, textwrap                                   # figures in memory, text wrapping (PDF)
import numpy as np                                    # arrays and math
import CoolProp.CoolProp as CP                        # sCO2 and water properties
import matplotlib.pyplot as plt                       # figures
from matplotlib.backends.backend_pdf import PdfPages  # results PDF

# Output: the cells print nothing (unless SHOW_IN_CELLS = True). Every message, the
# report, all figures and the references are written to ONE PDF at the end.
SHOW_IN_CELLS = False
OUT_PDF = 'Lunar_TES_sCO2_results.pdf'
RUN_LOG, REPORT = [], []

def say(*args, sep=' '):
    """Message from a cell: kept for the PDF run log (printed only if SHOW_IN_CELLS)."""
    text = sep.join(str(a) for a in args)
    RUN_LOG.extend(text.split('\n'))
    if SHOW_IN_CELLS:
        print(text)

def rep(*args, sep=' '):
    """Line of the report: kept for the PDF (printed only if SHOW_IN_CELLS)."""
    text = sep.join(str(a) for a in args)
    REPORT.extend(text.split('\n'))
    if SHOW_IN_CELLS:
        print(text)
