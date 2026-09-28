"""Regenerate ska_hi_paper_figures.ipynb from the live plotting scripts.

    /home/geoff/gibbs_venv_312/bin/python tools_regen_figures_notebook.py

Each figure cell is `inspect.getsource` of the real function, so the notebook
is guaranteed to match the scripts at the moment it is generated. Private
helpers are resolved transitively and inlined into the cell that needs them.

**DESTRUCTIVE.** This overwrites the notebook, discarding any tuning done in
it. The flow is one-way, script -> notebook. Once colours or labels have been
tuned in the notebook, paste them back into `ska_hi_plots.py` /
`ska_hi_design.py` BEFORE running this again.

Close the notebook in VSCode first, or VSCode will write its in-memory copy
back over the new file.
"""
import ast, builtins, inspect, json, os, sys, textwrap
sys.path.insert(0, "/home/geoff/limTOD/examples/SKA")
import ska_hi_plots as P
import ska_hi_design as D
import ska_hi_kmodes as K

FIGS = [
    ("hi_ladder_geometry", D, "fig_ladder_geometry", "fig_ladder_geometry(lad)",
     "Sec. 2.3 -- Scan strategies",
     "The ladder as sky geometry. Below one beam the strips merge and concentrate "
     "measurement; at one beam and beyond they separate and buy area instead."),
    ("hi_priorvar", D, "fig_priorvar", "fig_priorvar(pv)",
     "Sec. 3 -- Map-making, `Prior choice'",
     "Loosening the prior raises N_eff and worsens the residual at the same time."),
    ("hi_ladder", D, "fig_ladder", "fig_ladder(lad, lad_long, floors)",
     "Sec. 4 -- How much sky does a strategy measure",
     "The drift's only geometry lever, out to N = 20, with the two measured "
     "floors overlaid on the middle panel."),
    ("hi_patch_maps", P, "fig_patch_maps", "fig_patch_maps(exp)",
     "Sec. 6 -- HI signal recovery",
     "The field and the pixel sets everything is quoted on: raster-only, the "
     "277 px all-channel intersection, and the 119 px interior."),
    ("hi_maps", P, "fig_hi_maps", "fig_hi_maps(exp)",
     "Sec. 6 -- HI signal recovery",
     "Map domain. True HI, both map-made responses on one scale, and the "
     "cleaned data needing a scale 82x wider. PLACEHOLDER -- see the TODO in "
     "the paper caption."),
    ("hi_transfer_function", P, "fig_transfer_function", "fig_transfer_function(exp)",
     "Sec. 6 -- HI signal recovery",
     "T(k_par) per mode count, +/-1 sigma over 20 mocks."),
    ("hi_ablation", P, "fig_ablation", "fig_ablation(abl)",
     "Sec. 6 -- HI signal recovery",
     "Floor vs total vs noise. Common resolution since 2026-08-26 (tag='cr_')."),
    ("hi_cleaning_depth", D, "fig_cleaning_depth", "fig_cleaning_depth()",
     "Sec. 6 -- HI signal recovery",
     "Three bands. The net never turns over: there is no optimal cleaning depth."),
    ("hi_realisations", D, "fig_realisations", "fig_realisations(vals)",
     "Sec. 6 -- HI signal recovery",
     "Nine HI realisations. The paired ratio is tighter than either absolute."),
    ("hi_chromatic", D, "fig_chromatic", "fig_chromatic(rk)",
     "Sec. 7 -- Why a different cleaning basis cannot help",
     "Rank, eigenvectors and the principal angles between floor and HI. "
     "Replaces the old hi_rank + hi_eigenvectors pair."),
    ("hi_design", D, "fig_design", "fig_design(d)",
     "Sec. 8 -- Survey design (headline)",
     "The trade space: two regimes, and the 9x step at matched mode density. "
     "Split from the lever chart 2026-09-28."),
    ("hi_levers", D, "fig_levers", "fig_levers(d)",
     "Sec. 8 -- Survey design",
     "Every lever as a factor on the floor, including the three that are "
     "1.00x exactly. The other half of the old two-panel figure."),
    ("hi_kmodes", K, "fig_kmodes", "fig_kmodes(kres)",
     "Sec. 8 -- Survey design",
     "tr(WA) decomposed by angular scale, and the sky area each strategy would "
     "need for sigma_P/P = 1 in auto and in cross-correlation."),
    ("hi_frequency_ratio", D, "fig_frequency_ratio",
     "fig_frequency_ratio(bands)", "Sec. 8 -- Survey design",
     "The quotient panel on its own. Same data as the third panel of "
     "hi_frequency_design, drawn through the same _freq_panel helper."),
    ("hi_frequency_design", D, "fig_frequency", "fig_frequency(bands)",
     "Sec. 8 -- Survey design",
     "Instrument, cosmology and their misleading quotient, over three bands."),
]

def md(src):  return {"cell_type":"markdown","metadata":{},"source":src.splitlines(True)}
def code(src):return {"cell_type":"code","execution_count":None,"metadata":{},
                      "outputs":[],"source":src.rstrip("\n").splitlines(True)}

cells = [md("""# Paper figures --- editable source

Every figure in `paper/ska_drift/main.tex`, with the **actual source** of its
plotting function inlined so colours, labels and layout can be tuned here.

**How this works.** Each cell below holds one `fig_*` function copied verbatim
from `ska_hi_plots.py` or `ska_hi_design.py`. Running a cell redefines that
function *in the notebook*, shadowing the module version, and the call at the
bottom of the cell writes the PNG to `figures/`. The paper picks it up with no
copy step --- `\\graphicspath` searches `../../figures/` first.

**Two sources of truth.** That is the cost of this notebook. Once a cell is
edited it diverges from the `.py`, and whichever ran last owns the PNG. When a
tuning is settled, paste it back into the script; until then, treat the
notebook as the live version and re-run it after any script change.

**VSCode.** Close this notebook before anything edits it on disk, or VSCode
will write its in-memory copy back over the change.

**Venv.** Run under `/home/geoff/gibbs_venv_312`, not `limTOD/.venv` ---
`ska_hi_plots` builds an HI cube and needs `pyccl`.
"""),
 md("## Setup"),
 code('''# Pick up edits to ska_hi_*.py without restarting the kernel. Without this,
# a name ADDED to a module after the kernel first imported it raises
# "ImportError: cannot import name ..." from the cached module object, and the
# only cure is a restart. Re-run this cell after editing a script.
%load_ext autoreload
%autoreload 2

import os, sys
sys.path.insert(0, "/home/geoff/limTOD/examples/SKA")
for _m in [m for m in list(sys.modules) if m.startswith("ska_hi")]:
    del sys.modules[_m]          # force a clean re-import on re-run

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import ska_hi_plots as P
import ska_hi_design as D
import ska_hi_kmodes as K
import ska_hi_analysis as A
from ska_common import gdsm_equatorial_sky_model
# brings in the palette, NMODES, RESDIR, collect, load_bands, load_realisations,
# freq_mode_density, BANDS, SEED_RUNS ... i.e. everything the pasted bodies use.
from ska_hi_design import *
# private helpers, which `import *` will not bring across
from ska_hi_plots import _tidy, _save, _title, _figsize, _load
# the style registries. ska_hi_design refers to these as P.<name>, but the
# ska_hi_plots functions pasted below use the bare names, so import them too.
# anything else the pasted bodies need by bare name is added automatically
# below by the generator -- do not maintain a list here, it drifts

%matplotlib inline
# Inline DISPLAY size only. Saved resolution is set by _save(), which uses
# STYLE["dpi"] (default 200) -- change that, not this, to alter the PNGs:
#     P.STYLE["dpi"] = 300
plt.rcParams["figure.dpi"] = 110
print("figures ->", P.FIGDIR, "| save dpi =", P.STYLE.get("dpi", 200))'''),
 md("""## Palette and series styles

`set_palette()` pushes colour changes into this notebook **and** into both
modules, so helpers like `_tidy` (grid colour) and `_title` follow.

**Identity is carried by the MARKER first and hue second**, so a series stays
readable in greyscale and for colour-vision-deficient readers. Two channels,
two registries:

- `P.STRAT_ALL[name]` -> `{color, marker, ls, label}` for every strategy
  (`drift`, `drift8`, `drift12`, `drift3t`, `raster`, `rasternarrow`). Edit a
  marker here and every figure that plots that strategy follows.
- `P.MONO` / `P.MONO_STRAT` -> the black-only styles used where a plot has just
  **two** series, so hue is dropped entirely and linestyle + marker carry it.

`P.STRAT` is deliberately only the two headline arms, because several `fig_*`
functions *iterate* it and index each strategy's own results keys --- iterating
the full registry would ask a run for data it does not contain.

Linestyle is the third channel and is spent on a second dimension (band, prior,
pipeline variant), not on identity --- see `hi_cleaning_depth`, where colour +
marker is the strategy and linestyle is the band."""),
 code("PALETTE = dict(\n"
      # emitted from the module's CURRENT values rather than hardcoded, so a
      # colour tuned in one place can never silently disagree with the other
      + "".join(f"    {n:<8s} = \"{getattr(P, n)}\",{c}\n" for n, c in (
          ("DRIFT",    "   # strategy 1"),
          ("RASTER",   "   # strategy 2"),
          ("THIRD",    "   # third arm (drift8 / drift12 / 1-f)"),
          ("SURFACE",  "   # marker fill / figure background"),
          ("INK",      "   # primary text"),
          ("INK2",     "   # panel titles"),
          ("MUTED",    "   # annotations"),
          ("GRID",     ""),
          ("BASELINE", ""),
      ))
      + ")\n" + '''
def set_palette(**kw):
    """Update the palette everywhere: notebook globals and both modules."""
    PALETTE.update(kw)
    for name, value in PALETTE.items():
        globals()[name] = value
        for mod in (P, D):
            if hasattr(mod, name):
                setattr(mod, name, value)
    return PALETTE

set_palette()'''),
 md("""## Data

Read once from `results/`. Nothing here recomputes --- these are all cached
`.npz` files, so the whole notebook runs in seconds."""),
 code('''exp, abl, rk = _load()          # experiment 006, its ablation, its rank study

lad      = np.load(os.path.join(RESDIR, "hi_ladder_f350_ns64.npz"))
lad_long = np.load(os.path.join(RESDIR, "hi_ladder_long_ns64.npz"))
pv       = np.load(os.path.join(RESDIR, "hi_priorvar_tk_f350_400_nc32_ns64.npz"))
matched  = np.load(os.path.join(RESDIR, "hi_ablation_matched_f350_400_nc32_ns64.npz"))
mech     = np.load(os.path.join(RESDIR, "hi_ablation_mech_f350_400_nc32_ns64.npz"))
modes    = np.load(os.path.join(RESDIR, "hi_modes_f350_400_nc32_ns64.npz"))

d = collect(matched, lad, modes, freq_mode_density(), extra_abl=mech)

# the published drift is 3 strips, drift12 is 12 -- the only two geometries
# with a measured HI floor, so the only overlay points fig_ladder can take
floors = {n: d[k]["floor"] for n, k in ((3, "drift"), (12, "drift12")) if k in d}
bands  = load_bands()
vals   = load_realisations()
kres   = K.compute()   # ~1 min: two eigendecompositions, the only slow cell

for s_ in ("drift", "drift8", "drift12", "raster"):
    if s_ in d:
        print(f"{s_:8s} {d[s_]['dens']:6.2f} modes/100deg^2   floor {d[s_]['floor']:7.1f}x")'''),
]

# Private helpers that `import *` will not carry. Resolved transitively and
# inlined into their figure's cell, because they ARE the drawing code worth
# tuning. Anything public is already in the namespace via `from ... import *`.
import ast, builtins
_PUBLIC = {n for n in dir(D) if not n.startswith("_")} | set(dir(builtins)) | {
    "_tidy", "_save", "_title", "_figsize", "_load",
    "np", "plt", "Line2D", "Rectangle", "os", "sys", "P", "D"}

def deps_of(mod, fn, seen=None):
    """Private module-level callables `fn` needs, depth first, no duplicates."""
    seen = [] if seen is None else seen
    tree = ast.parse(inspect.getsource(getattr(mod, fn)))
    names = sorted({n.id for n in ast.walk(tree)
                    if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)})
    for n in names:
        if n in _PUBLIC or n in seen or n == fn or not hasattr(mod, n):
            continue
        if not callable(getattr(mod, n)):
            continue
        deps_of(mod, n, seen)
        if n not in seen:
            seen.append(n)
    return seen

for name, mod, fn, call, sec, blurb in FIGS:
    helpers = deps_of(mod, fn)
    if helpers:
        print(f"  {fn}: inlining {helpers}")
    src = "\n\n".join([inspect.getsource(getattr(mod, h)) for h in helpers]
                       + [inspect.getsource(getattr(mod, fn))])
    cells.append(md(f"## `{name}`\n\n*{sec}.* {blurb}\n\n"
                    f"Source: `{mod.__name__}.{fn}`. "
                    f"Every `fig_*` takes `text=` to override the title/subtitle "
                    f"without touching the body, e.g.\n"
                    f"`{call[:-1]}, text=('My title', 'My subtitle'))`."))
    cells.append(code(src + "\n\n" + call))

# --- resolve the pasted bodies' module-level names automatically -----------
# The figure sources are copied verbatim, so any module CONSTANT they use
# (TF_PANELS, MONO_STRAT, FLOOR_AMP ...) has to exist in the notebook too.
# Maintaining that list by hand has broken the notebook twice; compute it.
def _free_names(src):
    # cells may carry IPython magics, which are not valid Python
    src = "\n".join(l for l in src.split("\n")
                    if not l.lstrip().startswith(("%", "!")))
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return set()
    loaded = {n.id for n in ast.walk(tree)
              if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    bound = {n.id for n in ast.walk(tree)
             if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}
    bound |= {f.name for f in ast.walk(tree) if isinstance(f, ast.FunctionDef)}
    for f in ast.walk(tree):
        if isinstance(f, (ast.FunctionDef, ast.Lambda)):
            a = f.args
            bound |= {x.arg for x in list(a.args) + list(a.kwonlyargs)
                      + list(a.posonlyargs)}
            for x in (a.vararg, a.kwarg):
                if x:
                    bound.add(x.arg)
    return loaded - bound

_provided = ({n for n in dir(D) if not n.startswith("_")}   # from ... import *
             | set(dir(builtins))
             | {"np", "plt", "Line2D", "Rectangle", "os", "sys", "P", "D", "K",
                "hp", "A", "gdsm_equatorial_sky_model",
                "_tidy", "_save", "_title", "_figsize", "_load"})
_need, _seen = {}, set()
for c in cells:
    if c["cell_type"] != "code":
        continue
    for n in sorted(_free_names("".join(c["source"]))):
        if n in _provided or n in _seen:
            continue
        _seen.add(n)
        for mod in (P, K):
            if hasattr(mod, n):
                _need.setdefault(mod.__name__, []).append(n)
                break
if _need:
    extra = "\n".join(f"from {m} import {', '.join(sorted(v))}"
                       for m, v in sorted(_need.items()))
    print("  auto-imported into the notebook:",
          {m: sorted(v) for m, v in _need.items()})
    setup = next(c for c in cells if c["cell_type"] == "code"
                 and "import ska_hi_plots as P" in "".join(c["source"]))
    body = "".join(setup["source"])
    marker = "%matplotlib inline"
    assert marker in body
    body = body.replace(marker,
                        "# resolved automatically by the generator from the "
                        "pasted figure bodies\n" + extra + "\n\n" + marker)
    setup["source"] = body.splitlines(True)

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 5}
out = "/home/geoff/limTOD/examples/SKA/ska_hi_paper_figures.ipynb"
json.dump(nb, open(out, "w"), indent=1)
print("wrote", out, len(cells), "cells")
