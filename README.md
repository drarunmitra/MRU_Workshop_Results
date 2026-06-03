# MRU Workshop — Results Section (Session 2)

Quarto Reveal.js deck for *Writing the Results Section: Text, Tables & Figures*
(Research Methodology Workshop, Day 2).

## Files

| File | Purpose |
|---|---|
| `index.qmd` | The slide deck (source) |
| `style.css` | Tiny house stylesheet (logo + maroon `.important` + `.inverse`) |
| `aiims_bibinagar_logo.png` | Deck logo |
| `generate_figures.R` | Regenerates the K–M, ROC and forest figures into `images/` |
| `images/` | Pre-rendered figures used by the deck |
| `_extensions/` | `tabset` Quarto plugin (needed to render) |
| `Results_Section_Workshop_Session2.md` | Full content outline the deck is built from |
| `index.html` | Last rendered deck (self-contained) |

## Render (Mac or Windows)

```bash
quarto render index.qmd      # produces a self-contained index.html
```

The deck itself needs no R (mermaid + pre-made PNGs only).

## Regenerate the figures (optional)

Needs R with `ggplot2`, `pROC`, `survival`, `here`:

```bash
Rscript generate_figures.R   # writes km_curve.png, roc_curve.png, forest_plot.png
```

Data are simulated and **seeded**, so the figures reproduce exactly.
