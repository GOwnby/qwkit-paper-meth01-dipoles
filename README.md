[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23444554.svg)](https://doi.org/10.5281/zenodo.23444554)

# Dipoles of nonparabolic quantum-well states: why orthonormalizing the bound states is not enough

**Author:** Gracen Ownby (University of Oklahoma)  
**Date:** October 2026  
**Manuscript PDF:** [`METH01_paper.pdf`](./METH01_paper.pdf)  
**LaTeX Source:** [`paper_meth01/`](./paper_meth01/)  
**Study Code & JSON Proofs:** [`studies/meth01/`](./studies/meth01/)  

## Overview
This repository is the permanent, self-contained mathematical and computational archive for **METH01**: *Dipoles of nonparabolic quantum-well states: why orthonormalizing the bound states is not enough*. It includes:
1. The complete reviewed paper (`METH01_paper.pdf`) and its LaTeX source (`paper_meth01/`).
2. The complete `qwkit` (v0.7.0) Python package (`qwkit/`) and its 69 automated validation tests (`tests/`), which verify analytic limits, BenDaniel-Duke finite-well convergence, exact mirror-symmetry oddness $\chi^{(2)}(-F) = -\chi^{(2)}(F)$ to $10^{-11}$, two-band Kane orthogonality ($10^{-13}$), energy-conserving RK4 coupled-mode propagation, implicit-differentiation tangents, and RPA sum rules.
3. The study scripts and precomputed JSON datasets (`studies/meth01/`) that generate every table, bound, and figure in the manuscript.

## Quick Verification (Mathematical & Numerical Proof)

Install `qwkit` in editable mode (requires Python >= 3.10, `numpy`, `scipy`, `matplotlib`):
```bash
pip install -e .
```

Run the single-command verification suite (runs all 69 unit tests, verifies the study numbers, and regenerates the manuscript figures):
```bash
python prove_paper.py
```

## License
Code is released under the MIT License (`LICENSE`). The manuscript and figures are (c) 2026 Gracen Ownby.
