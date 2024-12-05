# evolve_binary
Evolve a binary through merger using PN equations and a NR remnant surrogate; see example.py for an example!

Tips for Installation:
- instal Julia via `curl -fsSL https://install.julialang.org | sh`
- create a mcfacts-dev conda environment via `make setup` in your McFACTS repository
  - but change the Makefile so you install python=3.12.0
- activate the mcfacts-dev conda environment
- install `sxs` via `conda install -c conda-forge sxs numba::numba numba::llvmlite`
- install the other dependencies in `requirements.txt` with `conda`
- download `surrogate.joblib` from [this Google Drive](https://www.dropbox.com/scl/fo/p33rqfjew5vu5qzksu32w/AEr4moWujITfl46ezybjE1Q?rlkey=1lladw82d8twlpt2xi5hidscv&st=xctpnkyj&dl=0)
- try running `python example.py` (ane make sure `surrogate.joblib` is in this directory)
