# hyperTDA
Piece-wise Linear curves converted to point-clouds, analysed with Persistent Homology, represented as HyperGraphs.

## INSTALL
- Julia. CNN was tested under Julia v1.8.3 (CUDA and Flux had issues on linux for v1.8.5, v1.9.3, and v1.9.4).
- `./install.sh` which sets a git alias `git root`, initializes submodules and 
  runs `./install.jl`.
- Julia packages are installed to local environment by `./install.sh` but they 
  can also be installed globally with `./install.jl --global`
- Python packages are listed in `requirements.txt` and can be installed with 
  e.g. `pip install -r requirements.txt` or `conda create --name=hyperTDA 
  --file=requirements.txt` requiring either `pip` or `conda`.

