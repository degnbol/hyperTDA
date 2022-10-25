# hyperTDA
Piece-wise Linear curves converted to point-clouds, analysed with Persistent Homology, represented as HyperGraphs.
Associated with preprint "Hypergraphs for multiscale cycles in structured data" at https://doi.org/10.48550/arXiv.2210.07545

## REQUIREMENTS
- Julia (tested on v1.7.3)
- Python3 (tested on v3.9.7)
- unix shell (tested on zsh v5.8.1)

## INSTALL
- `./install.sh` which sets a git alias `git root`, initializes submodules and 
  runs `./install.jl`.
- Julia packages are listed in `./install.jl` and tested versions in `Manifest.toml` and `Project.toml`. Running `./install.sh` installs the tested versions to a local environment and `./install.jl --global` installs the newest available versions to the global environment.
- Python packages are listed in `requirements.txt` and can be installed with 
  e.g. `pip install -r requirements.txt` or `conda create --name=hyperTDA 
  --file=requirements.txt` requiring either `pip` (pip3) or `conda` (e.g. miniconda3 or miniforge).

## EXAMPLE USE
Example input files, commands and output files are provided in subfolders of `examples/`. See `RUNME.sh` for each example.

## PAPER RESULTS
Results discussed in the paper are provided in `results/`, see e.g. jupyter notebooks in subfolders.

## TIMING
- Git clone ~ 3 min
- Julia package install ~ 3 min
- Python package install ~ 1 min
- running examples
  - default pipeline (two curves of 132 and 101 points) ~ 30 sec
  - with interpolation (two curves each of 500 points) ~ 1 min 30 sec
