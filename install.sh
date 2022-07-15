#!/usr/bin/env zsh
# make a command "git root" that gives the root folder of the repo.
git config alias.root 'rev-parse --show-toplevel'
# make sure submodules of other publications are initialized.
git submodule update --init
# install julia packages
./install.jl
