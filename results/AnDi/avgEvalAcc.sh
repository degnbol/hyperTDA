#!/usr/bin/env zsh
grep -i 'avg eval acc' */CNN/*.out | sed 's,_.*slurm,\t,' | sed 's,\.out:avg eval Acc = ,\t,' | grep -v 123
