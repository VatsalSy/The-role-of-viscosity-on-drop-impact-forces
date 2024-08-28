#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=192
#SBATCH --time=100:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=v.sanjay@utwente.nl

source ~/.bash_shell

MAXlevel="13"
tmax="10"
We="20"
Ohd="0.05"
Ohs="1e-4"
Bo="0.25"
Ldomain="6"
a0="1e0"

srun --mpi=pmi2 -n 192 --mem-per-cpu=1750mb bounce_noSphere $MAXlevel $tmax $We $Ohd $Ohs $Bo $Ldomain $a0