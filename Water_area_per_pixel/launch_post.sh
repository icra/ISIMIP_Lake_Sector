#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1
###new version of R

spack load r@3.6.3%gcc@9.4.0
spack load r-raster
spack load r-sf
spack load r-rgdal@1.5-19%gcc@9.4.0
spack load r-lwgeom


Rscript Calculate.R



