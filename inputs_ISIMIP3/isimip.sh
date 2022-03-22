#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=1

spack load r@3.6.3%gcc@9.4.0 arch=linux-centos7-skylake
spack load r-raster@3.4-5%gcc@9.4.0 arch=linux-centos7-skylake
spack load r-rgdal@1.5-19%gcc@9.4.0 arch=linux-centos7-skylake

Rscript Calculate_rasters.R
