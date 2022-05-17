#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=512M
#SBATCH --nodelist=compute-12

spack load r@4
spack load r-raster@3.4-5%gcc@9.4.0 arch=linux-centos7-x86_64
spack load r-sf@0.9-7%gcc@9.4.0 arch=linux-centos7-x86_64
spack load r-rgdal@1.5-19%gcc@9.4.0 arch=linux-centos7-x86_64
spack load r-lwgeom@0.2-5%gcc@9.4.0 arch=linux-centos7-x86_64

Rscript area_pixel_year_rest.R
Rscript area_pixel_year_final.R

wait

