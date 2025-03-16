
# Set the working directory.
cd /omics/groups/OE0436/internal/Linh/Code

# Load R so that the R script can be executed.
module load R/4.3.0
module load gdal/3.0.2
module load hdf5/1.8.18
module load gcc/11.1.0
module load binutils/2.34
module load libpng/1.6.37
module load jags/4.3.0
module load freetype/2.10.0
module load imagemagick/6.9.12



# Run the R script.
Rscript ST_LIANA_TLS_combine.R
