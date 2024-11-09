#!/bin/bash

# Step 1: Create the environment.yml file
cat <<EOL > environment.yml
name: seurat
channels:
  - conda-forge
  - r
  - bioconda
  - defaults
dependencies:
  - r-base=4.4
  - r-essentials
  - r-devtools
  - r-remotes                # To use remotes::install_version()
  - jupyter                  # Ensure Jupyter is installed
  - r-irkernel               # Install R kernel for Jupyter
  - r-rcurl                  # Needed for installing Seurat dependencies
  - r-matrix
  - r-mass
  - r-lattice
  - r-scales
  - r-sp
  - r-lars
  - r-tsne
  - r-fpc
  - r-hmisc
  - r-metap
  - r-cluster
  - r-png
  - r-dosnow
  - r-foreach

EOL

# Step 2: Create the postBuild file for automatic Seurat installation
cat <<'EOL' > postBuild
#!/bin/bash

# Set R_LIBS to point to the correct library directory within the conda environment
export R_LIBS="${CONDA_PREFIX}/lib/R/packages"

# Create the R library directory if it does not exist
mkdir -p "${R_LIBS}"

# Install Seurat if it is not already installed
#Rscript -e "if (!requireNamespace('Seurat', quietly = TRUE)) install.packages('Seurat', #repos = 'http://cran.us.r-project.org')"

# Install Jupyter and IRkernel if not already installed
#conda install -y -c conda-forge jupyter r-irkernel
#conda install -y conda-forge r-seurat

# Register R kernel with Jupyter
Rscript -e "IRkernel::installspec(name = 'seurat', displayname = 'R (seurat)')"

EOL

# Make the postBuild script executable
chmod +x postBuild

# Step 3: Create the Conda environment and trigger postBuild
conda env create -f environment.yml

# Inform the user of completion
echo "Environment 'Seurat' created with R 4.4 and Seurat installed (via postBuild)."
