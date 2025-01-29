#!/bin/bash

#SBATCH --job-name="SP_v1"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6                # Allocate 6 CPUs to match the number of scripts
#SBATCH --partition=cpuq                  # Use the appropriate partition
#SBATCH --mem=150GB                        # Adjust memory as needed
#SBATCH --time=5-23:45:13                  # Maximum time allocation
#SBATCH --error=test_%j.sterr              # Standard error log
#SBATCH --output=test_%j.log               # Standard output log

# Load required modules
module load R singularity gcc/10.2.0

# Define base directory and Singularity image
BASE_DIR="SP_models"
SIF_IMAGE="$BASE_DIR/synthos.sif"
REPORT_SCRIPT="SP_models/scripts/report_v2.qmd"

# Array of scripts to run
SP_SCRIPTS=(
    "SP_models/pipeline_SP_v1.R"
    "SP_models/pipeline_SP_v2.R"
    "SP_models/pipeline_SP_v3.R"
    "SP_models/pipeline_SP_v4.R"
    "SP_models/pipeline_SP_v5.R"
    "SP_models/pipeline_SP_v6.R"
)

# Run SP scripts in parallel using SLURM job arrays
for i in $(seq 0 $((${#SP_SCRIPTS[@]} - 1))); do
    srun -n 1 singularity exec -B "$BASE_DIR:/home/Project" "$SIF_IMAGE" Rscript "${SP_SCRIPTS[$i]}" &
done

# Wait for all background tasks to finish
wait

# Run Quarto report generation **after** all scripts finish
singularity exec -B "$BASE_DIR:/home/Project" "$SIF_IMAGE" quarto render "$REPORT_SCRIPT"

echo "All SP scripts and report generation completed."