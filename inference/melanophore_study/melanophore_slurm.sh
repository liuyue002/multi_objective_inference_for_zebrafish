!/bin/bash

#SBATCH --account=username
#SBATCH --qos=standby
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --job-name=alpha_beta_param_scan
#SBATCH --output=/scratch/username/
#SBATCH --error=/scratch/username/
#SBATCH --mail-user=your.email@example.com
#SBATCH --mail-type=NONE

# Loads Matlab and sets the application up
module load matlab/R2024a
module load conda/2025.02
# make sure relevant python packages are installed in my_env
conda activate my_env
# If Python complains about libstdc++.so.6 missing, try to manually locate it, and do
export LD_PRELOAD=/path/to/libstdc++.so.6

cd /home/liu4194/zebrafish_inference/zebrafish_abm_tda

matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "melanophore_sim(${1},${2},${3},${4},${5}); exit"
