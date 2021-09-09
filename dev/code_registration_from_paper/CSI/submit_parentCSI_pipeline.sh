#!/bin/bash

#SBATCH -p jic-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH -t 2-00:00  # max allowed time before kicked off (D-H:m)
#SBATCH -o ./slurm_output/slurm.%j.out # STDOUT write path
#SBATCH -e ./slurm_output/slurm.%j.err # STDERR write path
#SBATCH --mail-type=END,FAIL # notifications by email in the event of these events
#SBATCH --mail-user=email@address # emails sent to address
#SBATCH --mem 16000 # memory pool
#SBATCH --constraint=intel # require intel processor as TensorFlow library was compiled to use SSE4.1 instructions
# end of slurm job manager parameters

set -eux  # report as failiure + quit if any step fails

export TF_CPP_MIN_LOG_LEVEL=2  # suppress tensorflow binary warning

# get parameters from "do_run_parent_master_arabidopsis.sh"
all_args=("$@")
RDS=${all_args[0]} # gene expression data rds file name
maxNumParents=${all_args[1]} # max number of simultaneous regulators of each gene
target=${all_args[2]}  # id of target, regulated gene
amalgamate=${all_args[3]}  # whether to sum brassica paralogues, or treat individually
symbols=${all_args[@]:4} # all other all_args from (4 to end) are genes to include in analysis

# specify input & output directories:
# directory where gene expression data tables are saved as .rds files
RDS_DIR='./data/'
# specify directory where intermediate pre-processed gene expression data is saved
DATA_DIR='./intermediate_data/'$RDS'/T='$target'/'


# preprocess expression data:
# - filter to only include genes of interest.
# - apply standard scaling
# - reformat from long to wide format table
Rscript ./parent_prepare_data.R $RDS $RDS_DIR $DATA_DIR $amalgamate "${symbols[@]}"

CURR_DATA_DIR=$DATA_DIR

# make candidate parent combinations
COMBOS=$(python ./parent_make_parent_combos.py "$CURR_DATA_DIR" "$maxNumParents" "$target")
echo COMBOS: "$COMBOS"

for C in ${COMBOS}
do
    echo 'calculating parent likelihood!'
    # CALC. BEST Log Marginal Lieklihood (LML) FOR TARGET GENE ASSUMING IS REGULATED BY PARENT COMBO C
    python ./parent_CSI.py "$CURR_DATA_DIR" "$target" "$C" "$RDS"
done # end iterating over C
