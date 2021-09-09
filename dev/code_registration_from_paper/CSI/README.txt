CSI supporting information README:
The directory "CSI/" contains scripts and gene expression data used as input to the Causal Structure Inference algorithm. 

NB this is provided as a means to document methods, and is not intended to be used as software!


R_env_details.txt lists the R library versions used. py36_env_details.txt details the conda python environment used.

- run_submit_parentCSI_pipeline.sh: sets up the environment, defines genes of interest,  and submits submit_parent_CSI_pipeline.sh to a slurm job manager.

- submit_parent_CSI_pipeline.sh: a slurm job submission script which calls scripts to pre-process the expression data, define the considered candidate regulators of each gene, and use CSI to calculate the likelihood of regulation by each candidate regulatory set. It calls "parent_prepare_data.R", "parent_make_parent_combos.py" and "parent_CSI.py".

- parent_prepare_data.R : preprocesses and normalises expression data

- parent_make_parent_combos.py: defines the possible parent regulatory combinations of each gene (for which likelihood will be calculated)

- parent_CSI.py: apply the CSI algorithm to calculate the likelihoods of each regulatory set.
