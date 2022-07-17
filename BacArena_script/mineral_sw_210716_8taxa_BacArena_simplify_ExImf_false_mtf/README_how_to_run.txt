###################### How to make your own R script to run BacArena?
# Step 1.
# Open the following R script:
Rscript_BacArena__Best_models_STY_8taxa_20220615__mtf_rmRate_tp_ExInf_replenish_growlimit_parL_diffModel_addSubsF_auto1.R
# Modify line 71-103, 109, and 127-165.

# Step 2.
# Open the following R script:
Rscript_BacArena__Best_models_STY_8taxa_20220615__mtf_rmRate_tp_ExInf_replenish_growlimit_parL_diffModel_addSubsF_auto2.R
# Modify line 100-115 according to the previous step.

# Step 3. Make your own diet.

# Step 4. run Cmds_BioSAK_prepare_qsub_job_files_auto_example.txt
# To prepare files:
cmds_200iter_J16_1v1.txt
# and
js_header_100G_2h.sh

# Step 5. Get qsub job script in batch and run them on katana:
