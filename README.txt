Examples folder contains 3 subfolder for 1dfj, lcn2 and mmp9

1. To start, open Codes/bashScripts/makeCrosslinks.sh and set the proteinName="1dfj" or "lcn2" or "mmp9", it is already set to "1dfj" for now. This file will search for pdb files present in Examples/proteinName/ 
For now I have kept chain_A.pdb and chain_B.pdb as the name for 2 pdb files.

1.2) Running the above script will generate a file named "chain_A_crosslink30_chain_B_LYS_Ca.csv" which are generated crosslinks, to add experimental crosslinks, we can skip the first step and add our crosslinks to this file directly. This csv file is available in Codes/crosslinks/proteinNameCLs/ 
for example Codes/crosslinks/1dfjCLs/

2. Now to extract the equality constraints, go to Codes/bashScripts/makeCrosslinks2.sh, again we have to add the proteinName in this script, and this will use previously generated file i.e "chain_A_crosslink30_chain_B_LYS_Ca.csv" to generate the equality constraints and will save it into eq_dists_prots_model1.csv in the same location.

3. Moving to 3rd script which is separateEqCrossLinks.sh, we have to provide proteinName and path of previously generated files, in this script and it will separate the extracted equality bounds into 2 different files named "eq_dists_chain_A.csv" and "eq_dists_chain_B.csv"

(Please run tree final1/ to see location of different files)

4. There is another script named remove13col.sh, when we run this, it will remove the residue names and updates the above 2 generated files. Now those 2 files contains 
i,j,dist(i,j)

5. Now we have to solve the localization problem, for that give the generated equality constraints file generated in step 3 and crosslink file from step 1 in Codes/bashMatlab/runLocalization.m
To run this go to Codes/bashScripts/localization.sh and execute it, this will solve the localization problem and will save the generated files in 
Codes/LocalizationBamdev/1dfjResult/

6. Next step comes which is generating files, which will be used in making x_n_index for registration part. To do this, run the 6th bash script named localization.sh, it will generate 6 files, 3 coordinate files for chain_A, chain_B, crosslinks and 3 index files for the same. These files will be saved in Codes/Registration/proteinNameData/

7. Registration Step - run the runRegiste.sh script to do the registration, final files will be saved in Registration/proteinName/


Visualizing 1dfj structure
Go to Codes/Registration/1dfjData/
1. pymol ../../../Examples/1dfj/chain_A.pdb ../../../Examples/1dfj/chain_B.pdb
2. PYMOL> run ../showCoords.py
3. PYMOL> showRegPDBs('chain_A,chain_B','X_noref_run1Y30.csv', 'run1_Y30_chain_A_indx.txt,run1_Y30_chain_B_indx.txt')

Visualizing lcn2 structure
Go to Codes/Registration/lcn2Data/
1. pymol ../../../Examples/lcn2/chain_A.pdb ../../../Examples/lcn2/chain_B.pdb
2. PYMOL> run ../showCoords.py
3. PYMOL> showRegPDBs('chain_A,chain_B','X_noref_run1Y30.csv', 'run1_Y30_chain_A_indx.txt,run1_Y30_chain_B_indx.txt')

Visualizing mmp9 structure
Go to Codes/Registration/mmp9Data/
1. pymol ../../../Examples/mmp9/chain_A.pdb ../../../Examples/mmp9/chain_B.pdb
2. PYMOL> run ../showCoords.py
3. PYMOL> showRegPDBs('chain_A,chain_B','X_noref_run1Y30.csv', 'run1_Y30_chain_A_indx.txt,run1_Y30_chain_B_indx.txt')

For fast testing purposes, do not change anything, just put 2 pdb files in Examples/1dfj/ , name of the pdb files should be chain_A.pdb and chain_B.pdb, and run all scripts in order. If we are not generating synthetic crosslinks, leave first script, do step 1.2 and run remaining scripts. In this way, time taken will be less per iteration.
I will also combine all these scripts into one in a while, then we will implement logic which will run all this with different subsets of
 "chain_A_crosslink30_chain_B_LYS_Ca.csv"
