% automate_script.m

% Read the protein name from the file
filename = '../bashScripts/protein_name.txt'; % Update the path as necessary
if isfile(filename)
    fileID = fopen(filename, 'r');
    protein_name = fscanf(fileID, '%s');
    fclose(fileID);
else
    error('Protein name file not found.');
end

% Change directory to the parent directory and add paths
cd('../Registration');
addpath(genpath(pwd));
cd('../bashMatlab');

% Run the file and save results
x_n_indx = readFileMakeXnIndx_v2('run1_Y30_chain_A_xyz.txt,run1_Y30_chain_B_xyz.txt', ...
                                  'run1_Y30_chain_A_indx.txt,run1_Y30_chain_B_indx.txt', ...
                                  'run1_Y30crosslink_xyz.txt', 'run1_Y30crosslink_indx.txt');
save('x_n_index_run1_Y30.mat', 'x_n_indx');

% Change directory and add paths
cd('/data2/nmr/our_algo_final_test/code/lib/Factorize/Factorize/');
addpath(genpath(pwd));
cd('/data2/nmr/Ayush/registration_test');
addpath(genpath(pwd));

% Perform registration and save results
[X_noref, atom_map, rot, L_inv, B] = doRegForGroup([1,2,3], x_n_indx);

% Change directory and write results to file
output_dir = ['/data2/nmr/complex_modelling/final1/Codes/Registration/', protein_name, 'Data/'];
cd(output_dir);
writematrix(X_noref', 'X_noref_run1Y30.csv');

