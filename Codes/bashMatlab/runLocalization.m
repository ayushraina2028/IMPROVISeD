function runLocalization(protein_name)


    cd('../')
    addpath(genpath(pwd))
    cd('bashMatlab')

    % Define the paths using the protein name argument
    chain_A_file = ['../crosslinks/' protein_name 'CLs/eq_dists_chain_A.csv'];
    chain_B_file = ['../crosslinks/' protein_name 'CLs/eq_dists_chain_B.csv'];
    crosslink_file = ['../crosslinks/' protein_name 'CLs/chain_A_crosslink30_chain_B_LYS_Ca.csv'];

    % Call the MATLAB function with the necessary arguments
    [Y, info] = callLowrankBamdev(chain_A_file, chain_B_file, crosslink_file, 30);

    % Optionally, save results or print information
    disp('MATLAB tasks completed.');
end

