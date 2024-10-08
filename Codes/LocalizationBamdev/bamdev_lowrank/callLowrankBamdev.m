
%%
%  callLowrankBamdev('eq1.txt', 'eq2.txt', '1dfj_DSSO_9.csv', 30)
%%
function [Y, infos] = callLowrankBamdev(eq1file, eq2file, crosslinksfile, crosslinks_val, Y0, noise_sigma, embed_dim)
%%
% I/p:  eq1file:   see eq1.txt
%       eq2file:    
%       crosslinksfile
%       crosslinks_val
%       Y0
%       noise_sigma
%       embed_dim
%      crosslinks:  no. of links x 2 ... (row, col)
%%
    if nargin <5
        Y0 = [];
    end
    if nargin <6
        noise_sigma = 0.1;
    end
    if nargin <7
        embed_dim = 3;
    end

	%%
	eq1_arr     = table2array(makeArrayFromEq(eq1file));
	eq2_arr     = table2array(makeArrayFromEq(eq2file));
	crosslinks_ = table2array(makeArrayFromCrosslinks(crosslinksfile));
    crosslinks  = [double(crosslinks_(:,1)), double(crosslinks_(:,3))];
	[distmat, dict_indx1, dict_indx2]   = fillDistMat_Eq(eq1_arr, eq2_arr, crosslinks, crosslinks_val);
    
 	%% add noise to crosslinks %%
	std_crosslinks = crosslinks_val;
	for i = 1:size(crosslinks)
		distmat(dict_indx1(crosslinks(i,1)), dict_indx2(crosslinks(i,2))) = distmat(dict_indx1(crosslinks(i,1)), dict_indx2(crosslinks(i,2))) + noise_sigma * std_crosslinks;
	end	    
    
        %%
	[I, J, ~] = find(~isnan(distmat));
    dist_ = NaN(length(I),1);
    for i = 1:length(I)
       dist_(i) = distmat(I(i), J(i)); 
    end
        

	%% initialize Y0 (if not supplied) %%
    % Starting approximation rank
    p = 1;
    if isempty(Y0)
        Y0 = randn(length(distmat), p) ;
    end       
	
	%% 

    % Set randstream based on clock
    s = RandStream.create('mt19937ar','seed',sum(100*clock));
    %RandStream.setDefaultStream(s);
    RandStream.setGlobalStream(s);
    % Select optimization algorithm
    %methodName = 'Gradient Descent';
    methodName = 'Gradient Descent';

    % Parameters
    params.pmax = embed_dim;
    params.tol = 1e-3;
    params.vtol = 1e-3;
    params.verb = true; % "Bark" only when asked
    
    switch methodName,
        case 'Gradient Descent',
            methodFun = @gd_dist_completion;
        case 'Trust Region',
            methodFun = @tr_dist_completion;
        otherwise,
            error('Unknown method %s',methodName);
    end

    
	[Y, infos] = lowrank_dist_completion(methodFun,I,J,dist_,Y0,params);
    opfile = sprintf('Y_%d.csv',crosslinks_val)
    writematrix(Y,opfile)
end

function eq1 = makeArrayFromEq(filename, dataLines)
%IMPORTFILE Import data from a text file
%  EQ1 = IMPORTFILE(FILENAME) reads data from text file FILENAME for the
%  default selection.  Returns the data as a table.
%
%  EQ1 = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  eq1 = importfile("/data2/nmr/complex_modelling/localization/bamdev_low_rank_matrix_compeltion/distcompletion_25jan_2013/distcompletion_25jan_2013/eq1.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 19-Aug-2024 13:09:32

	%% Input handling

	% If dataLines is not specified, define defaults
	if nargin < 2
	    dataLines = [2, Inf];
	end

	%% Setup the Import Options
	opts = delimitedTextImportOptions("NumVariables", 3);

	% Specify range and delimiter
	opts.DataLines = dataLines;
	opts.Delimiter = ",";

	% Specify column names and types
	opts.VariableNames = ["resii", "resij", "dist"];
	opts.VariableTypes = ["double", "double", "double"];
	opts.ExtraColumnsRule = "ignore";
	opts.EmptyLineRule = "read";

	% Import the data
	eq1 = readtable(filename, opts);
    
end

function dfjDSSO9 = makeArrayFromCrosslinks(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DFJDSSO9 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  DFJDSSO9 = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  dfjDSSO9 = importfile("/data2/nmr/complex_modelling/localization/bamdev_low_rank_matrix_compeltion/distcompletion_25jan_2013/distcompletion_25jan_2013/1dfj_DSSO_9.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 20-Aug-2024 12:00:04

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["res1", "prot1", "res2", "prot2"];
opts.VariableTypes = ["double", "string", "double", "string"];
opts = setvaropts(opts, [2, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dfjDSSO9 = readtable(filename, opts);

end

function [dict_indx1, dict_indx2, tmp_indx1, tmp_indx2] = getTmpIndx(eq1_arr, eq2_arr, crosslinks_arr)
%%
%
%%
	tmp_indx1 = unique(union(eq1_arr(:,1), union(eq1_arr(:,2), crosslinks_arr(:,1))))
	tmp_indx2 = unique(union(eq2_arr(:,1), union(eq2_arr(:,2), crosslinks_arr(:,2))))

    % -- dictionary not availabe in 2019b -- %
	%dict_indx1 = dictionary(tmp_indx1, 1:size(tmp_indx1));
	%dict_indx2 = dictionary(tmp_indx2, 1:size(tmp_indx1));
    
    dict_indx1 = containers.Map(tmp_indx1, 1:length(tmp_indx1));
    
    dict_indx2 = containers.Map(tmp_indx2, length(tmp_indx1)+[1:length(tmp_indx2)]);
    
    writematrix(tmp_indx1,'tmp_indx1.txt')
    writematrix(tmp_indx2,'tmp_indx2.txt')

end

function [distmat, dict_indx1, dict_indx2] = fillDistMat_Eq(eq1_arr, eq2_arr, crosslinks_arr, crosslinks_val)                     
%%
%
%%

	[dict_indx1, dict_indx2, tmp_indx1, tmp_indx2] = getTmpIndx(eq1_arr, eq2_arr, crosslinks_arr);
	
	size_ = length(union(tmp_indx1, tmp_indx2));
    distmat = NaN(size_,size_);

        %% ------- fill the distance matrix ---- %%
		
		for i = 1:size(eq1_arr)
			distmat(dict_indx1(eq1_arr(i,1)), dict_indx1(eq1_arr(i,2))) = eq1_arr(i,3);
			%distmat[ dict_indx1{eq1_arr[i,2]}, dict_indx1{eq1_arr[i,1]}] = eq1_arr[i,3]
		end
		for i = 1:size(eq2_arr)
        	        distmat(dict_indx2(eq2_arr(i,1)), dict_indx2(eq2_arr(i,2))) = eq2_arr(i,3);
                	%distmat[ dict_indx2{eq2_arr[i,2]}, dict_indx2{eq1_arr[i,1]}] = eq2_arr[i,3]
	    end

		for i = 1:size(crosslinks_arr)
			distmat(dict_indx1(crosslinks_arr(i,1)), dict_indx2(crosslinks_arr(i,2))) = crosslinks_val;
			%dist[dict_indx1{crosslinks_arr[i,2]}, dict_indx2{crosslinks_arr[i,1]}] = crosslinks_arr[i,3]
		end
	%%

end


	
