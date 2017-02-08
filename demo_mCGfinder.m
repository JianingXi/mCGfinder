% If you want to analyze the demo data (TCGA BRCA somatic mutation data) through mCGfinder on with
% the default configures, please run `./demo_mCGfinder.m` and then the result file will be saved
% in the directory `./output/` as `./output/Result.mat`.
% 
% If you want to analyze a user-specific data, the sample-gene mutation binary matrix, the sample names 
% and the gene symbols of the samples must be provided. These three variable are then saved to be MATLAB 
% variables 'mutation_data' (binary matrix), 'sample_names' and 'gene_names' in the `./data` directoryas
% as the demo data. 
% Then update the variable 'dir_data' in file `./demo_mCGfinder.m` with the new data.
% 
% The configurations of mCGfindercan be changed in script file `./demo_mCGfinder.m`, and the descriptions
% of these parameters are provided below:
% 
%         =================================================================================================
%         | PARAMETER NAME       | DESCRIPTION                                                            |
%         =================================================================================================
%         |CompLeastProportion   |Least sample proportion included in each components, which represents   |
%         |                      |minimum proportion of the samples in every components given by the      |
%         |                      |mCGfinder. The default proportion is set to 15%.                        |
%         -------------------------------------------------------------------------------------------------
%         |maxCompoent           |Maximum number of components, which denotes the number of components    |
%         |                      |given components given by mCGfinder at most. The default number is 5.   |
%         -------------------------------------------------------------------------------------------------
%         |NetConf.lambda_T      |The tuning parameter of network regularization, which is used to balance|
%         |                      |the fitness of the model (first term) and the smoothness of the scores  |
%         |                      |of connected genes (second term). The default number is 0.1.            |
%         -------------------------------------------------------------------------------------------------
% 
% 
% The descriptions of output variables of mCGfinder are provided below:
% 
%         =================================================================================================
%         | VARIABLE NAME        | DESCRIPTION                                                            |
%         =================================================================================================
%         |detected_genes        |Genes detected by mCGfinder as significantly mutated cancer genes.      |
%         -------------------------------------------------------------------------------------------------
%         |S_sample_indicator    |The sample indicator vectors of all component, which indicates the      |
%         |                      |assignment of tumour samples to the every components. The i-th          |
%         |                      |coefficient being 1 represents that the i-th samples are included in the|
%         |                      |component, and 0 otherwise.                                             |
%         -------------------------------------------------------------------------------------------------
%         |Symbol_Net            |The investigated gene list in the gene interaction network.             |
%         -------------------------------------------------------------------------------------------------
%         |G_gene_score          |The gene score vectors of all components, of which the coefficients are |
%         |                      |related to the gene lists variable 'Symbol_Net', and a higher value of  |
%         |                      |a certain coefficient presents a larger potential of the gene to be     |
%         |                      |cancer gene candidate.                                                  |
%         -------------------------------------------------------------------------------------------------
%         |Q_values              |The q-values of all investigated genes in variable 'Symbol_Net', which  |
%         |                      |are obtained by Benjamini-Hochberg false discovery rates control of the |
%         |                      |p-values of the investigated genes.                                     |
%         -------------------------------------------------------------------------------------------------
% 

bin_path = './bin';
addpath(genpath(bin_path));

disp('Starting analyzing ...');

if ~exist('./network','dir')
    error('-- Directory "./network" is not found.');
end

% --- loading network ---
GeneNodeFileDir = './network/index_genes.txt';
NetworkFileDir = './network/edge_list.txt';

if ~exist(GeneNodeFileDir,'file')
    error(['-- File "' GeneNodeFileDir  '" is not found.']);
end
if ~exist(NetworkFileDir,'file')
    error(['-- File "' NetworkFileDir  '" is not found.']);
end

disp('Loading network files ...');
[net_map,Lap_mat] = PreprocessNetwork(GeneNodeFileDir,NetworkFileDir);

% --- configure ---
CompLeastProportion = 0.15;
maxComponent = 5;

NetConf.Lap_mat = Lap_mat;
NetConf.lambda_T = 0.1;
verbose = 1;
clear Lap_mat NetworkFileDir GeneNodeFileDir

if ~exist('./data','dir')
    error('-- Directory "./data" is not found.');
end


fid_list = dir('./data/*.txt');
num_File = length(fid_list);
if num_File == 0
    disp('-- No input data file found.');
else
    disp(['-- Totally ' num2str(num_File) ' data files found.']);
end
mkdir('./output');

for i_file = 1:num_File
    file_name_t = fid_list(i_file).name;
    disp(['-- -- File No.' num2str(num_File) ': ' file_name_t]);
    output_file_dir = ['./output/' file_name_t(1:(end-4))];
    disp(['Creating directory "' output_file_dir '" ...']);
    mkdir(output_file_dir);
    
    % read data
    [mutation_mat, sample_id, gene_id_symbol] = R01_read_gene_mat(['./data/' file_name_t]);

    [X_input,~,Symbol_Net] = A00_00_InputToNetMat(net_map,mutation_mat,gene_id_symbol,[bin_path '/GeneIDPreprocess'],1);
    clear dir_data mutation_mat

    % --- mCGfinder ---
    disp('Starting mCGfinder ...');
    tStart = tic;

    [S_sample_indicator,G_gene_score] = mCGfinder(X_input,NetConf,CompLeastProportion,maxComponent,verbose);

    disp('Significance Test...');
    [Q_values, ~] = SignificLayerTest(X_input,S_sample_indicator,G_gene_score,NetConf,verbose);
    disp('Significance Test done!');

    running_time = toc(tStart);
    disp(['Total Time of mCGfinder: ' num2str(running_time/60,'%2.2f') ' minutes' char(10)])
    clear X_input bkg_gene_ids CompLeastProportion tStart ...
        NetConf CompLeastProportion maxComponent verbose
    
    all_genes = net_map.Node2Gene_map.values;
    significance_scores = min(Q_values,[],2);
    [sort_signif,ind_sort] = sort(significance_scores);
    ind_th = find(sort_signif>=0.05,1) - 1;
    detected_genes = all_genes(ind_sort(1:ind_th));
    clear all_genes net_map ind_th sort_signif ind_sort significance_scores
    
    
    R02_00_output_vars(output_file_dir,detected_genes,S_sample_indicator,Symbol_Net,G_gene_score,Q_values,sample_id);
    
    save([output_file_dir '/Result.mat'],'detected_genes','S_sample_indicator','Symbol_Net','G_gene_score','Q_values');
    
end % end of file list for

rmpath(genpath(bin_path));
clear bin_path fid_list file_name_t gene_id_symbol i_file num_File output_file_dir sample_id ans

