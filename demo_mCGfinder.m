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
% The configurations of mCGfindercan be changed in script file `./demo_mCGfinder.m`, and the discriptions
% of these parameters are provided below:
% 
%         ========================================================================================
%         | PARAMETER NAME          | DESCRIPTION                                                |
%         ========================================================================================
%         |CompLeastProportion      |Least sample proportion included in each components. The    |
%         |                         |default proportion is set to 15%.                           |
%         ----------------------------------------------------------------------------------------
%         |maxCompoent              |Maximum number of components. The default number is 5.      |
%         ----------------------------------------------------------------------------------------
%         |NetConf.lambda_T         |Tuning parameter . The default number is 5.                 |
%         ----------------------------------------------------------------------------------------


bin_path = './bin';
addpath(genpath(bin_path));

% --- loading network ---
GeneNodeFileDir = './network/index_genes.txt';
NetworkFileDir = './network/edge_list.txt';
[net_map,Lap_mat] = PreprocessNetwork(GeneNodeFileDir,NetworkFileDir);

% --- configure ---
CompLeastProportion = 0.15;
maxComponent = 5;

NetConf.Lap_mat = Lap_mat;
NetConf.lambda_T = 0.1;
verbose = 1;


    
% --- format ---
dir_data = './data/somatic_data_BRCA.mat';   % 'BLCA', 'GBM', 'HNSC'
load(dir_data);
X_in = gene_indiv_mat;

[X_Net,~,Net_Idx,~,Symbol_Net,~,bkg_gene_ids] = ...
    A00_00_InputToNetMat(net_map,X_in,gene_id_symbol,[bin_path '/GeneIDPreprocess'],1);


% --- mCGfinder ---
disp('Run mCGfinder ...')
tStart = tic;

[S_vec_final,G_vec_final] = mCGfinder(X_Net,NetConf,CompLeastProportion,maxComponent,verbose);

disp('Significance Test...')
[Q_values, P_values] = SignificLayerTest(X_Net,S_vec_final,G_vec_final,NetConf);
disp('Significance Test done!')

running_time = toc(tStart);
disp(['Total Time of NetSSVD: ' num2str(running_time/60,'%2.2f') ' minutes' char(10)])


all_genes = net_map.Node2Gene_map.values;
detected_genes = all_genes((min(Q_values,[],2)<0.05));
rmpath(genpath(bin_path));

save('./output/Result.mat',...
    'X_Net','S_vec_final','G_vec_final','Q_values','P_values',...
    'detected_genes','Symbol_Net','bkg_gene_ids','running_time');




