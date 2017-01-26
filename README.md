# mCGfinder
A novel network regularized matrix decomposition method to detect mutated cancer genes in tumour samples with inter-patient heterogeneity

=======================
Instructions to mCGfinder software (version 1.0.0)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China


Requirement
------------------------
* 4GB memory
* MATLAB 2013a or later

Gene interaction network information
------------------------
The default 
If you want to use the default gene interaction network [iRefIndex 9] (http://irefindex.org), just ignore this step. Otherwise, if you want to use a user-defined network, replace the two files below with the files of the user-defined network:

* `./network/index_genes`
* `./network/edge_list`

The first text file `index_genes` is a table of the connection between the node id and the related gene names. The second text file `edge_list` is the edges table that assigns the source nodes to the target nodes.


Run mCGfinder on somatic mutation data files
------------------------
In the `./data` directory, we provide example data files of somatic mutations in breast invasive carcinoma (BRCA) samples, bladder urothelial carcinoma (BLCA), glioblastoma multiforme (GBM), head and neck squamous cell carcinoma (HNSC) from [UCSC Cancer Genomics Browser] (https://genome-cancer.soe.ucsc.edu/proj/site/hgHeatmap/). If you want to analyze the example data files with the default configurations, please run `./demo_mCGfinder.m` and then the result file will be saved in the directory `./output/` as `./output/Result.mat`.

If you want to analyze a user-specific data, the mutation binary matrix (samples x genes), the sample names and the gene symbols of the samples must be provided. Save the three variables as 'mutation_mat', 'sample_id' and 'gene_id_symbol' as the example data files in the `./data` directory . Then update the string variable 'dir_data' in file `./demo_mCGfinder.m` with the directory of the new data file.


Configurations of mCGfinder 
------------------------

The configurations of mCGfinder can be changed in script file `./demo_mCGfinder.m`, and the descriptions of these parameters are provided below:

        =================================================================================================
        | PARAMETER NAME       | DESCRIPTION                                                            |
        =================================================================================================
        |CompLeastProportion   |Least sample proportion included in each components, which represents   |
        |                      |minimum proportion of the samples in every components given by the      |
        |                      |mCGfinder. The default proportion is set to 15%.                        |
        -------------------------------------------------------------------------------------------------
        |maxCompoent           |Maximum number of components, which denotes the number of components    |
        |                      |given components given by mCGfinder at most. The default number is 5.   |
        -------------------------------------------------------------------------------------------------
        |NetConf.lambda_T      |The tuning parameter of network regularization, which is used to balance|
        |                      |the fitness of the model (first term) and the smoothness of the scores  |
        |                      |of connected genes (second term). The default number is 0.1.            |
        -------------------------------------------------------------------------------------------------


Output variables of mCGfinder 
------------------------

The descriptions of output variables of mCGfinder are provided below:

        =================================================================================================
        | VARIABLE NAME        | DESCRIPTION                                                            |
        =================================================================================================
        |detected_genes        |Genes detected by mCGfinder as significantly mutated cancer genes.      |
        |                      |minimum proportion of the samples in every components given by the      |
        |                      |mCGfinder. The default proportion is set to 15%.                        |
        -------------------------------------------------------------------------------------------------
        |S_sample_indicator    |The sample indicator vectors of all component, which indicates the      |
        |                      |assignment of tumour samples to the every components. The i-th          |
        |                      |coefficient being 1 represents that the i-th samples are included in the|
        |                      |component, and 0 otherwise.                                             |
        -------------------------------------------------------------------------------------------------
        |Symbol_Net            |The investigated gene list in the gene interaction network.             |
        -------------------------------------------------------------------------------------------------
        |G_gene_score          |The gene score vectors of all components, of which the coefficients are |
        |                      |related to the gene lists variable 'Symbol_Net', and a higher value of  |
        |                      |a certain coefficient presents a larger potential of the gene to be     |
        |                      |cancer gene candidate.                                                  |
        -------------------------------------------------------------------------------------------------
        |Q_values              |The q-values of all investigated genes in variable 'Symbol_Net', which  |
        |                      |are obtained by Benjamini-Hochberg false discovery rates control of the |
        |                      |p-values of the investigated genes.                                     |
        -------------------------------------------------------------------------------------------------

Contact
------------------------
Please feel free to contact us if you need any help: <xjn@mail.ustc.edu.cn>.