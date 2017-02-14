# mCGfinder
A novel network regularized matrix decomposition method to detect mutated cancer genes in tumour samples with inter-patient heterogeneity

=======================
Instructions to mCGfinder software (version 1.0.0)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China


Requirement
------------------------
* 4GB memory
* MATLAB R2015a or later

Gene interaction network information
------------------------
The default 
If you want to use the default gene interaction network [iRefIndex 9] (http://irefindex.org), just ignore this step. Otherwise, if you want to use a user-defined network, replace the two files below with the files of the user-defined network:

* `./network/index_genes.txt`
* `./network/edge_list.txt`

The first text file `index_genes.txt` is a table of the connection between the node id and the related gene names. The second text file `edge_list.txt` is the edges table that assigns the source nodes to the target nodes.


Run mCGfinder on somatic mutation data file
------------------------

We provide an example data file of somatic mutations in breast invasive carcinoma (BRCA) samples from [UCSC Cancer Genomics Browser] (https://genome-cancer.soe.ucsc.edu/proj/site/hgHeatmap/) in `./data/example_data.zip`. To analyze this data, please extract the .txt file `somatic_data_BRCA.txt` from `./example_data.zip` and put the .txt file in the `./data` folder. If you want to analyze the example data file with the default configurations, please run `./demo_mCGfinder.m` and then the result file will be automatically saved in the directory `./output/`.

If you want to analyze a user-specific data, a .txt file of the mutation binary table (samples x genes) of the sample names and the gene symbols must be provided as the format of the example data file. Put the txt file of the user-specific data in the `./data` folder and run `./demo_mCGfinder.m`.


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


mCGfinder standalone version
------------------------

For users without MATLAB licenses, we also offer mCGfinder standalone version for Windows `./mCGfinder_standalone.zip`.

* Step 1: MATLAB runtime installer: verify the MATLAB runtime is installed and ensure you have installed version 8.5 (R2015a). If the MATLAB runtime is not installed, please download the Windows 64-bit version of the MATLAB runtime for R2015a from the MathWorks Web site by navigating to
 
   [http://www.mathworks.com/products/compiler/mcr/index.html] (http://www.mathworks.com/products/compiler/mcr/index.html)
   
  Or run the installation file `./mCGfinder_standalone/MyAppInstaller_web.exe` provided in the .zip file. 

* Step 2: Extract the `./mCGfinder_standalone.zip` file, and locate the executable file `./mCGfinder_standalone/mCGfinder.exe` and the two folders `./mCGfinder_standalone/data` and `./mCGfinder_standalone/network` in the same directory `./mCGfinder_standalone/`. Put the .txt files of mutation data in folder `./mCGfinder_standalone/data/` and network files `index_genes.txt` and `edge_list.txt` in folder `./mCGfinder_standalone/network/`.

* Step 3: Run `./mCGfinder_standalone/mCGfinder.exe`. Please wait until the current program is finished, and the output variables are automatically saved as .txt files in folder `./mCGfinder_standalone/output`.


Contact
------------------------
Please feel free to contact us if you need any help: <xjn@mail.ustc.edu.cn>.