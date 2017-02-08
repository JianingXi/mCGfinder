function R02_00_output_vars(output_file_dir,detected_genes,S_sample_indicator,Symbol_Net,G_gene_score,Q_values,sample_id)

RowN = size(Symbol_Net,1);
if RowN == 1
    Symbol_Net = Symbol_Net';
end
RowN = size(detected_genes,1);
if RowN == 1
    detected_genes = detected_genes';
end

comp_n = size(S_sample_indicator,2);

% --Q_values--
num_vec = num2cell(1:comp_n);
col_name = cellfun(@(x) ['Comp' num2str(x,'%.2d')],num_vec,'UniformOutput',0);
row_name = Symbol_Net;
R02_01_write_table_vars([output_file_dir '/' 'Q_values' '.txt'],Q_values,'GeneSymbol',col_name,row_name,'%.2e');

% --G_gene_score--
num_vec = num2cell(1:comp_n);
col_name = cellfun(@(x) ['Comp' num2str(x,'%.2d')],num_vec,'UniformOutput',0);
row_name = Symbol_Net;
R02_01_write_table_vars([output_file_dir '/' 'G_gene_score' '.txt'],G_gene_score,'GeneSymbol',col_name,row_name,'%.2e');

% --S_sample_indicator--
num_vec = num2cell(1:comp_n);
col_name = cellfun(@(x) ['Comp' num2str(x,'%.2d')],num_vec,'UniformOutput',0);
row_name = sample_id;
R02_01_write_table_vars([output_file_dir '/' 'S_sample_indicator' '.txt'],S_sample_indicator,'SampleID',col_name,row_name,'%.1d');


% --Symbol_Net--
T_print = table(Symbol_Net,'VariableNames',{'GeneSymbol'});
writetable(T_print,[output_file_dir '/' 'Symbol_Net' '.txt'],'Delimiter','\t')

% --detected_genes--
T_print = table(detected_genes,'VariableNames',{'GeneSymbol'});
writetable(T_print,[output_file_dir '/' 'detected_genes' '.txt'],'Delimiter','\t')

end

