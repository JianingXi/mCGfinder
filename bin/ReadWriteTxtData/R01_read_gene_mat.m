function [mutation_mat, sample_id, gene_id_symbol] = R01_read_gene_mat(input_txt_file_str)

try
    result_struct = importdata(input_txt_file_str);
    mutation_mat = result_struct.data;
    [n_gene,n_sample] = size(mutation_mat);

    [n_gene_1,n_sample_1] = size(result_struct.textdata);
    if n_gene + 1 ~= n_gene_1
        error('Mutation data file: rows mismatch.');
    end
    if n_sample + 1 ~= n_sample_1
        error('Mutation data file: columns mismatch.');
    end

    sample_id = (result_struct.textdata(1,2:end))';
    gene_id_symbol = result_struct.textdata(2:end,1);
    mutation_mat = mutation_mat';
catch
    error('Mutation data file: format error!');
end

end