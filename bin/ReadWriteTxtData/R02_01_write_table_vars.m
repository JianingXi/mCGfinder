function R02_01_write_table_vars(output_file_str,d_print_variable,head_str,col_name,row_name,format_str)

if isnumeric(d_print_variable)
    x_cell = num2cell(d_print_variable);
    x_print = cellfun(@(x) num2str(x,format_str),x_cell,'UniformOutput',0);
else
    x_print = d_print_variable;
end

comp_n = size(d_print_variable,2);
num_vec = num2cell(1:comp_n);

str_print_in = cellfun(@(x) ['x_print(:,' num2str(x) ')'],num_vec,'UniformOutput',0);
eval(['T_print = table(row_name,' strjoin(str_print_in,', ') ',''VariableNames'',[head_str col_name]);']);

writetable(T_print,output_file_str,'Delimiter','\t')
end

