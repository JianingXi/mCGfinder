function [Q_values,P_values] = SignificLayerTest(X,S_vec_final,G_vec_final,NetConf,verbose)

if ~exist('verbose','var')
    verbose = 0;
end
GeneLen = size(X,2);


[geneNum,k] = size(G_vec_final);
Q_values = zeros(geneNum,k);
P_values = Q_values;

for i = 1:k
    s1_cur = S_vec_final(:,i);
    cur_idx_u = (s1_cur~=0);
    
    A_Mat_all = sum(s1_cur.^2)*eye(GeneLen) + (NetConf.lambda_T)*(NetConf.Lap_mat);
    X_in = (A_Mat_all\X')';
    
    if verbose
        disp(['Permutation test of Comp' num2str(i) ' ...']);
    end
    [Q_values(:,i),P_values(:,i)] ...
        = Perm_Test(full(X_in(cur_idx_u,:)),full(s1_cur(cur_idx_u)),full(G_vec_final(:,i)),verbose);
end

end 




function [Q_values,P_values] = Perm_Test(X_hat,S_cur,G_cur,verbose)

GeneNum = length(G_cur);

if sum(G_cur~=0) == 0
    if verbose
        disp('All coefficients in v is 0.');
    end
    % v_cur contains only zero
    P_values = ones(GeneNum,1);
    Q_values = ones(GeneNum,1);
else

    P_values = Signif_Test_Weighted_conv(X_hat,S_cur,G_cur,verbose);
    Q_values = mafdr(P_values, 'BHFDR', 1);
end

end



function P_values = Signif_Test_Weighted_conv(X_hat,S_cur,G_cur,verbose)

Num_samp = length(S_cur);
X_all_raw = diag(S_cur)*X_hat;
X_all_abs = abs(X_all_raw);
C_bin = min(min(X_all_abs(X_all_abs>0)));
Max_all = max(max(X_all_abs));

dom_num = ceil(Max_all/C_bin);

if dom_num > 10^5
    dom_num = 10^5;
    C_bin = Max_all/dom_num;
end

for i_s = 1:Num_samp
    
    if ~mod(i_s,10) && verbose
        disp('Permuting ...');
    end
    
    X_hat_i_s = X_all_raw(i_s,:);
    
    max_dom = max(X_hat_i_s);
    dom_num = ceil(max_dom/C_bin);
    N_domain = dom_num+1;
    idx_domain = C_bin*(0:dom_num);
    
    if sum(X_hat_i_s ~= 0) == 0
        % all zeros case
        h2 = 1;
        d2 = 0;
    else
        h2_raw = hist(X_hat_i_s,idx_domain);
        
        p1 = find(h2_raw~=0,1);
        p2 = find(fliplr(h2_raw)~=0,1);
        idx2 = p1:(N_domain-p2+1);
        h2 = h2_raw(idx2);
        h2 = h2/sum(h2);
        d2 = idx_domain(idx2);
    end
    
    if i_s == 1
        h1 = h2;
        d1 = d2;
    else
        h3 = conv(h1,h2);
        nd3 = length(h3);
        
        if nd3 > 1
            maxd1 = max(d1);
            maxd2 = max(d2);
            mind1 = min(d1);
            mind2 = min(d2);
            slope = (maxd1+maxd2-mind1-mind2)/(nd3-1);
            d3 = (mind1+mind2)+slope*(0:(nd3-1));
        else
            d3 = d2;
        end
        
        h1 = h3;
        d1 = d3;
    end
end

cdf_Up = fliplr(cumsum(fliplr(h1)));
d_Up = d1;

G_ind = ceil(G_cur/C_bin);
G_ind(G_ind == 0) = 1;
G_ind(G_ind>length(d_Up)) = length(d_Up);
P_values = cdf_Up(G_ind)';

end

