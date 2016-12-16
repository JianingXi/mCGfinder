function [S_vec,G_vec,detectComponent] = mCGfinder(X,NetConf,CompLeastProportion,maxComponent,verbose)
% The discriptions of the configurations of mCGfindercan of these parameters are provided below:
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



[TotalNum, TotalGene] = size(X);

DefaultCLP = 0.15;
DefaultMC = 5;

if ~exist('NetConf','var')
    error('No network information input.')
end

if ~exist('CompLeastProportion','var')
    CompLeastProportion = DefaultCLP;
elseif CompLeastProportion <= 0 || CompLeastProportion >= 1
    disp('CompLeastProportion: out of domain. S_vecse default value instead.');
    CompLeastProportion = DefaultCLP;
end

if ~exist('maxComponent','var')
    maxComponent = DefaultMC;
elseif maxComponent <= 0
    disp('Max Component number: out of domain. S_vecse default value instead.');
    maxComponent = DefaultMC;
end

if ~exist('verbose','var')
    verbose = 0;
elseif verbose ~= 1
    disp('verbose number: out of domain. S_vecse default value instead.');
    verbose = 0;
end

leastGroupNum = max(ceil(CompLeastProportion*TotalNum),10);

S_vec = zeros(TotalNum,DefaultMC);
G_vec = zeros(TotalGene,DefaultMC);

selected_idx_sample = zeros(TotalNum,1);
% cur_idx_s = ones(TotalNum,1);

Remain_sample_num = TotalNum;

detectComponent = 0;
% && sum(cur_idx_s) >= leastGroupNum 
while sum(selected_idx_sample) < TotalNum && ...
        detectComponent < (maxComponent-1) && Remain_sample_num >= leastGroupNum
    % if all sample are selected,
    % or sample number of current comp is too small
    % or remain sample number is too small
    % or comp number detected is too large,
    % then break
    
    [s_vec_1,g_vec_1] = CurCompDecomp(X,NetConf,leastGroupNum,verbose);
    
    detectComponent = detectComponent + 1;
    
    % --- component samples --- %
    cur_idx_s = (s_vec_1~=0);
    S_vec_cur = zeros(TotalNum,1);
    S_vec_cur(selected_idx_sample==0)=s_vec_1;
    selected_idx_sample(selected_idx_sample==0) = cur_idx_s;

    G_vec_cur = g_vec_1;
    
    S_vec(:,detectComponent) = S_vec_cur;
    G_vec(:,detectComponent)  = G_vec_cur;
    X = X(~cur_idx_s,:);
    
    Remain_sample_num = sum(~selected_idx_sample);
    if verbose
        disp([char(10) 'Component ' num2str(detectComponent) ' done.']);
        disp(['Remain_sample_num: ' num2str(sum(~selected_idx_sample))]);
        disp(['Include_sample_num: ' num2str(sum(cur_idx_s))]);
    end
    
    if Remain_sample_num == 0
        return;
    end
end

% --- If remain samples exist---
if Remain_sample_num > 0
    Remain_sample_num = size(X,1);
    s_vec_1 = ones(Remain_sample_num,1);
    G_vec_cur = g_vec_Est(X,s_vec_1,NetConf);
    
    detectComponent = detectComponent + 1;
    
    cur_idx_s = (s_vec_1~=0);
    S_vec_cur = zeros(TotalNum,1);
    S_vec_cur(selected_idx_sample==0) = s_vec_1;
    selected_idx_sample(selected_idx_sample==0) = cur_idx_s;

    S_vec(:,detectComponent) = S_vec_cur;
    G_vec(:,detectComponent)  = G_vec_cur;

    if verbose
        disp([char(10) 'Component ' num2str(detectComponent) ' done.']);
        disp(['Remain_sample_num: ' num2str(sum(~selected_idx_sample))]);
        disp(['Include_sample_num: ' num2str(sum(cur_idx_s))]);
    end
end

ind_sample = (sum(abs(S_vec),1)>0);
ind_gene = (sum(abs(G_vec),1)>0);
S_vec = S_vec(:,ind_sample);
G_vec = G_vec(:,ind_gene);

end






function [s_vec,g_vec,iter] = CurCompDecomp(X,NetConf,leastGroupNum,verbose)

[n_samp,p_gene] = size(X);

if sum(sum(abs(X))) == 0
    s_vec = ones(n_samp,1)/n_samp;
    g_vec = zeros(p_gene,1);
    iter = 0;
    return;
end

try
    [s_vec_0,~,~] = svds(X,1);
    if sum(s_vec_0<-(10^-4)) > sum(s_vec_0>(10^-4))
        % ensure most element of u0 is positive
        s_vec_0 = -s_vec_0;
    end
    s_vec_0(s_vec_0<=0)=0;
    s_vec_0 = s_vec_0/max(s_vec_0);
    g_vec_0 = X'*s_vec_0;
catch
    s_vec_0 = ones(n_samp,1);
    g_vec_0 = sum(X,1)';
end


Sample_Number_in = 0;

relative_merr = 0.01;
niter = 100;

relative_s_vec_d = 2;
relative_g_vec_d = 2;
iter = 0;
while (relative_s_vec_d > relative_merr || relative_g_vec_d > relative_merr) ...
        || (Sample_Number_in < leastGroupNum) && iter < niter
    iter = iter+1;
    
    if verbose
        disp([char(10) 'Interation: ' num2str(iter)]);
    end
    
    % ---- G_vec updating ----
    g_vec_1 = g_vec_Est(X,s_vec_0,NetConf);
    
    
    % ---- S_vec updating ----
    s_vec_1 = S_vec_binary(X,g_vec_1,leastGroupNum);
    Sample_Number_in = sum(s_vec_1~=0);
    
    
    % ---- Residual ----
	relative_s_vec_d = sum((s_vec_0-s_vec_1).^2)/sum(s_vec_0.^2);
    relative_g_vec_d = sum((g_vec_0-g_vec_1).^2)/sum(g_vec_0.^2);
    
    if verbose ~= 0
        disp([char(9) 'Residuel s_vec diff: ' num2str(relative_s_vec_d,'%2.2e')]);
        disp([char(9) 'Residuel g_vec diff: ' num2str(relative_g_vec_d,'%2.2e')]);
    end
    s_vec_0 = s_vec_1;
	g_vec_0 = g_vec_1;
end

s_vec = s_vec_1;
g_vec = g_vec_1;

end

function g_vec_1_hat = g_vec_Est(X,s_vec_fix,NetConf)

GeneLen = size(X,2);

A_Mat_all = sum(s_vec_fix.^2)*eye(GeneLen) + (NetConf.lambda_T)*(NetConf.Lap_mat);

% y = A\x = inv(A)*x
g_vec_1_hat = A_Mat_all\(X'*s_vec_fix);

end



function s1_hat = S_vec_binary(X,g_vec_fix,leastGroupNum)
a = sum(g_vec_fix.^2);
b_vec = X*g_vec_fix;

vec_x = b_vec/a;

th_candi = sort(vec_x,'descend');
g_vec_th = min(th_candi(leastGroupNum),0.5);

s1_hat = (vec_x >= g_vec_th);
end

