%% To run this code, change the data and the alpha cuts depending on how many required and their values
% Also change the number of matrices required to store the data from cell
% arrays to calculate the end points of the aggregated set

clc;
close all
clear

%% ************************* Define data *********************************

W_dom = 0.0:0.2:1.0;
A_dom = 0.0:1.0:5.0;

A_UMF = [1 1 0.5 0 0 0; 0 0.5 1 1 0.5 0; 0 0 0 0.5 1 1];
W_UMF = [1 0.5 0 0 0 0; 0 0 0.5 1 0 0; 0 0 0 0 0.5 1];

A_LMF = [0.8 0.8 0 0 0 0;0 0 0.8 0.8 0 0; 0 0 0 0 0.8 0.8];
W_LMF = [1 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1];

LT_1_xUMF=[0 1 2 3];
LT_1_uUMF=[1 1 0.5 0];
LT_1_xLMF=[0 0 1 2];
LT_1_uLMF=[0.8 0.8 0.8 0];

figure
subplot(3,1,1);
plotIT2(LT_1_xUMF,LT_1_uUMF,LT_1_xLMF,LT_1_uLMF,[0,5]); 
hold on;

LT_2_xUMF=[0 2 3 5];
LT_2_uUMF=[0 1 1 0];
LT_2_xLMF=[1 2  3 4];
LT_2_uLMF=[0 0.8  0.8 0];

subplot(3,1,1);
plotIT2(LT_2_xUMF,LT_2_uUMF,LT_2_xLMF,LT_2_uLMF,[0,5]); 
hold on;

LT_3_xUMF=[2 4 5 5];
LT_3_uUMF=[0 1 1 1];
LT_3_xLMF=[3 4 5 5];
LT_3_uLMF=[0 0.8 0.8 0.8];

subplot(3,1,1);
plotIT2(LT_3_xUMF,LT_3_uUMF,LT_3_xLMF,LT_3_uLMF,[0,5]); 
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_1_xUMF=[0 0 0.4];
w_1_uUMF=[1 1 0];
w_1_xLMF=[0 0 0.2];
w_1_uLMF=[1 1  0];

subplot(3,1,2);
plotIT2(w_1_xUMF,w_1_uUMF,w_1_xLMF,w_1_uLMF,[0,1]); 
hold on;

w_2_xUMF=[0.2 0.6 0.8];
w_2_uUMF=[0 1 0];
w_2_xLMF=[0.4 0.6  0.6];
w_2_uLMF=[0 1  0];

subplot(3,1,2);
plotIT2(w_2_xUMF,w_2_uUMF,w_2_xLMF,w_2_uLMF,[0,1]); 
hold on;

w_3_xUMF=[0.6 1 1];
w_3_uUMF=[0 1 1];
w_3_xLMF=[0.8 1 1];
w_3_uLMF=[0 1 1];

subplot(3,1,2);
plotIT2(w_3_xUMF,w_3_uUMF,w_3_xLMF,w_3_uLMF,[0,1]); 
hold on;

alpha = [0.0,0.5,1.0];
A_alpha_cut_sets = {};
W_alpha_cut_sets = {};
A_minus = [];
A_plus = [];

W_minus = [];
W_plus = [];

lr_int = {};
G_UMF = {}; % aggregated set points
mu_G_UMF = []; % aggregated set membership

G_LMF = {}; % aggregated set points
mu_G_LMF = []; % aggregated set membership

%% *********************************UMF***********************************
%************************* alpha-cuts *********************************

for k = 1:size(alpha,2)
    
    A_alpha_cell = {};
    W_alpha_cell = {};
    A_alpha_temp = [];
    W_alpha_temp = [];
    sigma = [];
    rho_minus = 0;
    rho_plus= 0;
    
    % A_alpha_count = [];
    % count = 0;
    
    for i = 1:size(A_UMF,1)
        for j = 1:length(A_UMF)
            if(A_UMF(i,j) >= alpha(k))
                A_alpha_temp = [A_alpha_temp,A_dom(j)];
                %             count = count+1;
            end
        end
        A_alpha_cell{end+1} = A_alpha_temp;
        A_minus = [A_minus;min(A_alpha_temp)];
        A_plus = [A_plus;max(A_alpha_temp)];
        %     A_alpha_count(i) = count;
        %     count = 0;
        A_alpha_temp = [];
    end
    
    for i = 1:size(W_UMF,1)
        for j = 1:length(W_UMF)
            if(W_UMF(i,j) >= alpha(k))
                W_alpha_temp = [W_alpha_temp,W_dom(j)];
                %             count = count+1;
            end
        end
        W_alpha_cell{end+1} = W_alpha_temp;
        W_minus = [W_minus;min(W_alpha_temp)];
        W_plus = [W_plus;max(W_alpha_temp)];
        %     A_alpha_count(i) = count;
        %     count = 0;
        W_alpha_temp = [];
    end
    
    % ***************** Calculation of left end point *********************
    [~,sigma] = sort(A_minus,'descend');
    temp = 1:length(A_minus);
    sigma(temp) = sigma;
    denom_sum = 0;
    numer_sum = 0;
    
    
    for l = 1:size(A_UMF,1) % for rho
        for m = l:size(A_UMF,1)
            denom_sum = denom_sum + W_plus(m);
            numer_sum = numer_sum + W_plus(m)*A_minus(sigma(m));
        end
        if(l>1)
            for n = 1:(l-1)
                denom_sum = denom_sum + W_minus(n);
                numer_sum = numer_sum + W_minus(n)*A_minus(sigma(n));
            end
        end
%         for o = 1:size(A,1)
%             numer_sum = numer_sum + W_plus(o)*A_minus(sigma(o));
%         end
        rho_minus = numer_sum / denom_sum;
        if(rho_minus >= A_minus(sigma(l)))
            break;
        else
            denom_sum = 0;
            numer_sum = 0;
            continue;
        end
        denom_sum = 0;
        numer_sum = 0;
    end
    
    sigma = [];
    temp = [];   
     
    % ***************** Calculation of right end point *********************
    [~,sigma] = sort(A_plus,'descend');
    temp = 1:length(A_plus);
    sigma(temp) = sigma;
    denom_sum = 0;
    numer_sum = 0;
    
    
    for l = 1:size(A_UMF,1) % for rho
        for m = l:size(A_UMF,1)
            denom_sum = denom_sum + W_minus(m);
            numer_sum = numer_sum + W_minus(m)*A_plus(sigma(m));
        end
        if(l>1)
            for n = 1:(l-1)
                denom_sum = denom_sum + W_plus(n);
                numer_sum = numer_sum + W_plus(n)*A_plus(sigma(n));
            end
        end
%         for o = 1:size(A,1)
%             numer_sum = numer_sum + W_minus(o)*A_plus(sigma(o));
%         end
        rho_plus = numer_sum / denom_sum;
        if(rho_plus >= A_plus(sigma(l)))
            break;
        else
            denom_sum = 0;
            numer_sum = 0;
            continue;
        end
        denom_sum = o;
        numer_sum = 0;
    end
    
    rho_minus = round(rho_minus);
    rho_plus = round(rho_plus);
    left_right_interval = [rho_minus:rho_plus];
    lr_int{end+1} = [rho_minus:rho_plus];
    G_UMF{end+1} = intersect(left_right_interval,A_dom);
    
    A_minus = [];
    A_plus = [];
    W_minus = [];
    W_plus = [];
    A_alpha_cell = {};
    W_alpha_cell = {};    
end

U = union(G_UMF{1},G_UMF{2});
for i = 1:length(G_UMF)
    U = union(U,G_UMF{i});
end

temp_mu = [];

for j = 1:length(U)
    for i = 1:length(alpha)
        if(ismembertol(U(j),G_UMF{i},0.001))
            temp_mu = [temp_mu alpha(i)];
        end
    end
    mu_G_UMF = [mu_G_UMF,max(temp_mu)];
    temp_mu = [];
end

alpha_LMF = [0.0,0.5,0.8];
A_alpha_cut_sets = {};
W_alpha_cut_sets = {};
A_minus = [];
A_plus = [];

W_minus = [];
W_plus = [];

G_LMF = {}; % aggregated set points
mu_G_LMF = []; % aggregated set membership


%% *********************************LMF***********************************
%************************* alpha-cuts *********************************

for k = 1:size(alpha,2)
    
    A_alpha_cell = {};
    W_alpha_cell = {};
    A_alpha_temp = [];
    W_alpha_temp = [];
    sigma = [];
    rho_minus = 0;
    rho_plus= 0;
    
    % A_alpha_count = [];
    % count = 0;
    
    for i = 1:size(A_LMF,1)
        for j = 1:length(A_LMF)
            if(A_LMF(i,j) >= alpha_LMF(k))
                A_alpha_temp = [A_alpha_temp,A_dom(j)];
                %             count = count+1;
            end
        end
        A_alpha_cell{end+1} = A_alpha_temp;
        A_minus = [A_minus;min(A_alpha_temp)];
        A_plus = [A_plus;max(A_alpha_temp)];
        %     A_alpha_count(i) = count;
        %     count = 0;
        A_alpha_temp = [];
    end
    
    for i = 1:size(W_LMF,1)
        for j = 1:length(W_LMF)
            if(W_LMF(i,j) >= alpha_LMF(k))
                W_alpha_temp = [W_alpha_temp,W_dom(j)];
                %             count = count+1;
            end
        end
        W_alpha_cell{end+1} = W_alpha_temp;
        W_minus = [W_minus;min(W_alpha_temp)];
        W_plus = [W_plus;max(W_alpha_temp)];
        %     A_alpha_count(i) = count;
        %     count = 0;
        W_alpha_temp = [];
    end
    
    % ***************** Calculation of left end point *********************
    [~,sigma] = sort(A_minus,'descend');
    temp = 1:length(A_minus);
    sigma(temp) = sigma;
    denom_sum = 0;
    numer_sum = 0;
    
    
    for l = 1:size(A_LMF,1) % for rho
        for m = l:size(A_LMF,1)
            denom_sum = denom_sum + W_plus(m);
            numer_sum = numer_sum + W_plus(m)*A_minus(sigma(m));
        end
        if(l>1)
            for n = 1:(l-1)
                denom_sum = denom_sum + W_minus(n);
                numer_sum = numer_sum + W_minus(n)*A_minus(sigma(n));
            end
        end
%         for o = 1:size(A,1)
%             numer_sum = numer_sum + W_plus(o)*A_minus(sigma(o));
%         end
        rho_minus = numer_sum / denom_sum;
        if(rho_minus >= A_minus(sigma(l)))
            break;
        else
            denom_sum = 0;
            numer_sum = 0;
            continue;
        end
        denom_sum = 0;
        numer_sum = 0;
    end
    
    sigma = [];
    temp = [];   
     
    % ***************** Calculation of right end point *********************
    [~,sigma] = sort(A_plus,'descend');
    temp = 1:length(A_plus);
    sigma(temp) = sigma;
    denom_sum = 0;
    numer_sum = 0;
    
    
    for l = 1:size(A_LMF,1) % for rho
        for m = l:size(A_LMF,1)
            denom_sum = denom_sum + W_minus(m);
            numer_sum = numer_sum + W_minus(m)*A_plus(sigma(m));
        end
        if(l>1)
            for n = 1:(l-1)
                denom_sum = denom_sum + W_plus(n);
                numer_sum = numer_sum + W_plus(n)*A_plus(sigma(n));
            end
        end
%         for o = 1:size(A,1)
%             numer_sum = numer_sum + W_minus(o)*A_plus(sigma(o));
%         end
        rho_plus = numer_sum / denom_sum;
        if(rho_plus >= A_plus(sigma(l)))
            break;
        else
            denom_sum = 0;
            numer_sum = 0;
            continue;
        end
        denom_sum = o;
        numer_sum = 0;
    end
    
    rho_minus = round(rho_minus);
    rho_plus = round(rho_plus);
    left_right_interval = [rho_minus:rho_plus];
    lr_int{end+1} = [rho_minus:rho_plus];
    G_LMF{end+1} = intersect(left_right_interval,A_dom);
    
    A_minus = [];
    A_plus = [];
    W_minus = [];
    W_plus = [];
    A_alpha_cell = {};
    W_alpha_cell = {};    
end

L = union(G_LMF{1},G_LMF{2});
for i = 1:length(G_LMF)
    L = union(L,G_LMF{i});
end

temp_mu = [];

for j = 1:length(L)
    for i = 1:length(alpha)
        if(ismembertol(L(j),G_LMF{i},0.001))
            temp_mu = [temp_mu alpha(i)];
        end
    end
    mu_G_LMF = [mu_G_LMF,max(temp_mu)];
    temp_mu = [];
end

subplot(3,1,3)
plotIT2(U,mu_G_UMF,L,mu_G_LMF,[0,5]);
