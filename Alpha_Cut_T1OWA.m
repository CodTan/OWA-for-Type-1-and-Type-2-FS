%% To run this code, change the data and the alpha cuts depending on how many required and their values
% Also change the number of matrices required to store the data from cell
% arrays to calculate the end points of the aggregated set

clc;
close all
clear

%% ************************* Define data *********************************

W_dom = 0.0:0.5:1.0;
A_dom = 0.0:1.0:2.0;
A = [];
W = [];
% A = [trimf(A_dom,[0.0 2.0 2.0]); trimf(A_dom,[0.0 0.0 2.0]); trimf(A_dom,[0.0 1.0 2.0])];
A = [trimf(A_dom,[0.0 0.0 2.0]); trimf(A_dom,[0.0 1.0 2.0])];
% W = [trimf(W_dom,[0.0 0.0 1.0]); trimf(W_dom,[0.0 0.5 1.0]); trimf(W_dom,[0.0 1.0 1.0])];
W = [trimf(W_dom,[0.0 0.0 1.0]); trimf(W_dom,[0.0 1.0 1.0])];

A_temp = [0.0 0.5 1.0];
W_temp = [0.0 1.0 0.0];

figure
subplot(2,1,1);
plot(A_dom,A,'b');
hold on;
plot(A_dom,A_temp,'b');
subplot(2,1,2);
plot(W_dom,W,'r');
hold on;
plot(W_dom,W_temp,'r');

A = [A_temp;[A]];
W = [trimf(W_dom,[0.0 0.0 1.0]); 0.0 1.0 0.0; trimf(W_dom,[0.0 1.0 1.0])];

alpha = [0.0,0.5,1.0];
A_alpha_cut_sets = {};
W_alpha_cut_sets = {};
A_minus = [];
A_plus = [];

W_minus = [];
W_plus = [];

lr_int = {};
G = {}; % aggregated set points
mu_G = []; % aggregated set membership

%% ************************* alpha-cuts *********************************

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
    
    for i = 1:size(A,1)
        for j = 1:length(A_dom)
            if(A(i,j) >= alpha(k))
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
    
    for i = 1:size(W,1)
        for j = 1:length(W_dom)
            if(W(i,j) >= alpha(k))
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
    
    
    for l = 1:size(A,1) % for rho
        for m = l:size(A,1)
            denom_sum = denom_sum + W_plus(m);
            numer_sum = numer_sum + W_plus(m)*A_minus(sigma(m));
        end
        if(l>1)
            for n = 1:(size(A,1)-l)
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
    
    
    for l = 1:size(A,1) % for rho
        for m = l:size(A,1)
            denom_sum = denom_sum + W_minus(m);
            numer_sum = numer_sum + W_minus(m)*A_plus(sigma(m));
        end
        if(l>1)
            for n = 1:(size(A,1)-l)
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
    
    left_right_interval = [rho_minus:rho_plus];
    lr_int{end+1} = [rho_minus:rho_plus];
    G{end+1} = intersect(left_right_interval,A_dom);
    
    A_minus = [];
    A_plus = [];
    W_minus = [];
    W_plus = [];
    A_alpha_cell = {};
    W_alpha_cell = {};    
end
