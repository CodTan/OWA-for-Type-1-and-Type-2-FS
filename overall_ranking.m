clear
clc

D_by_rank = [1 2 3 4;
             3 2 4 1;
             2 3 1 4;
             4 2 3 1;
             1 3 4 2;
             1 4 3 2];
         
D_w = [0.25 0.30 0.15 0.20 0.5 0.5];

num_of_ranks = 1:size(D_by_rank,2);
r = sort(num_of_ranks,'descend');
p = 1:length(num_of_ranks);
r(p) = r;
priority = r;

D_w_sum_final = [];
[r,c] = size(D_by_rank);
D_w_sum = zeros(1,c);
temp_sum = 0;

D_w_sum_matrix_element_wise = [];

for j = 1:c
    
    X = D_by_rank(:,j);
    unique_el = unique(X);
%     unique_el = unique_el';
    freq = [unique_el,histc(X(:),unique_el)];
    
    for i = 1:r
        D_w_sum(D_by_rank(i,j)) = D_w_sum(D_by_rank(i,j)) + D_w(i);
    end
    
    for k = 1:size(freq,1)
        D_w_sum(freq(k,1)) = D_w_sum(freq(k,1)) / freq(k,2);
    end
    
    D_w_sum_matrix_element_wise = [D_w_sum_matrix_element_wise; D_w_sum];
    
    D_w_sum_final = [D_w_sum_final; (priority(j)*D_w_sum)];
    
    D_w_sum = zeros(1,c);
end

element_overall_ranking = sum(D_w_sum_final,1);

element_wise_sum = sum(D_w_sum_matrix_element_wise);

element_overall_ranking = element_overall_ranking ./ element_wise_sum;

[~,p] = sort(element_overall_ranking,'descend');
r = 1:length(element_overall_ranking);
r(p) = r;
ranked_list = r;



