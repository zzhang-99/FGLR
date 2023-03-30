function [L] = FGLR_graph(img_all)
[M1  p1] = size(img_all);
edge_min  = 6;         % minimum degree of the graphs
c_eps = 50;      % control eps for constructing graphs
mask_diag = not(logical(eye(p1)));
eye_len = speye(p1);

%%
l1=1:1*p1;
F = l1'; % all feature functions
norm_node = sum(F .* F, 2);
A = repmat(norm_node, [1 p1]) + repmat(norm_node', [p1 1])-2 * (F * F');
A_sort = nth_element(A, edge_min + 1);
r = max(A_sort(edge_min + 1, :)) + 1e-6;
mask_tmp = triu(mask_diag&A <= r); % speed up with the upper-triangle characteristic
A(not(mask_tmp)) = 0;
A(mask_tmp) = exp(-A(mask_tmp) / (2 * c_eps^2));
A = sparse(A);
D = diag(sum(A) + sum(A, 2)'); % obtain the degree matrix
L = full(D - (A + A'));
end