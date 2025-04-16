function eff = get_efficiency(x,c,w_vect)
%
% efficiency of a design matrix
% x: design matrix
% c: contrast matrix
% w_vect: weight vector (how important each contrast is)
%
%

inv_xtx = inv(x'*x);
%c = [1 -1;1 0;0 1]
var_cBeta = c*inv_xtx*c';
%w_vect = [1 1 1];
diag_w = diag(w_vect);
weighted_var_cBeta = diag_w*var_cBeta;
% the diagonal elements of weighted_var_cBeta are the diagonal elements
% of var_cBeta weighted by w_vect
vari = trace(weighted_var_cBeta); % trace is the sum of the diagonal elements of a matrix
eff = 1/vari;