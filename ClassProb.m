%INPUT:
%   z: can be a vector of atributes (Lx1) or a matrix 
%       with the vector atributes of M observations (LxM)
%   beta: vector (Lx1) parameters
%OUTPUT: 
%   p: if z is (Lx1), p is (1x1)the ClassProbabiliry of obs k
%      if Z is (LxM), p is (Mx1) the ClassProb of M obs

function p = ClassProb(Z,beta)
dim = size(Z,2);
switch dim
    case 1 %vector
        aux = beta'*Z;
    otherwise %matrix Z
        aux = (beta'*Z')';
end
e_aux = exp(aux);
p = (e_aux)./(1+e_aux);

idx = find(isnan(p));  p(idx) = 1; %in case some probability is inf/inf
end