% INPUT:
%   beta: starting point of the coeficients
%   a: vector with the observed values of A
%   AL: (nxnxL)matrix with the powers of A or P in each dimension
%   nO: number of observed values
%   L: degree of the FIR filter
%   ind_o: indexes of the observed values
%   it_max: maximum iterations
% OUTPUT:
%   beta: beta_mle found with Newton Raphson algorithm
%   Z: nOxL matrix, feature vector of each observed edge in each row 

function [beta,Z] = MLE (beta,a,AL,nO,L,ind_o,it_max)

    Z = ones(nO,L);
    for j = 2:L
        AL_aux = AL(:,:,j);
        Z(:,j) = AL_aux(ind_o);
    end
    
    it=0;
    while it <= it_max  
        p = ClassProb(Z,beta);
        p_aux = p.*(1-p);
        W = spdiags(p_aux,0,nO,nO) ;

        z_pos = find(p_aux<=1e-20);
        p_aux(z_pos) = 1;
        inv_p_aux = 1./p_aux;
        inv_p_aux(z_pos) = 0;
        WW = spdiags(inv_p_aux,0,nO,nO) ;

        y = Z*beta+0.5.*(WW)*(a-p); %step size halving
        
        aux = sqrt(W)*Z;
        aux(isnan(aux)) = 0; aux(isinf(aux)) = 0; 
        
        beta_n = pinv(aux)*sqrt(W)*y;

        beta = beta_n;
        beta(1) = 0;
        it = it+1;
    end

end