clc
clear
close all
format long

% This code is associated with the paper
%       "Unconditionally stable space-time isogeometric discretization
%           for the wave equation in Hamiltonian formulation"
% by M. Ferrari, S. Fraschini, G. Loli and I. Perugia


% In this code we verify (for p = 2,...,25) the invertibility of the
% matrices in the family $\{C_n^p\}_n$ with C_h^p=C_n^p defined in 
% Equation (7). The invertibility is verified as in Theorem 3.4 ii)
% by checking the invertibility of an associated matrices with size
% not depending on n

% This code accompanies the proof of Proposition 3.15.

% C_h^p have been previously computed with 1000 digits precision


syms rho

digits(1000)

for p = 19
    for Nt = [2^7 2^8]

        m = p+1;
        ell = p-1;

        cd ./C_matrices
        load(strcat(num2str(p),'_',num2str(Nt),'_C_sym'))
        cd ..

        
        P = C_sym(2*p,p-1:p-1+2*p);
        z = roots(P);
        [I,M] = sort(abs(z));
        z = vpa(z(M));

        W = vpa(zeros(2*p,2*p));

        for i = 1:2*p
            for j = 1 : 2*p
                W(j,i) = vpa(z(i)^(j-1));
            end
        end

        Y2 = C_sym(1:m+ell,ell+1:m+2*ell);
        Y1 = C_sym(1:m+ell,1:ell);

        inv_W = inv(W);
        min_singular_value = min(double(svd(inv_W(m+1:m+ell,1:m+ell)*inv(Y2)*Y1)))

    end
end

