clc
clear
close all

% This code is associated with the paper
%       "Unconditionally stable space-time isogeometric discretization
%           for the wave equation in Hamiltonian formulation"
% by M. Ferrari, S. Fraschini, G. Loli and I. Perugia

% In this code we compute (for p = 2,...,17) the numbers N1(p,eps) and 
% N2(p,eps) such that the entries of the matrices
% 
%           B_n^p*C_n^p*(B_n^p)^(-1)*C_n^p-(C_n^p)^2 
% 
% in the first 2p+1 rows larger than eps are in the first N1(p,eps) columns
% in the last 2p-2 rows larger than eps are in the last N2(p,eps) columns

% This code accompanies the proof of Property 3.16.

% B_n^p and C_n^p have been previously computed with 1000 digits precision

eps = 10^(-13);

p_min = 2;
p_max = 17;

N1 = zeros(p_max-p_min+1,1);
N2 = zeros(p_max-p_min+1,1);

for p = p_min:p_max

    cont = 1;

    for Nt = [2^7 2^8]

        cd ./C_matrices
        load(strcat(num2str(p),'_',num2str(Nt),'_C_sym'))
        cd ..
        cd ./B_matrices
        load(strcat(num2str(p),'_',num2str(Nt),'_B_sym'))
        cd ..

        B = B_sym;
        C = C_sym;
        
        X = B*C*(B\C)-C^2;

        X1 = X(1:2*p+1,:);
        X2 = X(end-2*p+3:end,:);

        for i = 1 : Nt+p-1
            if max(max(abs(X1(:,i+1:end))))<eps
                N1(p-p_min+1,cont) = i;
                break
            end
        end

        for i = Nt+p-1 : -1 : 1
            if max(max(abs(X2(:,1:i))))<eps
                N2(p-p_min+1,cont) = Nt+p-1-i;
                break
            end
        end
        cont = cont+1;
    end
    
    deg_p = [p_min:p_max]';
    N1_128 = N1(:,1);
    N1_256 = N1(:,2);
    N2_128 = N2(:,1);
    N2_256 = N2(:,2);

    table(deg_p,N1_128,N1_256,N2_128,N2_256)
end
