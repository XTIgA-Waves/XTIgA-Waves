clc
clear
close all

% This code is associated with the paper
%       "Unconditionally stable space-time isogeometric discretization
%           for the wave equation in Hamiltonian formulation"
% by M. Ferrari, S. Fraschini, G. Loli and I. Perugia

% In this code we verifies (for p = 2,...,17) that the entries
% in the first 2p+1 rows and N1(p,10^(-13)) columns 
% and in the last 2p-1 rows and last N2(p,10^(-13)) columns 
% of the matrices
%
%           B_n^p*C_n^p*(B_n^p)^(-1)*C_n^p-(C_n^p)^2
%
% do not depend on the size of n (B_n^p and C_n^p defined in Equation (7))
% i.e. they are equal up to 10^(-13)

% This code accompanies the proof of Property 3.16.

% B_h^p and C_h^p have been previously computed with 1000 digits precision

%this numbers have been computed with the code Prop316_b.m
N1 = [0 20 31 39 47 55 61 68 73 79 83 87 91  94  97  99  101];
N2 = [0 23 34 44 53 61 69 76 82 88 93 98 102 106 110 112 115];

for p = 15:17

    p

    Nt = 2^8;

    cd ./C_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_C_sym'))
    cd ..
    cd ./B_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_B_sym'))
    cd ..

    B = B_sym;
    C = C_sym;

    X = B*C*(B\C)-C^2;

    X1_ref = X(1:2*p+1,1:N1(p));
    X2_ref = X(end-2*p+3:end,end-N2(p)+1:end);

    Nt = 2^7;

    cd ./C_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_C_sym'))
    cd ..
    cd ./B_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_B_sym'))
    cd ..

    B = B_sym/Nt;
    C = C_sym;

    X = B*C*(B\C)-C^2;

    X1 = X(1:2*p+1,1:N1(p));
    X2 = X(end-2*p+3:end,end-N2(p)+1:end);

    if norm(X1-X1_ref)>10^(-13) || norm(X2-X2_ref)>10^(-13)
        double(norm(X1-X1_ref))
        double(norm(X2-X2_ref))
        disp('problem')
    end
end