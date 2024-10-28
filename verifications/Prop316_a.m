clc
clear
close all

syms rho

digits(1000)

% This code is associated with the paper
%       "Unconditionally stable space-time isogeometric discretization
%           for the wave equation in Hamiltonian formulation"
% by M. Ferrari, S. Fraschini, G. Loli and I. Perugia


% In this code we verify (for p = 2,...,17) that there are no rho>0 such 
% that two of the four blocks outlined in Remark 3.5 of the matrices
%               rho*B_h^p*C_h^p*(B_h^p)^(-1)*C_h^p
% are simultaneously singular (B_h^p and C_h^p defined in Equation (7))

% This code accompanies the proof of Property 3.16.

% B_h^p and C_h^p have been previously computed with 1000 digits precision


for p = 19:19 
    Nt = 2^7;

    cd ./C_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_C_sym'))
    cd ..
    cd ./B_matrices
    load(strcat(num2str(p),'_',num2str(Nt),'_B_sym'))
    cd ..

    B = B_sym/Nt;
    C = C_sym;

    X1 = B*C/B*C;
    X2 = B^2;

    m = 2*p+2;
    ell = 2*p-2;

    %if (A*rho + B)x = 0 for some x not zero
    %then rho is an eigenvalue of -inv(A)*B

    %first set of critical rhos
    F11 = X1(1:m,ell+1:ell+m);
    F12 = X2(1:m,ell+1:ell+m);
    set_1 = -eig(F11\F12);

    %second set of critical rhos
    F21 = X1(m+1:m+ell,1:ell);
    F22 = X2(m+1:m+ell,1:ell);
    set_2 = -eig(F21\F22);

    %then we compute the persymmetrics
    X1J = exc(size(X1,1))*X1'*exc(size(X1,1));
    X2J = exc(size(X1,1))*X2'*exc(size(X1,1));

    %third set of critical rhos
    F11J = X1J(1:m,ell+1:ell+m);
    F12J = X2J(1:m,ell+1:ell+m);
    set_3 = -eig(F11J\F12J);

    %fouth set of critical rhos
    F21J = X1J(m+1:m+ell,1:ell);
    F22J = X2J(m+1:m+ell,1:ell);
    set_4 = -eig(F21J\F22J);

    bad_rho_set_1 = set_1(abs(imag(set_1))<1e-16 & real(set_1)>0);
    bad_rho_set_2 = set_2(abs(imag(set_2))<1e-16 & real(set_2)>0);
    bad_rho_set_3 = set_3(abs(imag(set_3))<1e-16 & real(set_3)>0);
    bad_rho_set_4 = set_4(abs(imag(set_4))<1e-16 & real(set_4)>0);

    if ~isempty(bad_rho_set_1)
        bad_rho_set_1   
        double_bad_rho_set_1 = double(bad_rho_set_1)
    end
    if ~isempty(bad_rho_set_2)
        bad_rho_set_2 
        double_bad_rho_set_2 = double(bad_rho_set_2)
    end
    if ~isempty(bad_rho_set_3)
        bad_rho_set_3
        double_bad_rho_set_3 = double(bad_rho_set_3)
    end
    if ~isempty(bad_rho_set_4)
        bad_rho_set_4   
        double_bad_rho_set_4 = double(bad_rho_set_4)
    end
end



