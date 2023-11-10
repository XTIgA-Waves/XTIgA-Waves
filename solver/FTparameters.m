function [U,Ds,Q,Z,Tr,A] = FTparameters(Ks,Kt,Ms,Wt)

%% Space diagonalization
if~iscell(Ks)
    dim=1;
    Ks={Ks};
    Ms={Ms};
else
    dim=numel(Ks);
end

for idim=1:dim
    [u,Ds{idim}] = eig(full(Ks{idim}+Ks{idim}')/2,full(Ms{idim}+Ms{idim}')/2,'vector');
    for j =1:size(u,2)
        v = u(:,j);
        U{idim}(:,j) = v/sqrt(v'*(Ms{idim}*v));
    end
end
if dim==2
    Ds=kron(Ds{2},ones(size(Ds{1})))+kron(ones(size(Ds{2})),Ds{1});
elseif dim==1
    Ds=Ds{1};
end

%% Time triangulation
[A,B,Q,Z]=qz(full(Wt),full(Kt),'complex');
Tr=A\B;

end