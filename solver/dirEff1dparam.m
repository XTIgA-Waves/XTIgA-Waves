function [Us,Ds,Ut,Tr] = dirEff1dparam(Ms,Ks,invWtKt)

%% Space diagonalization
Us = zeros(size(Ms));
[u,Ds] = eig(full(Ks+Ks')/2,full(Ms+Ms')/2,'vector');
for j =1:size(u,2)
    v = u(:,j);
    Us(:,j) = v/sqrt(v'*(Ms*v));
end

%% Time triangulation
[Ut,Tr] = schur(full(invWtKt),'complex');

end