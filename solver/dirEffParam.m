function [Ut,Tr] = dirEffParam(invWtKt)

%% Time triangulation
[Ut,Tr] = schur(full(invWtKt),'complex');

end