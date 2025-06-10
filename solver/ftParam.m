function [Q,Z,Tr,A] = ftParam(Kt,Wt)

[A,B,Q,Z]=qz(full(Kt),full(Wt),'complex');
Tr=A\B;

end