function y = FTapplication(x,U,Ds,Q,Z,Tr,A)

U=U{1};

nt = size(Tr,1);
ns = length(Ds);

Y = reshape(x,ns,nt);
Y = U'*Y*Q.';
Y = backwardSub(A,Y.').';
tmp = blocBackwardSub(Y(:),Tr,Ds);
Y = reshape(tmp,ns,nt);
Y = U*Y*Z.';
y = real(Y(:));

end