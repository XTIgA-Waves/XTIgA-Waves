function y = ftApply(x,Q,Z,Tr,A,Ks,Ms)

nT = size(Tr,1);
nS = length(Ks);

Y = reshape(x,nS,nT);
Y = Y*Q.';
Y = backwardSub(A,Y.').';
Y = blockBackSub(Y,Tr,Ks,Ms);
Y = Y*Z.';
y = real(Y(:));

end

function y = blockBackSub(x,Tr,Ks,Ms)

nT = size(Tr,1);
nS = size(Ms,1);
y = zeros(nS,nT);

y(:,end) = (Tr(end,end)*Ks-Ms)\x(:,end);
for iR = (nT-1):-1:1
    tmp = Ks*y(:,(iR+1):end)*Tr(iR,(iR+1):end).';
    y(:,iR) = (Tr(iR,iR)*Ks-Ms)\(x(:,iR)-tmp);
end

end