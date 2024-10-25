function y = dirEffApply(x,Ms,Ks,Ut,Tr,Kt,Mr)

nT = size(Ut,1);

Y = reshape(x,[],nT);
Y = Ut'*(Kt\Y.');
if nargin == 7
    Y = blockBackSub(Y.',Tr,inv(Tr),Ms,Ks,Mr);
else
    Y = blockBackSub(Y.',Tr,inv(Tr),Ms,Ks);
end
Y = Y*Ut.';
y = real(Y(:));

end

function y = blockBackSub(x,Tr,invTr,Ms,Ks,Mr)

nT = size(Tr,1);
nS = size(Ms,1);
y = zeros(nS,nT);

if nargin == 6
    y(:,end) = (invTr(end,end)*Ks+Tr(end,end)*Ms+Mr)\x(:,end);
    for iR = nT-1:-1:1
        tmp = sum(Ks*y(:,iR+1:end)*invTr(iR,iR+1:end).'+Ms*y(:,iR+1:end)*Tr(iR,iR+1:end).',2);
        y(:,iR) = (invTr(iR,iR)*Ks+Tr(iR,iR)*Ms+Mr)\(x(:,iR)-tmp);
    end
else
    y(:,end) = (invTr(end,end)*Ks+Tr(end,end)*Ms)\x(:,end);
    for iR = nT-1:-1:1
        tmp = sum(Ks*y(:,iR+1:end)*invTr(iR,iR+1:end).'+Ms*y(:,iR+1:end)*Tr(iR,iR+1:end).',2);
        y(:,iR) = (invTr(iR,iR)*Ks+Tr(iR,iR)*Ms)\(x(:,iR)-tmp);
    end
end


end