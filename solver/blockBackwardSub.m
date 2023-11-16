function y = blocBackwardSub(x,Tr,Ds)

nt = size(Tr,1);
ns = length(Ds);

y = zeros(nt*ns,1);

y(end-ns+1:end) = x(end-ns+1:end)./(Ds-Tr(end,end));

for ir=nt-1:-1:1
    tmp=zeros(ns,1);
    for iir=(ir+1):nt
        tmp=tmp-Tr(ir,iir).*y((iir-1)*ns+1:iir*ns);
    end
    y((ir-1)*ns+1:ir*ns) = (x((ir-1)*ns+1:ir*ns)-tmp)./(Ds-Tr(ir,ir));
end

end