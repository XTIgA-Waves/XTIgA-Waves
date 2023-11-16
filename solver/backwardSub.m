function[x]=backwardSub(U,B)

[r,c]=size(B);
x=zeros(r,c);

x(end,:)=B(end,:)/U(end,end);
for ir=r-1:-1:1
    x(ir,:)=(B(ir,:)-U(ir,ir+1:end)*x(ir+1:end,:))/U(ir,ir);
end

end