function Y = tsolve(A,X,d)

% This function multiplies by inv(A) the tensor X
% in the direction d = 1,2,3. 

[n1, n2, n3] = size(X);

if d == 1 
    Y = A\reshape(X,n1,n2*n3);
    Y = reshape(Y,n1,n2,n3);  
elseif d == 2
    Y = A\reshape(permute(X,[2 1 3]),n2,n1*n3);
    Y = permute(reshape(Y,n2,n1,n3),[2 1 3]);  
elseif d == 3 
    Y = A\reshape(permute(X,[3 2 1]),n3,n1*n2);
    Y = permute(reshape(Y,n3,n2,n1),[3 2 1]);   
end

end


