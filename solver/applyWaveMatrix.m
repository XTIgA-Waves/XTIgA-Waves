function Y = applyWaveMatrix(x,Kt,Wt,Ks,Ms,Mr,Wr)

if nargin == 5
    Y=reshape(x,[size(Ms,1),size(Wt,1)]);
    Y=Ks*Y*Wt'-Ms*Y*Kt';
    Y=Y(:);
else
    Y=reshape(x,size(Ms,1),size(Wt,1));
    Y=Ks*Y*Wt'-Ms*Y*Kt'+Mr*Y*Wr';
    Y=Y(:);
end

end