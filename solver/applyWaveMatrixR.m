function Y = applyWaveMatrixR(x,Ks,Kt,Ms,Mt,Mr,Wt)

Y=reshape(x,size(Ms,1),size(Mt,1));
Y=Ks*Y*Mt'-Ms*Y*Kt'+Mr*Y*Wt';
Y=Y(:);

end