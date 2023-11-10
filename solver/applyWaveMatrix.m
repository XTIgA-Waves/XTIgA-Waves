function Y = applyWaveMatrix(x,Ks,Kt,Ms,Mt)

Y=reshape(x,size(Ms,1),size(Mt,1));
Y=Ks*Y*Mt'-Ms*Y*Kt';
Y=Y(:);

end