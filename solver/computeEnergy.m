function E = computeEnergy(qn,qw,u,c,spaceST,geometryST,map_der,tPts)

E = zeros(numel(tPts),1);
dimS = numel(qn)-1;
pts = [cellfun(@(qn)qn(:),qn(1:end-1),'UniformOutput',false),linspace(0,1,2)];
sizePts = cellfun(@(p) numel(p),pts);
weightsS = cellfun(@(qw)qw(:),qw(1:end-1),'UniformOutput',false);
wS = 1;
for idim = dimS:-1:1
    wS = kron(wS,weightsS{idim});
end
grad_u_eval = sp_eval(u,spaceST,geometryST,pts,{'gradient'});
grad_u_eval = reshape(grad_u_eval,[dimS+1,prod(sizePts(1:end-1)),2]);
ndPts = cell(dimS,1);
[ndPts{:}] = ndgrid(pts{1:end-1});
ndPts = cellfun(@(pt)pt(:),ndPts,'UniformOutput',false);
for it = 1:2
    c_eval = c(ndPts{:},repmat(pts{end}(it),[prod(sizePts(1:end-1)),1]));
    Jac = map_der([pts(1:end-1),{0}]);
    JacS = squeeze(Jac(1,1,:));
    grad_x_l2_2 = sum(reshape(grad_u_eval(1:end-1,:,it),[],1).^2.*c_eval(:,it).^2.*wS.*abs(JacS))/2;
    grad_t_l2_2 = sum(reshape(grad_u_eval(end,:,it),[],1).^2.*wS.*abs(JacS))/2;
    E(it) = grad_x_l2_2 + grad_t_l2_2;
end

end