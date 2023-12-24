function geoST = geo_load_st(geoS,t0,t1)

map{1} = @(pts) evalGeoST(pts,geoS,t0,t1);
map{2} = @(pts) jacobianGeoST(pts,geoS,t0,t1);
map{3} = @(pts) hessianGeoST(pts,geoS);

geoST = geo_load(map);

end

function val = evalGeoST(pts,geoS,t0,t1)
if iscell(pts)
    t = reshape(t0+(t1-t0)*pts{end},1,[]);
    phyPts = geoS.map(pts(1:end-1));
    phyPts = repmat(phyPts,[1,numel(t)]);
    t = repmat(t,[1,prod(cellfun(@numel,pts(1:end-1)))]);
    val = cat(1,phyPts,t);
else
    phyPts = geoS.map(pts(1:end-1,:));
    t = t0+(t1-t0)*pts(end,:);
    val = cat(1,phyPts,t);
end
end

function Jac = jacobianGeoST(pts,geoS,t0,t1)
if iscell(pts)
    Jac = repmat(geoS.map_der(pts(1:end-1)),[1,1,numel(pts{end})]);
    Jac = cat(1,Jac,zeros([1,size(Jac,2),size(Jac,3)]));
    Jac = cat(2,Jac,repmat([zeros(size(Jac,1)-1,1);t1-t0],[1,1,size(Jac,3)]));
    Jac = reshape(Jac,[size(Jac,1),size(Jac,2),cellfun(@numel,pts)]);
else
    Jac = geoS.map_der(pts(1:end-1,:));
    Jac = cat(1,Jac,zeros(1,size(Jac,2:3)));
    Jac = cat(2,Jac,repmat([zeros(size(Jac,1)-1,1);t1-t0],[1,1,size(Jac,3)]));
    Jac = reshape(Jac,[size(Jac,1:2),size(pts,2:ndims(pts))]);
end
end

function Hes = hessianGeoST(pts,geoS)
if iscell(pts)
    Hes = repmat(geoS.map_der2(pts(1:end-1)),[1,1,1,numel(pts{end})]);
    Hes = cat(2,Hes,zeros([size(Hes,1),1,size(Hes,3:4)]));
    Hes = cat(3,Hes,zeros([size(Hes,1:2),1,size(Hes,4)]));
    Hes = cat(1,Hes,zeros([1,size(Hes,2:ndims(Hes))]));
    Hes = reshape(Hes,[size(Hes,1:3),cellfun(@numel,pts)]);
else
    Hes = geoS.map_der2(pts(1:end-1,:));
    Hes = cat(2,Hes,zeros([size(Hes,1),1,size(Hes,3:4)]));
    Hes = cat(3,Hes,zeros([size(Hes,1:2),1,size(Hes,4)]));
    Hes = cat(1,Hes,zeros([1,size(Hes,2:ndims(Hes))]));
end
end