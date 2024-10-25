function sol = solveWaveSTFO(problemData,methodData,varargin)

%% Unwrapping input structures
cellfun(@(x) assignin('caller',x,problemData.(x)),fieldnames(problemData));
clear problemData;
cellfun(@(x) assignin('caller',x,methodData.(x)),fieldnames(methodData));
clear methodData;

%% Chosing linear solver
p = inputParser;
addParameter(p,'Solver','dir');
addParameter(p,'ExactSol',false,@islogical);
addParameter(p,'Error',false,@islogical);
addParameter(p,'ErrorFine',false,@islogical);
addParameter(p,'Energy',false,@islogical);
p.parse(varargin{:});
cellfun(@(x) assignin('caller',x,p.Results.(x)),fieldnames(p.Results));
clear p varargin;

%% Spatial and Time structures
if~exist('geometryS','var')
    geometryS = geo_load(geoName);
end
dimS = geometryS.rdim;
nsubS = ceil((nsubS./(geometryS.nurbs.number-geometryS.nurbs.order+1))-1);
[rknotsS,zetaS] = kntrefine(geometryS.nurbs.knots,nsubS,degreeS,regularityS);

geometryT = geo_load(nrbline([0,0],[T,0]));
nsubT = ceil((nsubT./(geometryT.nurbs.number-geometryT.nurbs.order+1))-1);
[rknotsT,zetaT] = kntrefine(geometryT.nurbs.knots,nsubT,degreeT,regularityT);
clear nsubS nsubT;

%% Check for periodic conditions, and consistency with other boundary conditions
if ~isempty(perDir)
    tmp0 = kntunclamp([rknotsS(:)',{rknotsT}], [degreeS,degreeT], [regularityS,regularityT], perDir);
    clear rknotsS;
    [rknotsS{1:dimS}] = tmp0{1:dimS};
    rknotsT = tmp0{end};
    clear tmp0;
else
    perDir = [];
end

nquad = degreeS+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes(zetaS,quadRule);
mshS = msh_cartesian(zetaS,qn,qw,geometryS);
sol.hs = max(msh_precompute(mshS).element_size);
spaceS = sp_bspline(rknotsS,degreeS,mshS,[],perDir);

nquad = degreeT+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes(zetaT,quadRule);
mshT = msh_cartesian(zetaT,qn,qw,geometryT);
clear geometryT qn qw nquad;
sol.ht = max(msh_precompute(mshT).element_size);
spaceT = sp_bspline(rknotsT,degreeT,mshT);

%% Assembling Spatial matrices
switch dimS
    case 1
        Ks = op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x) c(x,0).^(2));
    case 2
        Ks = op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x,y) c(x,y,0).^(2));
    case 3
        Ks = op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x,y,z) c(x,y,z,0).^(2));
end
Ms = op_u_v_tp(spaceS,spaceS,mshS);

%% Assembling Time matrices
Kt = op_gradu_gradv_tp(spaceT,spaceT,mshT);
Wt = op_vel_dot_gradu_v_tp(spaceT,spaceT,mshT,@(t) ones(size(t)));

if ~isempty(rbnSides)
    nnzMr = sum(arrayfun(@(ibnd) numel(spaceS.boundary(ibnd).dofs),rbnSides))*(2*max(degreeS)+1);
    Mr = spalloc(size(Ms,1),size(Ms,2),nnzMr);
    for iside = rbnSides
        dofs = spaceS.boundary(iside).dofs;
        if mshS.ndim == 1
            Mr(dofs,dofs) = theta*c(mshS.map(iside-1),0);
        else
            mshSide = msh_eval_boundary_side(mshS,iside);
            spSide = sp_precompute(spaceS.boundary(iside), mshSide,'value',true,'gradient',false);
            x = cell(mshS.rdim,1);
            for idim = 1:mshS.rdim
                x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
            end
            Mr(dofs,dofs) = Mr(dofs,dofs) + op_u_v(spSide,spSide,mshSide,theta*c(x{:},zeros(size(x{1}))));
        end
    end
end
qnS = mshS.qn;
clear mshS mshSide spSide nnzMr x mshT spaceT;

%% Assembling Global structures
sol.geometryST = geo_load_st(geometryS,0,T);
nquad = [degreeS,degreeT]+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes([zetaS(:)',{zetaT}],quadRule);
sol.mshST = msh_cartesian([zetaS(:)',{zetaT}],qn,qw,sol.geometryST);
sol.spaceST = sp_bspline([rknotsS(:)',{rknotsT}],[degreeS,degreeT],sol.mshST,[],perDir);
clear quadRule qn qw zetaS zetaT rknotsT;
rhs2 = op_f_dtv_tp(sol.spaceST,sol.mshST,f);
sol.u = zeros(sol.spaceST.ndof,1);
sol.v = zeros(sol.spaceST.ndof,1);

%% Apply Dirichlet boundary conditions and initial conditions: u(0)=u_0, v(0)=u_1
switch dimS
    case 1
        if ExactSol
            gDrchlt = @(x,t,iside) u_ex(x,t);
            dt_gDrchlt = @(x,t,iside) dt_u_ex(x,t);
        else
            gDrchlt = @(x,t,iside) InitDrchlt1d(gDrchlt,gInit,x,t,dimS,iside);
            dt_gDrchlt = @(x,t,iside) InitDrchlt1d(gDrchltDer,gInitDer,x,t,dimS,iside);
        end
    case 2
        if ExactSol
            gDrchlt = @(x,y,t,iside) u_ex(x,y,t);
            dt_gDrchlt = @(x,y,t,iside) dt_u_ex(x,y,t);
        else
            gDrchlt = @(x,y,t,iside) InitDrchlt2d(gDrchlt,gInit,x,y,t,dimS,iside);
            dt_gDrchlt = @(x,y,t,iside) InitDrchlt2d(gDrchltDer,gInitDer,x,y,t,dimS,iside);
        end
    case 3
        if ExactSol
            gDrchlt = @(x,y,z,t,iside) u_ex(x,y,z,t);
            dt_gDrchlt = @(x,y,z,t,iside) dt_u_ex(x,y,z,t);
        else
            gDrchlt = @(x,y,z,t,iside) InitDrchlt3d(gDrchlt,gInit,x,y,z,t,dimS,iside);
            dt_gDrchlt = @(x,y,z,t,iside) InitDrchlt3d(gDrchltDer,gInitDer,x,y,z,t,dimS,iside);
        end
end
[uForced,forcedDofsU] = sp_drchlt_l2_proj(sol.spaceST,sol.mshST,gDrchlt,union(drchltSides,2*dimS+1));
[vForced,forcedDofsV] = sp_drchlt_l2_proj(sol.spaceST,sol.mshST,dt_gDrchlt,union(drchltSides,2*dimS+1));
sol.u(forcedDofsU) = uForced;
sol.v(forcedDofsV) = vForced;
intDofs = setdiff(1:sol.spaceST.ndof,forcedDofsU);
% intDofsV = setdiff(1:sol.spaceST.ndof,forcedDofsV);

drchltDofsS = [];
drchltDofsST = [];
for iside = drchltSides
    drchltDofsS = union(drchltDofsS,spaceS.boundary(iside).dofs);
    drchltDofsST = union(drchltDofsST,sol.spaceST.boundary(iside).dofs);
end
intDofsS = setdiff(1:spaceS.ndof,drchltDofsS);
sol.nDofS = numel(intDofsS);
testIntDofs = setdiff(1:sol.spaceST.ndof,union(drchltDofsST,sol.spaceST.boundary(end).dofs));
drchltDofsST = setdiff(drchltDofsST,sol.spaceST.boundary(end-1).dofs);

tmp0 = - Ms(intDofsS,:)*sol.u(sol.spaceST.boundary(end-1).dofs)*Kt(1:end-1,1)'...
    + Ms(intDofsS,:)*sol.v(sol.spaceST.boundary(end-1).dofs)*Wt(1,1:end-1);
tmpDrchlt = - Ms(intDofsS,drchltDofsS)*reshape(sol.u(drchltDofsST),numel(drchltDofsS),size(Wt,2)-1)*Kt(2:end,1:end-1)...
    + Ms(intDofsS,drchltDofsS)*reshape(sol.v(drchltDofsST),numel(drchltDofsS),size(Wt,2)-1)*Wt(2:end,1:end-1);
rhs1 = tmp0(:) + tmpDrchlt(:);

tmp0 = - Ks(intDofsS,:)*sol.u(sol.spaceST.boundary(end-1).dofs)*Wt(1,1:end-1)...
    - Ms(intDofsS,:)*sol.v(sol.spaceST.boundary(end-1).dofs)*Kt(1:end-1,1)';
tmpDrchlt = - Ks(intDofsS,drchltDofsS)*reshape(sol.u(drchltDofsST),numel(drchltDofsS),size(Wt,2)-1)*Wt(2:end,1:end-1)...
    - Ms(intDofsS,drchltDofsS)*reshape(sol.v(drchltDofsST),numel(drchltDofsS),size(Wt,2)-1)*Kt(2:end,1:end-1);

if ~isempty(rbnSides)
    tmp0 = tmp0 - Mr(intDofsS,:)*sol.u(sol.spaceST.boundary(end-1).dofs)*Kt(1,1:end-1);
    tmpDrchlt = tmpDrchlt - Mr(intDofsS,drchltDofsS)*reshape(sol.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Kt(2:end,1:end-1);
end
rhs2(testIntDofs) = rhs2(testIntDofs) + tmp0(:) + tmpDrchlt(:);
clear uForced vForced forcedDofsU forcedDofsV tmp0 tmpDrchlt;

%% Apply Neumann boundary conditions
for iside = nmnnSides
    dofs = sol.spaceST.boundary(iside).dofs;
    mshSide = msh_eval_boundary_side(sol.mshST,iside);
    spSide = sp_eval_boundary_side_dt(sol.spaceST,mshSide);
    spSide.shape_functions = spSide.shape_function_gradients(end,:,:,:);
    x=cell(sol.mshST.rdim,1);
    for idim = 1:sol.mshST.rdim
        x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    if ExactSol
        rhs2(dofs) = rhs2(dofs) + op_fdotn_v(spSide,mshSide,grad_u_ex_c2(x{:})); % spSide contains dt_v
    else
        rhs2(dofs) = rhs2(dofs) + op_f_v(spSide,mshSide,gNmnn(x{:}));
    end
end
clear spSide mshSide x dofs;

%% Apply Robin boundary conditions
for iside = rbnSides
    dofs = sol.spaceST.boundary(iside).dofs;
    mshSide = msh_eval_boundary_side(sol.mshST,iside);
    spSide = sp_eval_boundary_side_dt(sol.spaceST,mshSide);
    spSide.shape_functions = spSide.shape_function_gradients(end,:,:,:);
    x=cell(sol.mshST.rdim,1);
    for idim = 1:sol.mshST.rdim
        x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    if ExactSol
        rhs2(dofs) = rhs2(dofs) + op_fdotn_v(spSide,mshSide,grad_u_ex_c2(x{:})) + op_f_v(spSide,mshSide,theta*c_dt_u_ex(x{:})); % spSide contains dt_v
    else
        rhs2(dofs) = rhs2(dofs) + op_f_v(spSide,mshSide,gRbn(x{:}));
    end
end
clear mshSide spSide x dofs;

%% Solve the linear system
% case 'dir':    direct solver for the linear system
switch Solver
    case 'dir'
        K21 = kron(Wt(2:end,1:end-1)',Ks(intDofsS,intDofsS));
        K22 = kron(Kt(1:end-1,2:end),Ms(intDofsS,intDofsS));
        if ~isempty(rbnSides)
            K21 = K21 + kron(Kt(1:end-1,2:end),Mr(intDofsS,intDofsS));
        end
        K22_invK12i_K11i = kron(Kt(1:end-1,2:end)*(Wt(2:end,1:end-1)'\Kt(1:end-1,2:end)),Ms(intDofsS,intDofsS));
        clear Ks;
        tmp1 = reshape(rhs1,[numel(intDofsS),size(Wt,1)-1]);
        tmp1 = tsolve(Ms(intDofsS,intDofsS),tmp1,1);
        tmp1 = tsolve(Wt(2:end,1:end-1)',tmp1,2);
        sol.u(intDofs) = (K21 + K22_invK12i_K11i)\(K22*tmp1(:) + rhs2(testIntDofs));
        clear rhs2;
        tmp1 = tmprod(reshape(sol.u(intDofs),[numel(intDofsS),size(Wt,1)-1]),Ms(intDofsS,intDofsS),1);
        tmp1 = tmprod(tmp1,Kt(1:end-1,2:end),2);
        y = reshape(tmp1(:) - rhs1,[numel(intDofsS),size(Wt,1)-1]);
        clear tmp1 rhs1;
        y = tsolve(Ms(intDofsS,intDofsS),y,1);
        y = tsolve(Wt(2:end,1:end-1)',y,2);
        sol.v(intDofs) = y(:);
        clear K21 K22 K22_invK12i_K11i;
    case 'dirEff'
        tmp1 = reshape(rhs1,[numel(intDofsS),size(Wt,1)-1]);
        tmp1 = tsolve(Ms(intDofsS,intDofsS),tmp1,1);
        tmp1 = tsolve(Wt(2:end,1:end-1)',tmp1,2);
        invWtKt = (Wt(2:end,1:end-1)'\Kt(1:end-1,2:end));
        [Ut,Tr] = dirEffParam(invWtKt);
        clear rknotsS;
        tmp1 = Ms(intDofsS,intDofsS)*tmp1*Kt(1:end-1,2:end)';
        rhs2 = tmp1(:)+rhs2(testIntDofs);
        if ~isempty(rbnSides)
            sol.u(intDofs) = dirEffApply(rhs2,Ms(intDofsS,intDofsS),Ks(intDofsS,intDofsS),Ut,Tr,Kt(1:end-1,2:end),Mr(intDofsS,intDofsS));
        else
            sol.u(intDofs) = dirEffApply(rhs2,Ms(intDofsS,intDofsS),Ks(intDofsS,intDofsS),Ut,Tr,Kt(1:end-1,2:end));
        end
        clear rhs2;
        tmp1 = tmprod(reshape(sol.u(intDofs),[numel(intDofsS),size(Wt,1)-1]),Ms(intDofsS,intDofsS),1);
        tmp1 = tmprod(tmp1,Kt(1:end-1,2:end),2);
        y = reshape(tmp1(:)-rhs1,[numel(intDofsS),size(Wt,1)-1]);
        clear tmp1 rhs1;
        y = tsolve(Ms(intDofsS,intDofsS),y,1);
        y = tsolve(Wt(2:end,1:end-1)',y,2);
        sol.v(intDofs) = y(:);
end
sol.nDof = numel(intDofs);%+numel(intDofsV);
clear Wt Kt Ms Ks Mr y;

%% Errori in norma L2 e in norma H1
if Error
    [errH1,errL2,errH1s] = sp_h1_c_error_tp(sol.spaceST, sol.mshST, sol.u, u_ex, grad_u_ex,c);
    [normH1,normL2,normH1s] = sp_h1_c_error_tp(sol.spaceST, sol.mshST, zeros(size(sol.u)), u_ex, grad_u_ex,c);
    sol.relErrH1U = errH1/normH1;
    sol.relErrL2U = errL2/normL2;
    sol.relErrH1sU = errH1s/normH1s;
    [errH1,errL2,errH1s] = sp_h1_c_error_tp(sol.spaceST,sol.mshST,sol.v,dt_u_ex,grad_dt_u_ex,c);
    [normH1,normL2,normH1s] = sp_h1_c_error_tp(sol.spaceST,sol.mshST,zeros(size(sol.v)),dt_u_ex,grad_dt_u_ex,c);
    sol.relErrH1V = errH1/normH1;
    sol.relErrL2V = errL2/normL2;
    sol.relErrH1sV = errH1s/normH1s;
    clear errH1s errH1 errL2 normH1s normH1 normL2;
end

if ErrorFine
    solFine=load(solFine).solution;
    [errH1, errL2, errH1s] = sp_h1_error_spline_ex(sol.spaceST,sol.mshST,sol.u,solFine.spaceST,sol.geometryST,solFine.u,c);
    [normH1,normL2,normH1s] = sp_h1_error_spline_ex(sol.spaceST,sol.mshST,sol.u*0,solFine.spaceST,solFine.geometryST,solFine.u,c);
    sol.relErrH1U = errH1/normH1;
    sol.relErrL2U = errL2/normL2;
    sol.relErrH1sU=errH1s/normH1s;
    % [errH1,errL2,errH1s] = sp_h1_error_spline_ex(sol.spaceST,sol.mshST,sol.v,solFine.spaceST,sol.geometryST,solFine.v,c);
    % [normH1,normL2,normH1s] = sp_h1_error_spline_ex(sol.spaceST,sol.mshST,sol.v*0,solFine.spaceST,sol.geometryST,solFine.v,c);
    % sol.relErrH1V = errH1/normH1;
    % sol.relErrL2V = errL2/normL2;
    % sol.relErrH1sV = errH1s/normH1s;
end

%% Compute Energy at t=0 and t=T
if Energy
    sol.E0Et = computeEnergyFO(sol.mshST.qn,sol.mshST.qw,sol.u,sol.v,c,sol.spaceST,sol.geometryST,sol.mshST.map_der,[0,1]);
end

end


function gDrchlt = InitDrchlt3d(gDrchlt,gInit,x,y,z,t,dimS,iside)
if iside == 2*dimS+1
    gDrchlt = gInit(x,y,z);
else
    gDrchlt = gDrchlt(x,y,z,t);
end
end

function gDrchlt = InitDrchlt2d(gDrchlt,gInit,x,y,t,dimS,iside)
if iside == 2*dimS+1
    gDrchlt = gInit(x,y);
else
    gDrchlt = gDrchlt(x,y,t);
end
end

function gDrchlt = InitDrchlt1d(gDrchlt,gInit,x,t,dimS,iside)
if iside == 2*dimS+1
    gDrchlt = gInit(x);
else
    gDrchlt = gDrchlt(x,t);
end
end