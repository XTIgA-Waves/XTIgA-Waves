function solution = solveWaveSTstab(problemData,methodData,varargin)

%% Unwrapping input structures
cellfun(@(x) assignin('caller',x,problemData.(x)), fieldnames(problemData));
clear problemData;
cellfun(@(x) assignin('caller',x,methodData.(x)), fieldnames(methodData));
clear methodData;

%% Chosing linear solver
p = inputParser;
errorMsg = 'In presence of Robin boundary condition, the ''efficient'' version of the linear solver is not exact.';
validationFcn = @(x) assert((isempty(rbnSides)||~strcmp(x,'dirEff')),errorMsg);
addParameter(p,'Solver','dir',validationFcn);
addParameter(p,'ExactSol',false,@islogical);
addParameter(p,'Error',false,@islogical);
addParameter(p,'ErrorFine',false,@islogical);
addParameter(p,'Energy',false,@islogical);
p.parse(varargin{:});
cellfun(@(x) assignin('caller',x,p.Results.(x)), fieldnames(p.Results));
clear p varargin;

%% Spatial and Time structures
if ~exist('geometryS','var')
    geometryS = geo_load(geoName);
end
dimS = geometryS.rdim;
nsubS = ceil((nsubS./(geometryS.nurbs.number-geometryS.nurbs.order+1))-1); %#ok<*NODEF>
[rknotsS,zetaS] = kntrefine(geometryS.nurbs.knots,nsubS,degreeS,regularityS);

geometryT = geo_load(nrbline([0,0],[T,0]));
nsubT = ceil((nsubT./(geometryT.nurbs.number-geometryT.nurbs.order+1))-1);
[rknotsT,zetaT] = kntrefine(geometryT.nurbs.knots,nsubT,degreeT,regularityT);
clear nsubS nsubT;

%% Check for periodic conditions, and consistency with other boundary conditions
if ~isempty(perDir)
    tmp = kntunclamp ([rknotsS(:)',{rknotsT}], [degreeS,degreeT], [regularityS,regularityT], perDir);
    clear rknotsS;
    [rknotsS{1:dimS}] = tmp{1:dimS};
    rknotsT = tmp{end};
    clear tmp;
else
    perDir = [];
end

nquad = degreeS+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes(zetaS,quadRule);
mshS = msh_cartesian(zetaS,qn,qw,geometryS);
solution.hs = max(msh_precompute(mshS).element_size);
spaceS = sp_bspline(rknotsS,degreeS,mshS,[],perDir);

nquad = degreeT+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes(zetaT,quadRule);
mshT = msh_cartesian(zetaT,qn,qw,geometryT);
clear geometryT nquad quadRule qn qw;
solution.ht = max(msh_precompute(mshT).element_size);
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
Mt = op_u_v_tp(spaceT,spaceT,mshT);
Kt = op_gradu_gradv_tp(spaceT,spaceT,mshT);

if ~isempty(rbnSides)
    Wt = op_vel_dot_gradu_v_tp(spaceT,spaceT,mshT,@(t) ones(size(t)));
    nnzMr = sum(arrayfun(@(ibnd) numel(spaceS.boundary(ibnd).dofs),rbnSides))*(2*max(degreeS)+1);
    Mr = spalloc(size(Ms,1),size(Ms,2),nnzMr);
    for iside = rbnSides
        dofs = spaceS.boundary(iside).dofs;
        if mshS.ndim == 1
            Mr(dofs,dofs) = theta*c(mshS.map(iside-1),0);
        else
            mshSide = msh_eval_boundary_side (mshS,iside);
            spSide = sp_precompute (spaceS.boundary(iside), mshSide,'value',true,'gradient',false);
            x = cell(mshS.rdim,1);
            for idim = 1:mshS.rdim
                x{idim} = reshape (mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
            end
            Mr(dofs,dofs) = Mr(dofs,dofs) + op_u_v(spSide,spSide,mshSide, theta*c(x{1},x{2},zeros(size(x{1}))));
        end
    end
end
clear mshS mshSide spSide nnzMr dofs;

if ~isempty(stab_param)
    if stab_param>0
        shapeFunctionsDerT = sp_bspline_1d_param_all_derivatives (rknotsT,degreeT,mshT.qn{1});
        sizeReshapeT = size(shapeFunctionsDerT);
        sizeReshapeT = sizeReshapeT(2:end);
        shapeFunctionsDerT = reshape(shapeFunctionsDerT(degreeT+1,:,:,:),sizeReshapeT);
        clear sizeReshapeT;
        Pt = stab_param*solution.ht^(2*degreeT)*op_mat_hostab_tp(spaceT, spaceT, mshT, shapeFunctionsDerT,shapeFunctionsDerT,degreeT);
        clear shapeFunctionsDerT stab_param;
    end
end
clear mshT spaceT;

%% Assembling Global structures
solution.geometryST = geo_load_st(geometryS,0,T);
nquad = [degreeS,degreeT]+1;
quadRule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes([zetaS(:)',{zetaT}],quadRule);
solution.mshST = msh_cartesian([zetaS(:)',{zetaT}],qn,qw,solution.geometryST);
solution.spaceST = sp_bspline([rknotsS(:)',{rknotsT}],[degreeS,degreeT],solution.mshST,[],perDir);
clear geometryS quadRule qn qw zetaS zetaT rknotsT;
rhs = op_f_v_tp(solution.spaceST,solution.mshST,f);
solution.u = zeros (solution.spaceST.ndof,1);

%% Apply Dirichlet boundary conditions and initial condition (u(0)=u_0)
switch dimS
    case 1
        if ExactSol
            gDrchlt = @(x,t,iside) u_ex(x,t);
        else
            gDrchlt = @(x,t,iside) InitDrchlt1d(gDrchlt,gInit,x,t,dimS,iside);
        end
    case 2
        if ExactSol
            gDrchlt = @(x,y,t,iside) u_ex(x,y,t);
        else
            gDrchlt = @(x,y,t,iside) InitDrchlt2d(gDrchlt,gInit,x,y,t,dimS,iside);
        end
    case 3
        if ExactSol
            gDrchlt = @(x,y,z,t,iside) u_ex(x,y,z,t);
        else
            gDrchlt = @(x,y,z,t,iside) InitDrchlt3d(gDrchlt,gInit,x,y,z,t,dimS,iside);
        end
end
[uForced, trialForcedDofs] = sp_drchlt_l2_proj(solution.spaceST,solution.mshST,gDrchlt,union(drchltSides,2*dimS+1));
solution.u(trialForcedDofs) = uForced;
clear uForced;

testForcedDofs = [];
nent = 0;
for iside = union(drchltSides,numel(solution.spaceST.boundary))
    nent = nent + solution.mshST.boundary(iside).nel * solution.spaceST.boundary(iside).nsh_max^2;
    testForcedDofs = union (testForcedDofs, solution.spaceST.boundary(iside).dofs);
end

trialIntDofs = setdiff (1:solution.spaceST.ndof, trialForcedDofs);
testIntDofs = setdiff (1:solution.spaceST.ndof, testForcedDofs);

drchltDofsST = [];
drchltDofsS = [];
for iside = drchltSides
    drchltDofsST = union(drchltDofsST,solution.spaceST.boundary(iside).dofs);
    drchltDofsS = union(drchltDofsS,spaceS.boundary(iside).dofs);
end
drchltDofsST = setdiff(drchltDofsST,solution.spaceST.boundary(end-1).dofs);
intDofsS = setdiff(1:spaceS.ndof,drchltDofsS);
solution.nDofS = numel(intDofsS);
clear nent spaceS;

tmp0 = Ms*solution.u(solution.spaceST.boundary(end-1).dofs)*Kt(:,1)'...
    -Ks*solution.u(solution.spaceST.boundary(end-1).dofs)*Mt(:,1)'...
    +Ks*solution.u(solution.spaceST.boundary(end-1).dofs)*Pt(:,1)';

tmpDrchlt = Ms(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Kt(:,2:end)'...
    -Ks(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Mt(:,2:end)'...
    +Ks(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Pt(:,2:end)';

if ~isempty(rbnSides)
    tmp0 = tmp0 -Mr*solution.u(solution.spaceST.boundary(end-1).dofs)*Wt(:,1)';
    tmpDrchlt = tmpDrchlt -Mr(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Wt(:,2:end)';
end

rhs(testIntDofs) = rhs(testIntDofs) + tmp0(testIntDofs(:)) + tmpDrchlt(testIntDofs(:));
clear trialForcedDofs testForcedDofs tmp0 tmpDrchlt drchltDofsS drchltDofsST;

%% Apply Initial conditions (u'(0)=u_1)
dofs = solution.spaceST.boundary(2*dimS+1).dofs;
mshSide = msh_eval_boundary_side (solution.mshST,2*dimS+1);
spSide = sp_eval_boundary_side (solution.spaceST,mshSide);
x = cell(solution.mshST.rdim,1);
for idim = 1:solution.mshST.rdim
    x{idim} = reshape (mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
end
rhs(dofs) = rhs(dofs) + op_f_v(spSide,mshSide,gInitDer(x{:}));
clear mshSide spSide x dofs;

%% Apply Neumann boundary conditions
for iside = nmnnSides
    dofs = solution.spaceST.boundary(iside).dofs;
    mshSide = msh_eval_boundary_side (solution.mshST,iside);
    spSide = sp_eval_boundary_side (solution.spaceST,mshSide);
    x = cell(solution.mshST.rdim,1);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape (mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    if ExactSol
        rhs(dofs) = rhs(dofs) + op_fdotn_v(spSide,mshSide,grad_u_ex_c2(x{:}));
    else
        rhs(dofs) = rhs(dofs) + op_f_v(spSide,mshSide,gNmnn(x{:}));
    end
end
clear spSide mshSide x dofs;

%% Apply Robin boundary conditions
for iside = rbnSides
    dofs = solution.spaceST.boundary(iside).dofs;
    mshSide = msh_eval_boundary_side (solution.mshST,iside);
    spSide = sp_eval_boundary_side (solution.spaceST,mshSide);
    x = cell(solution.mshST.rdim,1);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    if ExactSol
        rhs(dofs) = rhs(dofs) + op_fdotn_v(spSide,mshSide,grad_u_ex_c2(x{:}))...
            + op_f_v(spSide,mshSide,theta*c_dt_u_ex(x{:}));
    else
        rhs(dofs) = rhs(dofs) + op_f_v(spSide,mshSide,gRbn(x{:}));
    end
end
clear mshSide spSide x dofs;

%% Solve the linear system
% case 'dir':    direct solver for the linear system
% case 'dirEff': a faster direct solver based on a triangularization of the
%                time matrices. In presence of Robin boundary condition,
%                this solver cannot be used.
% case 'precT':   GMRES solver preconditioned by dirEff.
switch Solver
    case 'dir'
        Mat = kron(Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Ks(intDofsS,intDofsS))-kron(Kt(1:end-1,2:end),Ms(intDofsS,intDofsS));
        if ~isempty(rbnSides)
            Mat = Mat+kron(Wt(1:end-1,2:end),Mr(intDofsS,intDofsS));
        end
        solution.u(trialIntDofs) = Mat\rhs(testIntDofs);
    case 'dirEff'
        [U,Ds,Q,Z,Tr,A] = FTparameters(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        clear Ks Kt Ms Mt intDofsS;
        solution.u(trialIntDofs) = FTapplication(rhs(testIntDofs),U,Ds,Q,Z,Tr,A);
        clear U Ds Q Z Tr A;
    case 'precT'
        [U,Ds,Q,Z,Tr,A] = FTparameters(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        Prec=@(x) FTapplication(x,U,Ds,Q,Z,Tr,A);
        if exist('Mr','var')
            Mat = @(x) applyWaveMatrix(x,Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Mr(intDofsS,intDofsS),Wt(1:end-1,2:end));
        else
            Mat = @(x) applyWaveMatrix(x,Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        end
        if~exist('tol','var')
            tol=10^-12;
        end
        [y,~,relres,iter] = gmres(Mat,rhs(testIntDofs),50,tol,numel(trialIntDofs),Prec);
        solution.u(trialIntDofs) = y;
        fprintf(1,'\n*Solver info: iter=[%d,%d]; Relres=%g*\n',iter(1),iter(2),relres);
        clear testIntDofs U Ds Q Z Tr A Ks Kt Ms Mt Mr Wt intDofsS iter relres;
end
solution.nDof = numel(trialIntDofs);
clear trialIntDofs rhs;

%% Errori in norma L2 e in norma H1
if Error
    [errH1,errL2,errH1s] = sp_h1_c_error_tp(solution.spaceST,solution.mshST,solution.u,u_ex,grad_u_ex,c);
    [normH1,normL2,normH1s] = sp_h1_c_error_tp(solution.spaceST,solution.mshST,zeros(size(solution.u)),u_ex, grad_u_ex,c);
    solution.relErrH1 = errH1/normH1;
    solution.relErrL2 = errL2/normL2;
    solution.relErrH1s = errH1s/normH1s;
end

if ErrorFine
    solFine=load(solFine).solution;
    [errH1, errL2, errH1s] = sp_h1_error_spline_ex(solution.spaceST,solution.mshST,solution.u,solFine.spaceST,solution.geometryST,solFine.u,c);
    [normH1,normL2,normH1s] = sp_h1_error_spline_ex(solution.spaceST,solution.mshST,solution.u*0,solFine.spaceST,solution.geometryST,solFine.u,c);
    solution.relErrH1=errH1/normH1;
    solution.relErrL2=errL2/normL2;
    solution.relErrH1s=errH1s/normH1s;
end

%% Compute Energy at t=0 and t=T
if Energy
    solution.E0Et = computeEnergy(solution.mshST.qn,solution.mshST.qw,solution.u,c,solution.spaceST,solution.geometryST,solution.mshST.map_der,[0,1]);
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