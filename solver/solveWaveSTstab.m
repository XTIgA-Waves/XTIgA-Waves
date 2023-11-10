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
addParameter(p,'computeError',false,@islogical);
addParameter(p,'computeEnergy',false,@islogical);
p.parse(varargin{:});
cellfun(@(x) assignin('caller',x,p.Results.(x)), fieldnames(p.Results));
clear p varargin;

%% Spatial and Time structures
if~exist('geometryS','var')
    geometryS=geo_load(geoName);
end
clear geoName;
dimS=geometryS.rdim;
nsubS = ceil((nsubS./(geometryS.nurbs.number-geometryS.nurbs.order+1))-1); %#ok<*NODEF>
[rknotsS,zetaS] = kntrefine(geometryS.nurbs.knots,nsubS,degreeS,regularityS);
clear nsubS;

nurbsT = nrbline([0,0],[T,0]);
geometryT = geo_load(nurbsT) ;
nsubT = ceil((nsubT./(geometryT.nurbs.number-geometryT.nurbs.order+1))-1);
[rknotsT,zetaT] = kntrefine(geometryT.nurbs.knots,nsubT,degreeT,regularityT);
clear nsubT;

%% Check for periodic conditions, and consistency with other boundary conditions
if ~isempty(perDir)
    tmp = kntunclamp ([rknotsS(:)',{rknotsT}], [degreeS,degreeT], [regularityS,regularityT], perDir);
    clear rknotsS;
    [rknotsS{1:dimS}]=tmp{1:dimS};
    rknotsT=tmp{end};
    clear tmp;
else
    perDir = [];
end

nquad = degreeS+1;
quad_rule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes(zetaS,quad_rule);
mshS = msh_cartesian(zetaS,qn,qw,geometryS);
solution.hs = max(msh_precompute(mshS).element_size);
spaceS = sp_bspline(rknotsS,degreeS,mshS,[],perDir);

nquad = degreeT+1;
quad_rule = msh_gauss_nodes(nquad);
clear nquad;
[qn,qw] = msh_set_quad_nodes(zetaT,quad_rule);
mshT = msh_cartesian(zetaT,qn,qw,geometryT);
clear geometryT;
solution.ht = max(msh_precompute(mshT).element_size);
spaceT = sp_bspline(rknotsT,degreeT,mshT);

%% Assembling Spatial matrices
switch dimS
    case 1
        Ks=op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x) c(x,0).^(2));
    case 2
        Ks=op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x,y) c(x,y,0).^(2));
end
Ms=op_u_v_tp(spaceS,spaceS,mshS);

%% Assembling Time matrices
Mt=op_u_v_tp(spaceT,spaceT,mshT);
Kt=op_gradu_gradv_tp(spaceT,spaceT,mshT);

if~isempty(rbnSides)
    Wt=op_vel_dot_gradu_v_tp(spaceT,spaceT,mshT,@(t) ones(size(t)));
    nnzMr=sum(arrayfun(@(ibnd) numel(spaceS.boundary(ibnd).dofs),rbnSides))*(2*max(degreeS)+1);
    Mr=spalloc(size(Ms,1),size(Ms,2),nnzMr);
    for iside = rbnSides
        dofs = spaceS.boundary(iside).dofs;
        if mshS.ndim==1
            Mr(dofs,dofs)=theta*c(mshS.map(iside-1));
        else
            msh_side = msh_eval_boundary_side (mshS,iside);
            sp_side = sp_precompute (spaceS.boundary(iside), msh_side,'value',true,'gradient',false);
            x=cell(mshS.rdim,1);
            for idim = 1:mshS.rdim
                x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
            end
            Mr(dofs,dofs) = Mr(dofs,dofs) + op_u_v(sp_side,sp_side,msh_side, theta*c(x{1},x{2},zeros(size(x{1}))));
        end
    end
end
clear mshS msh_side sp_side nnzMr;

if ~isempty(stab_param)
    if stab_param>0
        shapeFunctionsDerT = sp_bspline_1d_param_all_derivatives (rknotsT,degreeT,mshT.qn{1});
        sizeReshapeT = size(shapeFunctionsDerT);
        sizeReshapeT = sizeReshapeT(2:end);
        shapeFunctionsDerT = reshape(shapeFunctionsDerT(degreeT+1,:,:,:),sizeReshapeT);
        clear sizeReshapeT;
        Pt = op_mat_hostab_tp(spaceT, spaceT, mshT, shapeFunctionsDerT,shapeFunctionsDerT,degreeT);
        clear shapeFunctionsDerT;
        Mt = Mt - stab_param*solution.ht^(2*degreeT)*Pt;
        clear Pt stab_param;
    end
end
clear mshT spaceT;

%% Assembling Global structures
solution.geometryST = geo_load(nrbextrude(geometryS.nurbs,[zeros(1,dimS),T]));
clear geometryS;
nquad = [degreeS,degreeT]+1;
quad_rule = msh_gauss_nodes(nquad);
[qn,qw] = msh_set_quad_nodes([zetaS(:)',{zetaT}],quad_rule);
clear quad_rule;
solution.mshST = msh_cartesian([zetaS(:)',{zetaT}],qn,qw,solution.geometryST);
clear qn qw zetaS zetaT;
solution.spaceST = sp_bspline([rknotsS(:)',{rknotsT}],[degreeS,degreeT],solution.mshST,[],perDir);
clear rknotsS rknotsT degreeS degreeT;
rhs = op_f_v_tp(solution.spaceST,solution.mshST,f);
solution.u = zeros (solution.spaceST.ndof,1);

%% Apply Dirichlet boundary conditions and initial condition (u(0)=u_0)
switch dimS
    case 1
        gDrchlt=@(x,t,iside) InitDrchlt1d(gDrchlt,gInit,x,t,dimS,iside);
    case 2
        gDrchlt=@(x,y,t,iside) InitDrchlt2d(gDrchlt,gInit,x,y,t,dimS,iside);
end
[u_drchlt, drchlt_dofs_1] = sp_drchlt_l2_proj (solution.spaceST, solution.mshST, gDrchlt, union(drchltSides,numel(solution.spaceST.boundary)-1));
solution.u(drchlt_dofs_1) = u_drchlt;
clear u_drchlt;

drchlt_dofs_2 = [];
nent = 0;
for iside = union(drchltSides,numel(solution.spaceST.boundary))
    nent = nent + solution.mshST.boundary(iside).nel * solution.spaceST.boundary(iside).nsh_max^2;
    drchlt_dofs_2 = union (drchlt_dofs_2, solution.spaceST.boundary(iside).dofs);
end
clear nent;

int_dofs_1 = setdiff (1:solution.spaceST.ndof, drchlt_dofs_1);
clear drchlt_dofs_1;
int_dofs_2 = setdiff (1:solution.spaceST.ndof, drchlt_dofs_2);
clear drchlt_dofs_2;
tmp=Ks*solution.u(solution.spaceST.boundary(end-1).dofs)*(Mt(:,1))'-Ms*solution.u(solution.spaceST.boundary(end-1).dofs)*(Kt(:,1))';
tmp=tmp(:);
rhs(int_dofs_2) = rhs(int_dofs_2)-tmp(int_dofs_2);
drchltDofsST=[];
drchltDofsS=[];
for iside=drchltSides
    drchltDofsST=union(drchltDofsST,solution.spaceST.boundary(iside).dofs);
    drchltDofsS=union(drchltDofsS,spaceS.boundary(iside).dofs);
end
drchltDofsST=setdiff(drchltDofsST,solution.spaceST.boundary(end-1).dofs);
intDofsS=setdiff(1:spaceS.ndof,drchltDofsS);
solution.nDofS=numel(intDofsS);
clear spaceS;
tmp=Ks(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Mt(:,2:end)'-Ms(:,drchltDofsS)*reshape(solution.u(drchltDofsST),numel(drchltDofsS),size(Kt,2)-1)*Kt(:,2:end)';
clear drchltDofsS;
tmp=tmp(:);
rhs(int_dofs_2) = rhs(int_dofs_2)-tmp(int_dofs_2);
clear tmp drchltDofsST;

%% Apply Initial conditions (u'(0)=u_1)
dofs = solution.spaceST.boundary(2*dimS+1).dofs;
msh_side = msh_eval_boundary_side (solution.mshST,2*dimS+1);
sp_side = sp_eval_boundary_side (solution.spaceST,msh_side);
x=cell(solution.mshST.rdim,1);
for idim = 1:solution.mshST.rdim
    x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
end
rhs(dofs) = rhs(dofs) + op_f_v (sp_side,msh_side,gInitDer(x{:}));
clear sp_side dofs;

%% Apply Neumann boundary conditions
for iside = nmnnSides
    dofs = solution.spaceST.boundary(iside).dofs;
    msh_side = msh_eval_boundary_side (solution.mshST,iside);
    sp_side = sp_eval_boundary_side (solution.spaceST,msh_side);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
    end
    rhs(dofs) = rhs(dofs) + op_f_v (sp_side,msh_side,gNmnn(x{:}));
end
clear msh_side sp_side x dofs;

%% Apply Robin boundary conditions
for iside = rbnSides
    dofs = solution.spaceST.boundary(iside).dofs;
    msh_side = msh_eval_boundary_side (solution.mshST,iside);
    sp_side = sp_eval_boundary_side (solution.spaceST,msh_side);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
    end
    rhs(dofs) = rhs(dofs) + op_f_v (sp_side,msh_side,gRbn(x{:}));
end
clear msh_side sp_side x dofs;

%% Solve the linear system
% case 'dir':    direct solver for the linear system 
% case 'dirEff': a faster direct solver based on a triangularization of the
%                time matrices. In presence of Robin boundary condition,
%                this solver cannot be used.
% case 'prec':   GMRES solver preconditioned by dirEff.                
switch Solver
    case 'dir'
        Mat=kron(Mt(1:end-1,2:end),Ks(intDofsS,intDofsS))-kron(Kt(1:end-1,2:end),Ms(intDofsS,intDofsS));
        if~isempty(rbnSides)
            Mat=Mat+kron(Wt(1:end-1,2:end),Mr(intDofsS,intDofsS));
        end
        solution.u(int_dofs_1) = Mat\rhs(int_dofs_2);
    case 'dirEff'
        [U,Ds,Q,Z,Tr,A] = FTparameters(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end));
        clear Ks Kt Ms Mt intDofsS;
        solution.u(int_dofs_1) = FTapplication(rhs(int_dofs_2),U,Ds,Q,Z,Tr,A);
        clear U Ds Q Z Tr A;
    case 'prec'
        [U,Ds,Q,Z,Tr,A] = FTparameters(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end));
        Prec=@(x) FTapplication(x,U,Ds,Q,Z,Tr,A);
        if exist('Mr','var')
            Mat = @(x) applyWaveMatrixR(x,Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end),Mr(intDofsS,intDofsS),Wt(1:end-1,2:end));
        else
            Mat = @(x) applyWaveMatrix(x,Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end));
        end
        if~exist('tol','var')
            tol=10^-12;
        end
        [y,~,relres,iter] = gmres(Mat,rhs(int_dofs_2),50,tol,numel(int_dofs_1),Prec);
        solution.u(int_dofs_1) = y;
        fprintf(1,'\n*Solver info: iter=[%d,%d]; Relres=%g*\n',iter(1),iter(2),relres);
        clear int_dofs_2 U Ds Q Z Tr A Ks Kt Ms Mt Wt Mr intDofsS iter relres;
end
solution.nDof=numel(int_dofs_1);
clear int_dofs_1 int_dofs_2 rhs;

%% Compute Errors
if computeError
    [errH1,errL2,errH1s] = sp_h1_c_error_tp (solution.spaceST, solution.mshST, solution.u, u_ex, grad_u_ex, c);
    [normH1,normL2,normH1s] = sp_h1_c_error_tp (solution.spaceST, solution.mshST, zeros(size(solution.u)), u_ex, grad_u_ex, c);
    solution.relErrH1=errH1/normH1;
    solution.relErrL2=errL2/normL2;
    solution.relErrH1s=errH1s/normH1s;
end

%% Compute Energy at t=0 and t=T
if computeEnergy
    dofs = solution.spaceST.boundary(end-1).dofs;
    msh_side = msh_eval_boundary_side (solution.mshST,numel(solution.spaceST.boundary)-1);
    sp_side = sp_precompute (solution.spaceST.boundary(end-1), msh_side,'value',true,'gradient',true);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
    end
    Ecoefs=@(x,y,t) cat (1,repmat(eye(solution.mshST.rdim-1,solution.mshST.rdim),[1,1,size(x)]),cat(2,repmat(zeros(1,solution.mshST.rdim-1),[1,1,size(x)]),reshape(1./c(x,y,t).^2,[1,1,size(x)])));
    Emat=op_energy(sp_side,msh_side,Ecoefs(x{:}));
    solution.Einit=solution.u(dofs)'*Emat*solution.u(dofs);

    dofs = solution.spaceST.boundary(end).dofs;
    msh_side = msh_eval_boundary_side (solution.mshST,numel(solution.spaceST.boundary));
    sp_side = sp_precompute (solution.spaceST.boundary(end), msh_side,'value',true,'gradient',true);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:),msh_side.nqn,msh_side.nel);
    end
    Emat=op_energy(sp_side,msh_side,Ecoefs(x{:}));
    clear msh_side sp_side x;
    solution.Efinal=solution.u(dofs)'*Emat*solution.u(dofs);
    clear dofs;
end
end

function gDrchlt = InitDrchlt2d(gDrchlt,gInit,x,y,t,dimS,iside)
if iside ==2*dimS+1
    gDrchlt=gInit(x,y);
else
    gDrchlt=gDrchlt(x,y,t);
end
end

function gDrchlt = InitDrchlt1d(gDrchlt,gInit,x,t,dimS,iside)
if iside ==2*dimS+1
    gDrchlt=gInit(x);
else
    gDrchlt=gDrchlt(x,t);
end
end