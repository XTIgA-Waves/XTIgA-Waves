function solution = solveWaveSTMPstab(problemData,methodData,varargin)

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

%% Spatial structures
[geometryS,bnds,intfs,~,bndsIntfs] = mp_geo_load(geoName);
nPtc = numel(geometryS);

mshSptc = cell(nPtc,1);
spaceSptc = cell(nPtc,1);
dimS = geometryS.rdim;

for iPtc = 1:nPtc
    nsubSptc = ceil((nsubS./(geometryS(iPtc).nurbs.number-geometryS(iPtc).nurbs.order+1))-1);
    [rknotsS,zetaS] = kntrefine (geometryS(iPtc).nurbs.knots,nsubSptc,degreeS,regularityS);
    quadRule = msh_gauss_nodes(degreeS+1);
    [qn,qw] = msh_set_quad_nodes(zetaS,quadRule);
    mshSptc{iPtc} = msh_cartesian(zetaS,qn,qw,geometryS(iPtc));
    spaceSptc{iPtc} = sp_bspline(rknotsS,degreeS,mshSptc{iPtc});
end
mshS = msh_multipatch(mshSptc,bnds);
spaceS = sp_multipatch(spaceSptc,mshS,intfs,bndsIntfs);
clear spaceSptc mshSptc;
solution.hs = max(arrayfun(@(i)max(msh_precompute(mshS.msh_patch{i}).element_size),...
    1:numel(mshS.msh_patch)));

%% Time structures
geometryT = geo_load(nrbline([0,0],[T,0]));
nsubT = ceil((nsubT./(geometryT.nurbs.number-geometryT.nurbs.order+1))-1);
[rknotsT,zetaT] = kntrefine(geometryT.nurbs.knots,nsubT,degreeT,regularityT);
clear nsubS nsubT;
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
        Ks = op_gradu_gradv_mp(spaceS,spaceS,mshS,@(x) c(x,0).^(2));
    case 2
        Ks = op_gradu_gradv_mp(spaceS,spaceS,mshS,@(x,y) c(x,y,0).^(2));
    case 3
        Ks = op_gradu_gradv_tp(spaceS,spaceS,mshS,@(x,y,z) c(x,y,z,0).^(2));
end
Ms = op_u_v_mp(spaceS,spaceS,mshS);

%% Assembling Time matrices
Mt = op_u_v_tp(spaceT,spaceT,mshT);
Kt = op_gradu_gradv_tp(spaceT,spaceT,mshT);

if ~isempty(rbnSides)
    Wt = op_vel_dot_gradu_v_tp(spaceT,spaceT,mshT,@(t) ones(size(t)));
    nnzMr = sum(arrayfun(@(ibnd) spaceS.boundary.ndof_per_patch(ibnd),rbnSides))*(2*max(degreeS)+1);
    Mr = spalloc(size(Ms,1),size(Ms,2),nnzMr);
    for iSide = rbnSides
        dofs = spaceS.boundary.dofs(spaceS.boundary.gnum{iSide});
        if mshS.ndim == 1
            Mr(dofs,dofs) = theta*c(mshS.map(iSide-1),0);
        else
            mshSide = msh_eval_boundary_side (mshS.msh_patch{mshS.boundaries(iSide).patches},mshS.boundaries(iSide).faces);
            spSide = sp_precompute (spaceS.boundary.sp_patch{iSide}, mshSide,'value',true,'gradient',false);
            x = cell(mshS.rdim,1);
            for idim = 1:mshS.rdim
                x{idim} = reshape (mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
            end
            Mr(dofs,dofs) = Mr(dofs,dofs) + op_u_v(spSide,spSide,mshSide, theta*c(x{1},x{2},zeros(size(x{1}))));
        end
    end
end
clear mshSide spSide nnzMr dofs;

if ~isempty(stab_param)
    if stab_param>0
        shapeFunctionsDerT = sp_bspline_1d_param_all_derivatives (rknotsT,degreeT,mshT.qn{1});
        sizeReshapeT = size(shapeFunctionsDerT);
        sizeReshapeT = sizeReshapeT(2:end);
        shapeFunctionsDerT = reshape(shapeFunctionsDerT(degreeT+1,:,:,:),sizeReshapeT);
        clear sizeReshapeT;
        Pt = stab_param*solution.ht^(2*degreeT)*op_mat_hostab_tp(spaceT, spaceT, mshT, shapeFunctionsDerT,shapeFunctionsDerT,degreeT);
        clear shapeFunctionsDerT stab_param spaceT;
    end
end
clear mshT;

%% Assembling Global structures
solution.geometryST = arrayfun(@(geo) geo_load_st(geo,0,T),geometryS);
mshSTptc = cell(nPtc,1);
spaceSTptc = cell(nPtc,1);
for iPtc = 1:nPtc
    quadRule = msh_gauss_nodes([degreeS,degreeT]+1);
    [qn,qw] = msh_set_quad_nodes([mshS.msh_patch{iPtc}.breaks(:)',{zetaT}],quadRule);
    mshSTptc{iPtc} = msh_cartesian([mshS.msh_patch{iPtc}.breaks(:)',{zetaT}],qn,qw,solution.geometryST(iPtc));
    spaceSTptc{iPtc} = sp_bspline([spaceS.sp_patch{iPtc}.knots(:)',{rknotsT}],[degreeS,degreeT],mshSTptc{iPtc});
    ePtc(iPtc) = nrbextrude(geometryS(iPtc).nurbs,[0,0,1]);
end
[~,bndsST,intfsST,~,bndsIntfsST] = mp_geo_load(ePtc);
clear ePtc;
solution.mshST = msh_multipatch(mshSTptc,bndsST);
solution.spaceST = sp_multipatch(spaceSTptc,solution.mshST,intfsST,bndsIntfsST);
clear geometryS quadRule qn qw zetaS zetaT rknotsT;
rhs = op_f_v_mp(solution.spaceST,solution.mshST,f);
solution.u = zeros(solution.spaceST.ndof,1);

%% Finding boundary side for t=0 and t=T
initBndsST = zeros(1,nPtc);
finalBndsST = zeros(1,nPtc);
iptcI=1;
iptcF=1;
for ibnd = 1:numel(bndsST)
    if bndsST(ibnd).faces == 5
        initBndsST(iptcI) = ibnd;
        iptcI = iptcI + 1;
    elseif bndsST(ibnd).faces == 6
        finalBndsST(iptcF) = ibnd;
        iptcF = iptcF + 1;
    end
end
clear iptcF iptcI;

%% Apply Dirichlet boundary conditions and initial condition (u(0)=u_0)
switch dimS
    case 1
        if ExactSol
            gDrchlt = @(x,t,iside) u_ex(x,t);
        else
            gDrchlt = @(x,t,iside) InitDrchlt1d(gDrchlt,gInit,x,t,initBndsST,iside);
        end
    case 2
        if ExactSol
            gDrchlt = @(x,y,t,iside) u_ex(x,y,t);
        else
            gDrchlt = @(x,y,t,iside) InitDrchlt2d(gDrchlt,gInit,x,y,t,initBndsST,iside);
        end
    case 3
        if ExactSol
            gDrchlt = @(x,y,z,t,iside) u_ex(x,y,z,t);
        else
            gDrchlt = @(x,y,z,t,iside) InitDrchlt3d(gDrchlt,gInit,x,y,z,t,initBndsST,iside);
        end
end

drchltSidesST = zeros(size(drchltSides));
for iSide = 1:numel(drchltSides)
    drchltSidesST(iSide) = find(ismember([bndsST.patches],bnds(drchltSides(iSide)).patches) & ismember([bndsST.faces],bnds(drchltSides(iSide)).faces));
end
drchltSidesST = reshape(drchltSidesST,1,[]);

nmnnSidesST = zeros(size(nmnnSides));
for iSide = 1:numel(nmnnSides)
    nmnnSidesST(iSide) = find(ismember([bndsST.patches],bnds(nmnnSides(iSide)).patches) & ismember([bndsST.faces],bnds(nmnnSides(iSide)).faces));
end
nmnnSidesST = reshape(nmnnSidesST,1,[]);

rbnSidesST = zeros(size(rbnSides));
for iSide = 1:numel(rbnSides)
    rbnSidesST(iSide) = find(ismember([bndsST.patches],bnds(rbnSides(iSide)).patches) & ismember([bndsST.faces],bnds(rbnSides(iSide)).faces));
end
rbnSidesST = reshape(rbnSidesST,1,[]);

[uForced, trialForcedDofs] = sp_drchlt_l2_proj(solution.spaceST,solution.mshST,gDrchlt,union(drchltSidesST,initBndsST));
solution.u(trialForcedDofs) = uForced;
clear uForced;

testForcedDofs = [];
for iSide = union(drchltSidesST,finalBndsST)
    testForcedDofs = union (testForcedDofs,solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{iSide}));
end
clear finalBndsST;

trialIntDofs = setdiff (1:solution.spaceST.ndof, trialForcedDofs);
testIntDofs = setdiff (1:solution.spaceST.ndof, testForcedDofs);

drchltDofsST = [];
drchltDofsS = [];
for iSide = 1:numel(drchltSides)
    drchltDofsST = union(drchltDofsST,solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{drchltSidesST(iSide)}));
    drchltDofsS = union(drchltDofsS,spaceS.boundary.dofs(spaceS.boundary.gnum{drchltSides(iSide)}));
end
intDofsS = setdiff(1:spaceS.ndof,drchltDofsS);
solution.nDofS = numel(intDofsS);

initDofsST = [];
for iSide = 1:numel(initBndsST)
    initDofsST = union (initDofsST,solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{initBndsST(iSide)}));
end
drchltDofsST = setdiff(drchltDofsST,initDofsST);
clear spaceS;

%% Reorder global Dofs for exploting Kronecker products
gnumIntfsST = arrayfun(@(iInt) reshape(solution.spaceST.gnum{solution.spaceST.interfaces(iInt).patch1}(solution.spaceST.sp_patch{solution.spaceST.interfaces(iInt).patch1}.boundary(solution.spaceST.interfaces(iInt).side1).dofs),[],solution.spaceST.sp_patch{1}.ndof_dir(end)),1:numel(intfs),'UniformOutput',false);
intfsDofsST = cell2mat(gnumIntfsST);
gnumIntST = arrayfun(@(iPtc) reshape(setdiff(solution.spaceST.gnum{iPtc},intfsDofsST),[],solution.spaceST.sp_patch{1}.ndof_dir(end)),1:nPtc,'UniformOutput',false);
newDofs= [];
for iT=1:solution.spaceST.sp_patch{1}.ndof_dir(end)
    for iPtc=1:nPtc
        newDofs = cat(1,newDofs,gnumIntST{iPtc}(:,iT));
    end
    for iInt=1:numel(intfs)
        newDofs = cat(1,newDofs,gnumIntfsST{iInt}(:,iT));
    end
end
clear gnumIntfsST intfsDofsST gnumIntST;
newDofs = unique(newDofs,'stable');
oldDofs = (1:numel(newDofs))';
oldDofs(newDofs) = oldDofs;
% newInitDofs = arrayfun(@(x) find(sort(newDofs(initDofsST)) == x, 1), newDofs(initDofsST));
% oldInitDofs = (1:numel(newInitDofs))';
% oldInitDofs(newInitDofs) = oldInitDofs;
% newDrchltDofs = arrayfun(@(x) find(sort(newDofs(drchltDofsST)) == x, 1), newDofs(drchltDofsST));
% oldDrchltDofs = (1:numel(newDrchltDofs))';
% oldDrchltDofs(newDrchltDofs) = oldDrchltDofs;
% newTestIntDofs = arrayfun(@(x) find(sort(newDofs(testIntDofs)) == x, 1), newDofs(testIntDofs));
% oldTestIntDofs = (1:numel(newTestIntDofs))';
% oldTestIntDofs(newTestIntDofs) = oldTestIntDofs;
% newTrialIntDofs = arrayfun(@(x) find(sort(newDofs(trialIntDofs)) == x, 1), newDofs(trialIntDofs));
% oldTrialIntDofs = (1:numel(newTrialIntDofs))';
% oldTrialIntDofs(newTrialIntDofs) = oldTrialIntDofs;
%
[~,sortOldTest] = sort(oldDofs(testIntDofs));
sortOldTestT = (1:numel(sortOldTest))';
sortOldTestT(sortOldTest) = sortOldTestT;
[~,sortOldTrial] = sort(oldDofs(trialIntDofs));
sortOldTrialT = (1:numel(sortOldTrial))';
sortOldTrialT(sortOldTrial) = sortOldTrialT;
[~,sortOldDrchlt] = sort(oldDofs(drchltDofsST));
[~,sortOldInit] = sort(oldDofs(initDofsST));
%
u1 = solution.u(initDofsST);
u2 = solution.u(drchltDofsST);

tmp0 = Ms(intDofsS,:)*u1(sortOldInit)*Kt(1:end-1,1)'...
    -Ks(intDofsS,:)*u1(sortOldInit)*Mt(1:end-1,1)'...
    +Ks(intDofsS,:)*u1(sortOldInit)*Pt(1:end-1,1)';

tmpDrchlt = Ms(intDofsS,drchltDofsS)*reshape(u2(sortOldDrchlt),numel(drchltDofsS),size(Kt,2)-1)*Kt(1:end-1,2:end)'...
    -Ks(intDofsS,drchltDofsS)*reshape(u2(sortOldDrchlt),numel(drchltDofsS),size(Kt,2)-1)*Mt(1:end-1,2:end)'...
    +Ks(intDofsS,drchltDofsS)*reshape(u2(sortOldDrchlt),numel(drchltDofsS),size(Kt,2)-1)*Pt(1:end-1,2:end)';

if ~isempty(rbnSides)
    tmp0 = tmp0 -Mr(intDofsS,:)*u1(sortOldInit)*Wt(1:end-1,1)';
    tmpDrchlt = tmpDrchlt -Mr(intDofsS,drchltDofsS)*reshape(u2(sortOldDrchlt),numel(drchltDofsS),size(Kt,2)-1)*Wt(1:end-1,2:end)';
end

tmp0 = tmp0(:);
tmpDrchlt = tmpDrchlt(:);

rhs(testIntDofs) = rhs(testIntDofs) + tmp0(sortOldTestT) + tmpDrchlt(sortOldTestT);
clear trialForcedDofs testForcedDofs tmp0 tmpDrchlt drchltDofsS drchltDofsST;

%% Apply Initial conditions (u'(0)=u_1)
for iSide = 1:numel(initBndsST)
    dofs = solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{initBndsST(iSide)});
    mshSide = msh_eval_boundary_side(solution.mshST.msh_patch{solution.mshST.boundaries(initBndsST(iSide)).patches},solution.mshST.boundaries(initBndsST(iSide)).faces);
    spSide = sp_precompute(solution.spaceST.boundary.sp_patch{initBndsST(iSide)},mshSide,'value',true,'gradient',false);
    x = cell(solution.mshST.rdim,1);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    rhs(dofs) = rhs(dofs) + op_f_v(spSide,mshSide,gInitDer(x{:}));
end
clear spSide mshSide x initBndsST;

%% Apply Neumann boundary conditions
for iSide = 1:numel(nmnnSidesST)
    dofs = solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{nmnnSidesST(iSide)});
    mshSide = msh_eval_boundary_side(solution.mshST.msh_patch{solution.mshST.boundaries(nmnnSidesST(iSide)).patches},solution.mshST.boundaries(nmnnSidesST(iSide)).faces);
    spSide = sp_precompute(solution.spaceST.boundary.sp_patch{nmnnSidesST(iSide)}, mshSide,'value',true,'gradient',false);
    x = cell(solution.mshST.rdim,1);
    for idim = 1:solution.mshST.rdim
        x{idim} = reshape(mshSide.geo_map(idim,:,:),mshSide.nqn,mshSide.nel);
    end
    if ExactSol
        rhs(dofs) = rhs(dofs) + op_fdotn_v(spSide,mshSide,grad_u_ex_c2(x{:}));
    else
        rhs(dofs) = rhs(dofs) + op_f_v(spSide,mshSide,gNmnn(x{:}));
    end
end
clear spSide mshSide x dofs nmnnSidesST;

%% Apply Robin boundary conditions
for iSide = 1:numel(rbnSidesST)
    dofs = solution.spaceST.boundary.dofs(solution.spaceST.boundary.gnum{rbnSidesST(iSide)});
    mshSide = msh_eval_boundary_side(solution.mshST.msh_patch{solution.mshST.boundaries(rbnSidesST(iSide)).patches},solution.mshST.boundaries(rbnSidesST(iSide)).faces);
    spSide = sp_precompute(solution.spaceST.boundary.sp_patch{rbnSidesST(iSide)}, mshSide,'value',true,'gradient',false);
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
clear mshSide spSide x dofs rbnSidesST;

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
        rhs = rhs(testIntDofs);
        y = Mat\rhs(sortOldTest);
        solution.u(trialIntDofs) = y(sortOldTrialT);
    case 'dirEff'
        [U,Ds,Q,Z,Tr,A] = FTparameters(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        clear Ks Kt Ms Mt intDofsS;
        rhs = rhs(testIntDofs);
        y = FTapplication(rhs(sortOldTest),U,Ds,Q,Z,Tr,A);
        solution.u(trialIntDofs) = y(sortOldTrialT);
        clear U Ds Q Z Tr A;
    case 'precFDT'
        sTimeSetup = tic;
        [U,Ds,Q,Z,Tr,A] = fdtParam(Ks(intDofsS,intDofsS),Kt(1:end-1,2:end),Ms(intDofsS,intDofsS),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        solution.timeSetup = toc(sTimeSetup);
        Prec=@(x) fdtApply(x,U,Ds,Q,Z,Tr,A);
        if exist('Mr','var')
            Mat = @(x) applyWaveMatrix(x,Kt(1:end-1,2:end),Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Ks(intDofsS,intDofsS),Ms(intDofsS,intDofsS),Mr(intDofsS,intDofsS),Wt(1:end-1,2:end));
        else
            Mat = @(x) applyWaveMatrix(x,Kt(1:end-1,2:end),Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Ks(intDofsS,intDofsS),Ms(intDofsS,intDofsS));
        end
        if~exist('tol','var')
            tol=10^-12;
        end
        rhs = rhs(testIntDofs);
        sTimeSolver = tic;
        [y,~,solution.relres,solution.iter] = gmres(Mat,rhs(sortOldTest),100,tol,numel(trialIntDofs),Prec);
        solution.timeSolver = toc(sTimeSolver);
        solution.u(trialIntDofs) = y(sortOldTrialT);
        fprintf(1,'\n*Solver info: iter=[%d,%d]; Relres=%g*\n',solution.iter(1),solution.iter(2),solution.relres);
        clear testIntDofs U Ds Q Z Tr A Ks Kt Ms Mt Mr Wt intDofsS iter relres;
    case 'precFT'
        sTimeSetup = tic;
        [Q,Z,Tr,A] = ftParam(Kt(1:end-1,2:end),Mt(1:end-1,2:end)-Pt(1:end-1,2:end));
        solution.timeSetup = toc(sTimeSetup);
        Prec=@(x) ftApply(x,Q,Z,Tr,A,Ks(intDofsS,intDofsS),Ms(intDofsS,intDofsS));
        if exist('Mr','var')
            Mat = @(x) applyWaveMatrix(x,Kt(1:end-1,2:end),Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Ks(intDofsS,intDofsS),Ms(intDofsS,intDofsS),Mr(intDofsS,intDofsS),Wt(1:end-1,2:end));
        else
            Mat = @(x) applyWaveMatrix(x,Kt(1:end-1,2:end),Mt(1:end-1,2:end)-Pt(1:end-1,2:end),Ks(intDofsS,intDofsS),Ms(intDofsS,intDofsS));
        end
        if~exist('tol','var')
            tol=10^-12;
        end
        rhs = rhs(testIntDofs);
        sTimeSolver = tic;
        [y,~,solution.relres,solution.iter] = gmres(Mat,rhs(sortOldTest),100,tol,numel(trialIntDofs),Prec);
        solution.timeSolver = toc(sTimeSolver);
        solution.u(trialIntDofs) = y(sortOldTrialT);
        fprintf(1,'\n*Solver info: iter=[%d,%d]; Relres=%g*\n',solution.iter(1),solution.iter(2),solution.relres);
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
    [errH1, errL2, errH1s] = sp_h1_error_spline_ex_mp(solution.spaceST,solution.mshST,solution.u,solFine.spaceST,solution.geometryST,solFine.u,c);
    [normH1,normL2,normH1s] = sp_h1_error_spline_ex_mp(solution.spaceST,solution.mshST,solution.u*0,solFine.spaceST,solution.geometryST,solFine.u,c);
    solution.relErrH1=errH1/normH1;
    solution.relErrL2=errL2/normL2;
    solution.relErrH1s=errH1s/normH1s;
end

%% Compute Energy at t=0 and t=T
if Energy
    solution.E0Et = computeEnergy(solution.mshST.qn,solution.mshST.qw,solution.u,c,solution.spaceST,solution.geometryST,solution.mshST.map_der,[0,1]);
end

end


function gDrchlt = InitDrchlt3d(gDrchlt,gInit,x,y,z,t,initBndsST,iside)
if ismember(iside,initBndsST)
    gDrchlt = gInit(x,y,z);
else
    gDrchlt = gDrchlt(x,y,z,t);
end
end

function gDrchlt = InitDrchlt2d(gDrchlt,gInit,x,y,t,initBndsST,iside)
if ismember(iside,initBndsST)
    gDrchlt = gInit(x,y);
else
    gDrchlt = gDrchlt(x,y,t);
end
end

function gDrchlt = InitDrchlt1d(gDrchlt,gInit,x,t,initBndsST,iside)
if ismember(iside,initBndsST)
    gDrchlt = gInit(x);
else
    gDrchlt = gDrchlt(x,t);
end
end