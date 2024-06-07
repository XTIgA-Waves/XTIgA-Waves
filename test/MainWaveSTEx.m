%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution with a stabilized space-time method of the Wave Equation:    %
%                                                                         %
%  d_tt u(x,t) - div(c(x) grad(u(x,t))) = f(x,t), (x,t) in Omega x [0,T]  %
%  u(x,t) = gDrchlt(x,t)                          (x,t) in GammaD x [0,T] %
%  d_n u(x,t) = gNmnn(x,t)                        (x,t) in GammaN x [0,T] %
%  theta c(x) u'(x,t) + c(x)^2 d_n u = gRbn(x,t)  (x,t) in GammaR x [0,T] %
%  u(x,0) = u_0(x)                                x in Omega              %
%  u'(x,0) = u_1(x)                               x in Omega              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options
plotSolution = true;

%% Problem Data
syms x t;
sym_c = 1;
theta = 1;
sym_u_ex = sin(pi*x)*sin(5/4*pi*t)^2;
problemData = computeProblemData1d(sym_u_ex,sym_c,theta);
problemData.geometryS = geo_load(nrbline([0,0],[1,0]));
problemData.T = 2;
problemData.drchltSides = [1,2]; % only for space domain
problemData.nmnnSides = []; % only for space domain
problemData.rbnSides = []; % only for space domain
problemData.perDir = []; % only for space domain

%% Method data
methodData.nsubS = [1]*16; % a column for each spatial parametric direction
methodData.nsubT = 32;
methodData.degreeS = [1]*2;
methodData.degreeT = 2;
methodData.regularityS = (methodData.degreeS-1);
methodData.regularityT = (methodData.degreeT-1);
methodData.stab_param = 10^(-methodData.degreeT);
nrows = 3;
ncols = 3;

%% Solving the problem
[solution] = solveWaveSTstab(problemData,methodData,'Solver','dir','ExactSol',true,'Error',true,'Energy',false);

%% Plot of the solution
if plotSolution
    fig1 = plot2dSolutionTimeSection(solution.u,solution.spaceST,solution.geometryST,nrows,ncols);
%   fig2 = saveSolutionVideo(solution.u,solution.spaceST,solution.geometryST,'ExampleWaveAnimation',100,100);
end