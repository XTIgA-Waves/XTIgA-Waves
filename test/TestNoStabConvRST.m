%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution  of a Wave Equation  with a stabilized space-time method     %
%                                                                         %
%   d_tt u= c Lapl_x(u)   in Omega                                        %
%   u = g                 on GammaD                                       %
%   d_n u = h             on GammaN                                       %
%   u(0)=u_0                                                              %
%   u'(0)=u_1                                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Options
plotErrors = false;
saveTable = true;

%% Discretization Parameters
nsub = [8,16,32,64,128,256];
degree = [1:4]; %#ok<NBRAK2>

%% Problem Data
syms x t;
sym_c = 1;
theta = 1;
sym_u_ex = sin(pi*x)*sin(5/4*pi*t)^2;
problemData = computeProblemData1d(sym_u_ex,sym_c,theta);
problemData.geometryS = geo_load(nrbline(0,1));
problemData.T = 10;
problemData.drchltSides = [1,2]; % only for space domain
problemData.nmnnSides = []; % only for space domain
problemData.rbnSides = []; % only for space domain
problemData.perDir = []; % only for space domain

%% Preallocation
TabL2R0 = zeros(numel(nsub),numel(degree));
TabH1sR0 = zeros(numel(nsub),numel(degree));
TabL2R2 = zeros(numel(nsub),numel(degree));
TabH1sR2 = zeros(numel(nsub),numel(degree));
Tabht = zeros(numel(nsub),1);
TabL2R = zeros(numel(nsub),numel(degree));
TabH1sR = zeros(numel(nsub),numel(degree));

%% Loop
for ideg=1:numel(degree)
    methodData.degreeS = [degree(ideg)];
    methodData.degreeT = degree(ideg);
    methodData.regularityS = (methodData.degreeS-1);
    methodData.stab_param = 0;
    for isub=1:numel(nsub)
        methodData.nsubS = [nsub(isub)]; % a column for each spatial parametric direction
        methodData.nsubT = 2*nsub(isub);
        %
    	methodData.regularityS = 0;
        methodData.regularityT = 0;
        solutionR0 = solveWaveSTstabEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R0(isub,ideg) = solutionR0.relErrL2;
        TabH1sR0(isub,ideg) = solutionR0.relErrH1s;
        Tabht(isub,1) = solutionR0.ht;
        %
    	methodData.regularityS = max((methodData.degreeS-2),0);
        methodData.regularityT = max((methodData.degreeT-2),0);
        solutionR2 = solveWaveSTstabEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R2(isub,ideg) = solutionR2.relErrL2;
        TabH1sR2(isub,ideg) = solutionR2.relErrH1s;
        %
    	methodData.regularityS = (methodData.degreeS-1);
        methodData.regularityT = (methodData.degreeT-1);
        solutionR = solveWaveSTstabEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R(isub,ideg) = solutionR.relErrL2;
        TabH1sR(isub,ideg) = solutionR.relErrH1s;
    end
end

TabhpL2R = Tabht.^(degree+1).*TabL2R(end,:)./Tabht(end).^(degree+1);
TabhpH1sR = Tabht.^(degree).*TabH1sR(end,:)./Tabht(end).^(degree);

%% Save Table for Tikz
TabSumm = [Tabht,TabL2R0,TabH1sR0,TabL2R2,TabH1sR2,TabL2R,TabH1sR,TabhpL2R,TabhpH1sR];
nameVars = cell(1,8*numel(degree)+1);
nameVars{1} = 'ht';
for ideg = 1:numel(degree)
    nameVars{ideg+1+0*numel(degree)} = sprintf('L2P%dR0',degree(ideg));
    nameVars{ideg+1+1*numel(degree)} = sprintf('H1sP%dR0',degree(ideg));
    nameVars{ideg+1+2*numel(degree)} = sprintf('L2P%dR2',degree(ideg));
    nameVars{ideg+1+3*numel(degree)} = sprintf('H1sP%dR2',degree(ideg));
    nameVars{ideg+1+4*numel(degree)} = sprintf('L2P%dR',degree(ideg));
    nameVars{ideg+1+5*numel(degree)} = sprintf('H1sP%dR',degree(ideg));
    nameVars{ideg+1+6*numel(degree)} = sprintf('h%dL2R',degree(ideg)+1);
    nameVars{ideg+1+7*numel(degree)} = sprintf('h%dH1sR',degree(ideg));
end
ErrorTable = array2table(TabSumm,'VariableNames',nameVars);
if saveTable
    writetable(ErrorTable,strcat('ErrorTableNoStabConvRST.csv'));
end

%% Plot
if plotErrors
    CM = lines(numel(degree)); %#ok<*UNRCH>
    figL2 = figure('name','L^2 Rel','NumberTitle','off');
    xlabel('h');
    ylabel('L^2 Relative Error');
    hold on
    Leg1=[];
    Leg2=[];
    Leg3=[];
    for ideg = 1:numel(degree)
        ErrPlotR0(ideg) = loglog(Tabht,TabL2R0(:,ideg),'-.d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg1=cat(1,Leg1,cellstr(sprintf('p=%d, r=0',degree(ideg))));
        Leg2=cat(1,Leg2,cellstr(sprintf('O(h^%d)',degree(ideg)+1)));
        Leg3=cat(1,Leg3,cellstr(sprintf('p=%d, r=%d',degree(ideg),degree(ideg)-1)));
    end
    for ideg = 1:numel(degree)
        hPlotR(ideg) = loglog(Tabht,TabhpL2R(:,ideg),'--','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
    end
    for ideg = 1:numel(degree)
        ErrPlotR(ideg) = loglog(Tabht,TabL2R(:,ideg),':d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend([Leg1;Leg2;Leg3],'Location','northeast');
    %%
    figH1s = figure('name','H^1s Rel','NumberTitle','off');
    xlabel('h');
    ylabel('H^1 Relative Error');
    hold on
    Leg1=[];
    Leg2=[];
    Leg3=[];
    for ideg = 1:numel(degree)
        ErrPlotR0(ideg) = loglog(Tabht,TabH1sR0(:,ideg),'-.d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg1=cat(1,Leg1,cellstr(sprintf('p=%d, r=0',degree(ideg))));
        Leg2=cat(1,Leg2,cellstr(sprintf('O(h^%d)',degree(ideg))));
        Leg3=cat(1,Leg3,cellstr(sprintf('p=%d, r=%d',degree(ideg),degree(ideg)-1)));
    end
    for ideg = 1:numel(degree)
        hPlotR0(ideg) = loglog(Tabht,TabhpH1sR(:,ideg),'--','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
    end
    for ideg = 1:numel(degree)
        ErrPlotR(ideg) = loglog(Tabht,TabH1sR(:,ideg),':d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend([Leg1;Leg2;Leg3],'Location','northeast');
end
