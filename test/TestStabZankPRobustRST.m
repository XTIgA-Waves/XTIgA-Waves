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
plotErrors=false;
saveTable=true;

%% Discretization Parameters
nsubT=64;
nsubS=[4,8,16,32,64,128,256,512,1024];
degree=[1:4];

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
TabL2R0=zeros(numel(nsubS),numel(degree));
TabH1sR0=zeros(numel(nsubS),numel(degree));
Tabht=zeros(numel(nsubS),1);
Tabhs=zeros(numel(nsubS),1);
TabL2R2=zeros(numel(nsubS),numel(degree));
TabH1sR2=zeros(numel(nsubS),numel(degree));
TabL2R=zeros(numel(nsubS),numel(degree));
TabH1sR=zeros(numel(nsubS),numel(degree));

%% Loop
for ideg=1:numel(degree)
    methodData.degreeS=[degree(ideg)];
    methodData.degreeT=degree(ideg);
    methodData.regularityS = (methodData.degreeS-1);
    methodData.regularityT = (methodData.degreeT-1);
    for isubS=1:numel(nsubS)
        methodData.nsubS=nsubS(isubS); % a column for each spatial parametric direction
        methodData.nsubT=nsubT;
        %
    	methodData.regularityS = 0;
        methodData.regularityT = 0;
        solutionR0=solveWaveSTstabZankPEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R0(isubS,ideg)=solutionR0.relErrL2;
        TabH1sR0(isubS,ideg)=solutionR0.relErrH1s;
        Tabht(isubS,1)=solutionR0.ht;
        Tabhs(isubS,1)=solutionR0.hs;
        %
    	methodData.regularityS = max((methodData.degreeS-2),0);
        methodData.regularityT = max((methodData.degreeT-2),0);
        solutionR2=solveWaveSTstabZankPEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R2(isubS,ideg)=solutionR2.relErrL2;
        TabH1sR2(isubS,ideg)=solutionR2.relErrH1s;
        %
    	methodData.regularityS = (methodData.degreeT-1);
        methodData.regularityT = (methodData.degreeT-1);
        solutionR=solveWaveSTstabZankPEx(problemData,methodData,'Solver','dir','computeError',true);
        TabL2R(isubS,ideg)=solutionR.relErrL2;
        TabH1sR(isubS,ideg)=solutionR.relErrH1s;
    end
end

Tabhtfhs=Tabht./Tabhs;

%% Save Table for Tikz
TabSumm=[Tabhtfhs,TabL2R0,TabH1sR0,TabL2R2,TabH1sR2,TabL2R,TabH1sR];
nameVars=cell(1,6*numel(degree)+1);
nameVars{1}='htfhs';
for ideg=1:numel(degree)
    nameVars{ideg+1+0*numel(degree)}=sprintf('L2P%dR0',degree(ideg));
    nameVars{ideg+1+1*numel(degree)}=sprintf('H1P%dR0',degree(ideg));
    nameVars{ideg+1+2*numel(degree)}=sprintf('L2P%dR2',degree(ideg));
    nameVars{ideg+1+3*numel(degree)}=sprintf('H1P%dR2',degree(ideg));
    nameVars{ideg+1+4*numel(degree)}=sprintf('L2P%dR',degree(ideg));
    nameVars{ideg+1+5*numel(degree)}=sprintf('H1P%dR',degree(ideg));
end
ErrorTable = array2table(TabSumm,'VariableNames',nameVars);
if saveTable
    writetable(ErrorTable,strcat('ErrorTableStabZankPRobustRST.csv'));
end

%% Plot
if plotErrors
    CM = lines(numel(degree));
    figL2 = figure('name','L^2 Rel','NumberTitle','off');
    xlabel('h_t/h_s');
    ylabel('L^2 Relative Error');
    hold on
    Leg1=[];
    Leg2=[];
    for ideg = 1:numel(degree)
        ErrPlotR0(ideg) = plot(Tabhtfhs,TabL2R0(:,ideg),'-.d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg1=cat(1,Leg1,cellstr(sprintf('p=%d, r=0',degree(ideg))));
    end
    for ideg = 1:numel(degree)
        ErrPlotR(ideg) = plot(Tabhtfhs,TabL2R(:,ideg),':d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg2=cat(1,Leg2,cellstr(sprintf('p=%d, r=%d',degree(ideg),degree(ideg)-1)));
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend([Leg1;Leg2],'Location','SouthWest');
    %%
    figH1 = figure('name','H^1 seminorm Rel','NumberTitle','off');
    xlabel('h_t/h_s');
    ylabel('H^1 Relative Error');
    hold on
    Leg1=[];
    Leg2=[];
    for ideg = 1:numel(degree)
        ErrPlotR0(ideg) = plot(Tabhtfhs,TabH1sR0(:,ideg),'-.d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg1=cat(1,Leg1,cellstr(sprintf('p=%d, r=0',degree(ideg))));
    end
    for ideg = 1:numel(degree)
        ErrPlotR(ideg) = plot(Tabhtfhs,TabH1sR(:,ideg),':d','color',CM(ideg,:),'MarkerSize',6,'LineWidth',1.5);
        Leg2=cat(1,Leg2,cellstr(sprintf('p=%d, r=%d',degree(ideg),degree(ideg)-1)));
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend([Leg1;Leg2],'Location','SouthWest');
end
