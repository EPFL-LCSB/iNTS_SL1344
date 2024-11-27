%% script to gap-fill model

% go to the iNTS_SL1344 folder
cd('/PathTo/iNTS_SL1344/')

% add paths
addpath(genpath('./'));
addpath(genpath('/PathTo/mattfa/'));
addpath(genpath('/PathTo/NICEgame/'));

%load files
load('./matFiles/CompartmentData.mat')
load('./matFiles/thermo_data.mat')
load('./models/iNTS_SL1344_v1_0.mat')
sourceModel = iNTS_SL1344_v1_0;
% load DB model
% here as an example DB3
load('.matFiles/DB3.mat')
DBmodel = DB3;

% merge the two models
[GFmodel, conflict] = PrepareForGapFilling(sourceModel, {DBmodel},{''}, 0,1,{},[],thermo_data);
% set a lower bound to growth rate
GFmodel.var_lb(find(ismember(GFmodel.varNames,strcat('F_',model.rxns(find(model.c)))))) = 0.001;

% define conditions for which we want to gap-fill
EX = find(contains(GFmodel.varNames,'R_EX_'));
GFmodel.var_ub(EX) = 0;
media = {'EX_cl_e';'EX_ca2_e';'EX_cobalt2_e';'EX_mobd_e';...
    'EX_cu2_e';'EX_fe2_e';'EX_fe3_e';'EX_h2o_e';'EX_k_e';'EX_mg2_e';...
    'EX_mn2_e';'EX_na1_e';'EX_nh4_e';'EX_pi_e';'EX_so4_e';...
    'EX_zn2_e';'EX_his__L_e'};

f = find(ismember(GFmodel.varNames,strcat('R_',media)));
GFmodel.var_ub(f) = 25;

if aerobicTag
    GFmodel.var_ub(find(ismember(GFmodel.varNames,'R_EX_o2_e'))) = 20;
    GFmodel.var_ub(f) = 10^-6;
else
    GFmodel.var_ub(find(ismember(GFmodel.varNames,'R_EX_o2_e'))) = 0;
end

% choose compound to be uptaken
% we give some examples here
a = {'C02107';'C00033';'C01018';'C01040';'C03752';'C00711';'C00148';...
    'C05593';'C00349';'C00163';'C11625';'C00483';'C0064'};
i = 1;
GFmodeli = GFmodel;
f = find(ismember(GFmodel.varNames,strcat('R_EX_',a{i},'_e')));
GFmodeli.var_ub(f) = 25;
%make sure the carbon source is uptaken
GFmodeli.var_lb(f) = 10^-6;
f = find(ismember(GFmodel.varNames,strcat('NF_EX_',a{i},'_e')));
GFmodeli.var_lb(f) = -25;

% check if we can gap-fill for this compound
sol = solveTFAmodelCplex(GFmodeli);
% if yes then proceed to gap-filling
% and generate the alternative solutions
numAlt = 100;
time = 10*60;
tagMin = 1; % solutions of min size
tagMinP1 = 1; % solutions of min+1 size
[resultStat,ActRxns,DPsAll] = gfbiomass(GFmodeli,GFmodeli.indUSE,numAlt,time,1,tagMinP1);
