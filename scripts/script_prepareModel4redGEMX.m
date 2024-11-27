%% script to prepare model for redGEMX analysis

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

% save folder
path_save = 'myPath/';

% merge model with a database
[GFmodel, conflict] = MergeModels(sourceModel, models,'c');

% define the compounds that you want to connect
ExtracellularMedium_connect = {'EX_itacon_e'
'EX_ach_e'
'EX_srtn_e'
'EX_nrpphr_e'
'EX_adrnl_e'
'EX_34dhoxmand_e'
'EX_melatn_e'
'EX_dchac_e'
'EX_cholate_e'
'EX_tchola_e'
'EX_C02528_e'
'EX_C04643_e'
'EX_C20978_e'
'EX_fol_e'
'EX_but_e'
'EX_frulys_e'
'EX_acmum_e'
'EX_dopa_e'
'EX_taur_e'
'EX_cbl1_e'};

sourceModel.ExtracellularMedium_connect =  ExtracellularMedium_connect;
sourceModel.ExtracellularMedium = ExtracellularMedium_connect;


% split the rxns in the model in 3 subsystems
% 1. Database 2.Model metabolic 3.Model transport and exchange 
[e,i] = getExchangeRxns(sourceModel);
t = (findTransRxns(sourceModel));
t = find(ismember(sourceModel.rxns,t));
sourceModel.subSystems = sourceModel.rxnFrom;
sourceModel.subSystems([i;t]) = {'Exchange and Transport'};
clear sourceModel.metCharge;
model = sourceModel;

% save the model for the redGEMX analysis
save(strcat(path_save,'model4redGEMX.mat'),'model')