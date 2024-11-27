%% script to convert to thermo structure with thermodynamic constraints

% go to the iNTS_SL1344 folder
cd('/PathTo/iNTS_SL1344/')

% add paths
addpath(genpath('./'));
addpath(genpath('/PathTo/mattfa/'));

%load files
load('./matFiles/CompartmentData.mat')
load('./matFiles/thermo_data.mat')
load('./models/iNTS_SL1344_v1_0.mat')

% set media
%close all uptakes and sinks
idx_exchange = find(contains(model.rxns,{'EX_','sink_';'DM_'}));
model.lb(idx_exchange) = 0;

media = {'EX_cl_e';'EX_ca2_e';'EX_cobalt2_e';'EX_mobd_e';...
    'EX_cu2_e';'EX_fe2_e';'EX_fe3_e';'EX_h2o_e';'EX_k_e';'EX_mg2_e';...
    'EX_mn2_e';'EX_na1_e';'EX_nh4_e';'EX_pi_e';'EX_so4_e';...
    'EX_zn2_e';'EX_glc__D_e';'EX_his__L_e';'EX_sel_e';'EX_o2_e';'EX_co2_e'};

f = find(ismember(model.rxns,media));
model.lb(f) = -25;

%convert model to thermo 
model = iNTS_SL1344_v1_0;
model.CompartmentData = CompartmentData;
tmodel = prepModelforTFA(model,thermo_data,model.CompartmentData);
ttmodel = convToTFA(tmodel, thermo_data, {'AIRC3';'GARFT';'UDCPPtppi'}, 'DGo', {}, 0.001, 1, 1);

% add net flux variables (optional)
ttmodel = addNetFluxVariables(ttmodel);
indNF = getAllVar(ttmodel,{'NF'});