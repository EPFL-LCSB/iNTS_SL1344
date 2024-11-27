
% add paths
addpath(genpath('/yourPath/phenomapping/phenomapping/'))
addpath(genpath('/yourPath/mattfa/'))
addpath(genpath('/CPLEX_PATH'))

% save folder
path_save = 'YourPath';
mkdir(path_save)

load('./models/iNTS_SL1344_v1_0.mat')
model = iNTS_SL1344_v1_0;
model.description = 'iNTS_SL1344';
modeldescription = 'iNTS_SL1344';

%% generate iMMs
grRate = 0.12;          % optimal growth rate (this can be obtained optimizing for growth)
minObj = 0.99*grRate;    % minimal required growth
flagUptSec = 0;            % false for analysis of in silico minimal media (IMM). True for analysis of in silico minimal exchange (IME)
maxObj = 10;            % upper bound in growth (not applied here since 10 is much higher than the growth yield of any GEM)
drainsForiMM = {};      % names of drains/exchanges used for IMM analysis. Empty means by default it will search for minimal media accross all drains/exchanges
metabData = [];         % provide metabolomics data to integrate in this analysis. Empty if none.
NumAlt = 1000;             % number of alternative minimal media to find. We define 1 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
time = 600;              % time limit for optimization in seconds. If empty we do not set a time limit
tagMin = 1;             % additional constrain to avoid generating suboptimal solutions, meaning that if 1 we will not identify media that is not minimal
filename = strcat(modeldescription,'_PhenoMappingSubstrates');

% IMM analysis
% Define milp problem
modelmilp = analysisIME(model, flagUptSec, minObj, maxObj, drainsForiMM,metabData,time);
% Get alternative solutions for min size
[DPsimm,modelmilp] = findDPMax(modelmilp, NumAlt, modelmilp.indUSE,time, tagMin, strcat(path_save,filename));

%save results
save(strcat(path_save,filename,'.mat'))
r1 = load(strcat(path_save,filename,'.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfoIMMDPs_IME(r1.modelmilp, r1.DPsimm, ...
    r1.modelmilp.indUSE,  r1.modelmilp.rowDrain_F,...
    r1.modelmilp.rowDrain_R,r1.modelmilp.F_Flux,  r1.modelmilp.R_Flux) ;


immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
    num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];

writeData(strcat(path_save,filename,'.csv'), immOutput.summary,...
    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
    '%s\t%s\t%s\t%s\t%s');

save(strcat(path_save,'Active.mat'),'Active');
clear Active


% Get alternative solutions for min+1 size
f = find(contains(modelmilp.constraintNames,'CUT_0'));
modelmilp.rhs(f) = 0;
filename = strcat(modeldescription,'_PhenoMappingSubstrates_MinP1');

[DPsimm,modelmilp] = findDPMax(modelmilp, NumAlt, modelmilp.indUSE,time, tagMin, strcat(path_save,filename));

%save results
save(strcat(path_save,filename,'.mat'))
r1 = load(strcat(path_save,filename,'.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfoIMMDPs_IME(r1.modelmilp, r1.DPsimm, ...
    r1.modelmilp.indUSE,  r1.modelmilp.rowDrain_F,...
    r1.modelmilp.rowDrain_R,r1.modelmilp.F_Flux,  r1.modelmilp.R_Flux) ;


immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
    num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];

writeData(strcat(path_save,filename,'.csv'), immOutput.summary,...
    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
    '%s\t%s\t%s\t%s\t%s');

save(strcat(path_save,'ActiveMinP1.mat'),'Active');
clear Active


% Get alternative solutions for min+2 size
f = find(contains(modelmilp.constraintNames,'CUT_0'));
modelmilp.rhs(f) = 0;
filename = strcat(modeldescription,'_PhenoMappingSubstrates_MinP2');

[DPsimm,modelmilp] = findDPMax(modelmilp, NumAlt, modelmilp.indUSE,time, tagMin, strcat(path_save,filename));

%save results
save(strcat(path_save,filename,'.mat'))
r1 = load(strcat(path_save,filename,'.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfoIMMDPs_IME(r1.modelmilp, r1.DPsimm, ...
    r1.modelmilp.indUSE,  r1.modelmilp.rowDrain_F,...
    r1.modelmilp.rowDrain_R,r1.modelmilp.F_Flux,  r1.modelmilp.R_Flux) ;


immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
    num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];

writeData(strcat(path_save,filename,'.csv'), immOutput.summary,...
    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
    '%s\t%s\t%s\t%s\t%s');

save(strcat(path_save,'ActiveMinP2.mat'),'Active');
clear Active

% Get alternative solutions for min+3 size
f = find(contains(modelmilp.constraintNames,'CUT_0'));
modelmilp.rhs(f) = 0;
filename = strcat(modeldescription,'_PhenoMappingSubstrates_MinP3');

[DPsimm,modelmilp] = findDPMax(modelmilp, NumAlt, modelmilp.indUSE,time, tagMin, strcat(path_save,filename));

%save results
save(strcat(path_save,filename,'.mat'))
r1 = load(strcat(path_save,filename,'.mat'));

% Extract info about composition of the IMMs
[immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
    immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfoIMMDPs_IME(r1.modelmilp, r1.DPsimm, ...
    r1.modelmilp.indUSE,  r1.modelmilp.rowDrain_F,...
    r1.modelmilp.rowDrain_R,r1.modelmilp.F_Flux,  r1.modelmilp.R_Flux) ;


immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
    num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];

writeData(strcat(path_save,filename,'.csv'), immOutput.summary,...
    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
    '%s\t%s\t%s\t%s\t%s');

save(strcat(path_save,'ActiveMinP3.mat'),'Active');
clear Active

%% evaluate iMMs with respect to Biomass yield

model = iNTS_SL1344_v1_0;

%close all uptakes
idx_ex = find(contains(model.varNames,'R_EX_'));
model.var_ub(idx_ex) = 0;

% constraint carbon and nitrogen uptakes
% get carbon sources
f_ex = getCarbonSources(model,'tfa');
exchange = strrep(model.varNames(f_ex),'R_EX_','EX_');
for i = 1:length(exchange)
    f1 = find(ismember(model.rxns,exchange{i}));
    metsIdx(i,1) = find(model.S(:,f1));
end

%get carbon atoms
CarbonAtoms = getCarbonAtoms(model,metsIdx);
%
[NumCons, NumVars] = size(model.A);
NewCons = zeros(NumVars,1);
NewCons(f_ex) = CarbonAtoms;
model.A(NumCons+1,:) = NewCons;
model.rhs(NumCons+1) = 60; %10 mmol/gDW h of glucose
model.constraintNames{NumCons+1} = ['CUT_pUP_carbon'];model.constraintType{NumCons+1} = '<';

%get nitrogen sources
f_ex = getNitrogenSources(model,'tfa');
exchange = strrep(model.varNames(f_ex),'R_EX_','EX_');
for i = 1:length(exchange)
    f1 = find(ismember(model.rxns,exchange{i}));
    metsIdxN(i,1) = find(model.S(:,f1));
end

%get nitrogen atoms
NitrogenAtoms = getNitrogenAtoms(model,metsIdxN);

[NumCons, NumVars] = size(model.A);
NewCons = zeros(NumVars,1);
NewCons(f_ex) = NitrogenAtoms;
model.A(NumCons+1,:) = NewCons;
model.rhs(NumCons+1) =1000; %1 mmol nh4 / 5 mmol glucose
model.constraintNames{NumCons+1} = ['CUT_pU_Pnitrogen'];
model.constraintType{NumCons+1} = '<';


A_Active = [];
Sizes = {'';'MinP1';'MinP2';'MinP3'};
for i = 1:length(Sizes)
    dataPath = strcat(path_save,'/',Sizes{i},'Active.mat');
    load(dataPath)
    A_Active = [A_Active Active];
end

load(path_save,'modelmilp')

aux = sum(A_Active,1);
min_size = -max(aux);
max_size = -min(aux);
num_alts = size(A_Active,2);

M = zeros(max_size,num_alts);
Active_evaluated = cell(max_size,1);

%calculate the max growth rate for each iMM for the same total carbon and
%nitrogen upatke
for i = 1:num_alts

    modeli = model;

    alt = modelmilp.varNames(modelmilp.indUSE(find(A_Active(:,i)==-1)));

    iMM_R = strrep(alt,'BFUSE_R_','');
    A(1:length(alt),i) = printRxnFormula(modelmilp,iMM_R,1,1,1);
    iMM_R = strrep(alt,'BFUSE_','');
    idx = find(ismember(modeli.varNames,iMM_R));
    modeli.var_ub(idx) = 25;
    idx_his = find(ismember(modeli.varNames,'R_EX_his__L_e'));
    modeli.var_ub(idx_his) = 0.02;
    soli = solveTFAmodelCplex(modeli);
    if not(isempty(soli.val))
        sol(i,1) = soli.val;
    else
        sol(i,1) = nan;
        soli.x = zeros(length(model.varNames),1);
    end

    Active_evaluated(max_size+2,i) = num2cell(sol(i,1));
    Active_evaluated(max_size+3,i) = num2cell(sol(i,1)/(model.A(end,:)*soli.x));
end

save(strcat(path_save,'/iMMs_evaluated.mat'),Active_evaluated)
%% essentiality analysis for each one of the iMMs
indNF = getAllVar(ttmodel,{'NF'});
essThr = 0.1;
flagTasks = 0;
model = iNTS_SL1344_v1_0;

%close all upatkes
idx_ex = find(contains(model.varNames,'R_EX_'));
model.var_ub(idx_ex) = 0;

for i = 1:length(A_Active)

    modeli = model;

    alt = modelmilp.varNames(modelmilp.indUSE(find(A_Active(:,i)==-1)));
    iMM_R = strrep(alt,'BFUSE_','');
    f1 = find(ismember(modeli.varNames,iMM_R));
    modeli.var_ub(f1) = 25;

    idx_his = find(ismember(modeli.varNames,'R_EX_his__L_e'));
    modeli.var_ub(idx_his) = 0.02;

    solWT = solveTFAmodelCplex(modeli);

    [grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(modeli, 'TFA', ttmodel.genes, 0, 0, 0, essThr, indNF);
    if any(isnan(grRateKO_genetfa)) %keep track of the NaN KO by comparing essential_genes_tfaNaN with essential_genes_tfa
        grRateKO_genetfaNaN = grRateKO_genetfa(:,i);
        essential_genes_tfaNaN = ttmodel.genes(grRateKO_genetfaNaN(:,1)<essThr*solWT.val);
        yesgenetfaNaN = ismember(ttmodel.genes,essential_genes_tfaNaN);
    end

    grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0; %by default NaN is considered an essential KO
    essential_genes_tfa = ttmodel.genes(grRateKO_genetfa<essThr*solWT.val);
    yesgenetfa = ismember(ttmodel.genes,essential_genes_tfa);

    save(strcat(path_save,'/essentiality/iMM_',num2str(i),'_essetiality.mat'),'essential_genes_tfa','grRatio_genetfa','grRateKO_genetfa')
end
%% minimal products for each iMM
model = iNTS_SL1344_v1_0;
obj = find(model.f);

% close all uptakes
idx_ex = find(contains(model.varNames,'R_EX_'));
model.var_ub(idx_ex) = 0;

% constraint carbon and nitrogen uptakes
% get carbon sources
f_ex = getCarbonSources(model,'tfa');
exchange = strrep(model.varNames(f_ex),'R_EX_','EX_');
for i = 1:length(exchange)
    f1 = find(ismember(model.rxns,exchange{i}));
    metsIdx(i,1) = find(model.S(:,f1));
end

%get carbon atoms
CarbonAtoms = getCarbonAtoms(model,metsIdx);
%
[NumCons, NumVars] = size(model.A);
NewCons = zeros(NumVars,1);
NewCons(f_ex) = CarbonAtoms;
model.A(NumCons+1,:) = NewCons;
model.rhs(NumCons+1) = 60; %10 mmol/gDW h of glucose
model.constraintNames{NumCons+1} = ['CUT_pUP_carbon'];model.constraintType{NumCons+1} = '<';

%get nitrogen sources
f_ex = getNitrogenSources(model,'tfa');
exchange = strrep(model.varNames(f_ex),'R_EX_','EX_');
for i = 1:length(exchange)
    f1 = find(ismember(model.rxns,exchange{i}));
    metsIdxN(i,1) = find(model.S(:,f1));
end

%get nitrogen atoms
NitrogenAtoms = getNitrogenAtoms(model,metsIdxN);

[NumCons, NumVars] = size(model.A);
NewCons = zeros(NumVars,1);
NewCons(f_ex) = NitrogenAtoms;
model.A(NumCons+1,:) = NewCons;
model.rhs(NumCons+1) =1000; %1 mmol nh4 / 5 mmol glucose
model.constraintNames{NumCons+1} = ['CUT_pU_Pnitrogen'];
model.constraintType{NumCons+1} = '<';


for i = 1:length(A_Active)

    modeli = model;
    % get each iMM solutions
    alt = modelmilp.varNames(modelmilp.indUSE(find(A_Active(:,i)==-1)));

    iMM_R = strrep(alt,'BFUSE_','');
    A(1:length(alt),i) = iMM_R;

    f1 = find(ismember(modeli.varNames,iMM_R));
    modeli.var_ub(f1) = 25;
    idx_his = find(ismember(modeli.varNames,'R_EX_his__L_e'));
    modeli.var_ub(idx_his) = 0.02;

    soli = solveTFAmodelCplex(modeli);
    % identify products at optimal growth
    modeli.var_lb(obj) = soli.val*0.99;


    path_save = strcat(path_save,'/iMM_',num2str(i),'/');

    % Define milp problem
    modelmilpIMS = analysisIME(modeli,0,soli.val*0.99, soli.val, {},[],600);
    filename = strcat(modeldescription,'_PhenoMappingSecretions');
    % Get alternative solutions
    [DPsimm,modelmilpIMS] = findDPMax_IME(modelmilpIMS, 10000, modelmilpIMS.indUSE,600, 1, strcat(path_save,filename));
    save(strcat(path_save,filename,'.mat'));

    % test io
    filename = strcat(modeldescription,'_PhenoMappingSecretions');
    r1 = load(strcat(path_save,filename,'.mat'));

    %Extract info about composition of the IMMs
    [immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
        immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfoIMMDPs_IME(r1.modelmilpIMS, r1.DPsimm, ...
        r1.modelmilpIMS.indUSE,  r1.modelmilpIMS.rowDrain_F,...
        r1.modelmilpIMS.rowDrain_R,r1.modelmilpIMS.F_Flux,  r1.modelmilpIMS.R_Flux) ;


    immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
        num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];
    filename = strcat(modeldescription,'iMM_',num2str(i),'_PhenoMappingSecretions');

    writeData(strcat(path_save,filename,'.csv'), immOutput.summary,...
        '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
        'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
        '%s\t%s\t%s\t%s\t%s');
    save(strcat(path_save,'Active.mat'),'Active');
    %
    elimIntermPMFile(path_save,filename)
end

