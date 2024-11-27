function [RedModel, activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, rxns] = redGEM_salmonella(varargin)
% this function is based on the redGEM function by Ataman et al. (2017), and Masid
% et al. (2020) 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTING UP THE OPTIONS OF THE REDUCTION  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If there is no input to the function, assign a string name to varargin
if isempty(varargin)
    varargin = {'NoInput'};
end

% Set all non-specified flags to empty. In this way, whatever has not been
% specified by the "varargin", will be later on asked with individual
% questions. All the variables that have been specified in "varargin" will
% overwrite these empty settings.
[paramRedGEM, paramNames] = redGEMQuestionnaire(varargin{1});


% Assign the value names to all the parameters of the current redGEM run
eval(cell2mat(cellfun(@(x) [x,'= getfield(paramRedGEM,''',x,''');'],paramNames,'UniformOutput', false)'));

if ~exist([output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname],'dir')
    mkdir([output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname])
end

if ~exist([output_PATH,'/UserOutputs/RunTimes'],'dir')
    mkdir([output_PATH,'/UserOutputs/RunTimes'])
end

connectingpaths_folder = [output_PATH,'/UserOutputs/ConnectingPaths'];
if ~exist(connectingpaths_folder,'dir')
    mkdir(connectingpaths_folder)
end

if ~exist([output_PATH,'/UserOutputs/Models/',Organism],'dir')
    mkdir([output_PATH,'/UserOutputs/Models/',Organism])
end

% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET NECESSARY THE PATHS AND CHECK THE VERSIONS USED  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore the default path
restoredefaultpath

% Add required paths:
% - Cobra and the tFBA path
addpath(genpath(TFA_PATH))

% - Cplex path
addpath(genpath(CPLEX_PATH))

% - redGEM path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(path))

if isempty(which('cplex.p'))
    error('You need to add CPLEX to the Matlab-path!!')
end

% Set the solver
changeCobraSolver('cplex_direct','LP')

%% %%%%%%%%%%%%%%%%%%%%%%
%% LOAD NECESSARY FILES %
%% %%%%%%%%%%%%%%%%%%%%%%

% Load the thermodynamic database used
load(thermo_data_PATH)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the properties of the GEM that is reduced %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['[OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem] = '...
    case_filename,'(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate);'])

%%
% Model Reaction Redundancy Check
redundant_table = ReacRedundancy(GEMmodel);

%% For the process of lumping, we can use the entire genomescale model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GSM_ForLumping = GEMmodel;

if strcmp(ImposeThermodynamics, 'no')
    % Remove unnecessary thermo fields (if they exist). If not, functions
    % like the extractSubnetwork might have an error.
    GSM_ForLumping = removeThermoFields(GSM_ForLumping);
end

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_1.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

%% Connecting Subsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we only want to perform the lumping on a predefined subsystem, or set
% of subsystems (i.e. we do not want to connect these subsystem, and/or
% include more reactions) then we can directly jump to the lumping part.
% The above implies that D=0, as we do not wish to expand the

% To find the adjacency matrices and perform all the graph search to obtain
% the desired reaction network we need a genomescale model without biomass,
% and without the metabolite pairs:
GSM_ForAdjMat = GEMmodel;
% Find and remove the biomass reactions
[~, BMIds] = ismember(Biomass_rxns, GSM_ForAdjMat.rxns);

GSM_ForAdjMat   = removeRxns(GSM_ForAdjMat, GSM_ForAdjMat.rxns(BMIds));

% We directly remove cofactors and small ionrganic metabolites from the S matrix
Smatrix = removeMetPairsFromS(GSM_ForAdjMat, met_pairs_to_remove, InorgMetSEEDIDs);
GSM_ForAdjMat.S = Smatrix;

rxns_ss = cell(length(core_ss),1);
% For each of the subsystems:
for i=1:length(core_ss)
    % - find the reactions of this subsystem, and their indices
    [~, rxns_ss{i}]=ismember(GSM_ForAdjMat.subSystems, core_ss{i});
    rxns_ss{i}=find(rxns_ss{i});
end
connecting_reactions = [];

all_core_rxns = unique([connecting_reactions;vertcat(rxns_ss{:})]);

submodel = extractSubNetwork(GEMmodel, GSM_ForAdjMat.rxns(all_core_rxns));
% We go through all GS-reactions, and check if there exists a reactions
% that involves ONLY core metabolites. If yes, we keep it, and add it as a
% core (i.e. tightening)
k = 1;
for i=1:length(GEMmodel.rxns)
    mets = GEMmodel.mets(find(GEMmodel.S(:, i)));
    [~, ba] = ismember(mets, submodel.mets);
    if isempty(find(ba==0))
        core_rxns(k,1) = i;
        k=k+1;
    end
end
clear k;


otherReactions = setdiff(1:length(GSM_ForLumping.rxns),core_rxns);

% Transform indices to GEM model (GSM_ForLumping) indexing
[~, otherReactionsGSMForLump_idx] = ismember(GSM_ForLumping.rxns(otherReactions), GSM_ForLumping.rxns);

fprintf('Connecting the metabolites from the extracellular medium to the core\n')
% connect extracellular subsystem to core
[ConnectExtrCell_rxns_all, ConnectExtrCell_id_all, sol_all] = ...
    redGEMX(rxns_ss, GSM_ForLumping, GSM_ForAdjMat, OriginalGEM, ...
    GEMmodel, Organism, GEMname, NumOfConnections, CplexParameters, DB_AlbertyUpdate, ImposeThermodynamics,output_PATH);

otherReactionsGSMForLump_idx = setdiff(otherReactionsGSMForLump_idx,unique(ConnectExtrCell_id_all));

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_5a.mat;'])   %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <


RedModel.info_Econnect = [OriginalGEM.ExtracellularMedium_connect ConnectExtrCell_rxns_all'];

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_7.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

% Store the core reduced model stoichiometry in the output (without the lumped reactions and without the transports)
RedModel.core_model = extractSubNetwork(GEMmodel, GEMmodel.rxns(core_rxns));

% > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
[dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
eval(['save ',output_PATH,'/TEMP/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_8.mat;'])    %
% < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <

RedModel.ReductionProperties = paramRedGEM;

% Store in the model the core/initial subsystems that are selected for the
% reduction
RedModel.initialSubSystems = core_ss;

% Add a description field. ATTENTION: this is very important for the ORACLE
% workflow later on!!
RedModel.description = [RedModelName,'_Date_',dateStr,'_',timeStr,'_Generated_by_redGEMv1'];

% workflow later on!!
RedModel.solverPath = CPLEX_PATH;

% Add a field with the name of the user that generated this model
user_str = regexp(pwd,'[\\/]Users[\\/](\w+)[\\/]','tokens');
user_str = user_str{1}{1};
RedModel.GeneratedByUser = user_str;

%% Saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the model in the output folder
eval([model_name,' = RedModel;']);
eval(['save ',output_PATH,'/UserOutputs/Models/',Organism,'/',model_name,'.mat ', model_name,';']);
end
