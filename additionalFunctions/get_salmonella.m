
function get_salmonella
cd('./../')

paramSalmonella   = struct (...%model parameters
                        'Organism',                           'salmonella'                                ,... % human, ecoli, putida
                        'GEMname',                            'model4redGEMX'              ,... % name of the GEM used (as it is saved in GEMs folder)
                        'RedModelName',                       'redSalmonella'                      ,... % choose a name for the reduction
                        'SelectedSubsystems',                 {{'sourceModel'}}                       ,... % default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'x';'y';'z'}}
                        'AddETCAsSubsystem',                  'no'                                  ,... % yes, no
                        'AddExtracellularSubsystem',          'no'                                   ,... % no, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'DM_x';'DM_y';'DM_z'}}
                        'AerobicAnaerobic',                   'anaerobic'                              ,... % aerobic, anaerobic
                        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
                        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
                        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100      
                        'case_filename',                      'case_salmonella'                           ,... % name of the matlab function for the specific organism
                         ...%redGEM parameters
                        'L',                                  0                                    ,... % 1,2,...
                        'D',                                  0                                      ,... % 1,2,...
                        'startFromMin',                       'no'                                   ,... % yes, no
                        'ThrowErrorOnDViolation',             'error'                                ,... % error/continue If any two subsystems cannont connect to the level D                           
                        'OnlyConnectExclusiveMets',           'no'                                   ,... % yes, no
                        'ConnectIntracellularSubsystems',     'yes'                                  ,... % yes, no
                        'ApplyShortestDistanceOfSubsystems',  'bothways'                             ,... % bothways,eachdirection'
                        ...%redGEMX parameters
                        'performREDGEMX',                     'yes'                                  ,... % yes, no
                        'NumOfConnections',                   'SminMetE'                             ,... % OnePerMetE ,SminMetE
                        ...%lumpGEM parameters
                        'performLUMPGEM',                     'yes'                                  ,... % yes, no (do we want to perform lumping or not?)   
                        'PercentOfmuMaxForLumping',            100                                   ,... % Please specify the percentage of muMax that we should impose for the lumped reactions (100, 90, ..., 10, 0)
                        'addGAM',                             'yes'                                  ,... % Would you like to extract and add to the reduced model a growth associated maintenance (GAM) reaction? yes, no 
                        'PreventBBBuptake',                   'no'                                   ,... % yes/no: if yes, do not allow uptake through any bbb drains.     
                        'NumOfLumped',                        'OnePerBBB'                            ,... % OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
                        'AlignTransportsUsingMatFile',        'yesautomatic'                         ,... % yesusingmatfile, yesusingReducedmatfile, yesautomatic,no
                        'ImposeThermodynamics',               'yes'                                  ,... % would you like to impose thermodynamic cnstraints? yes, no
                        ...%postprocessing parameters
                        'performPostProcessing',              'yes'                                  ,... % yes, no
                        ...%solver parameters
                        'TimeLimitForSolver',                 'yes'                                  ,... % yes, no
                        'CplexParameters',                    'LCSBScaling'                          ,... % CplexDefault, LCSBDefault
                        ....%paths to folders
                        'CPLEX_PATH',                         'YourPathToCPLEX'                        ,... % provide path to CPLEX         
                        'TFA_PATH',                           'YourPath/mattfa'                       ,... % provide path to matTFA
                        'thermo_data_PATH',                   'YourPath/iNTS_SL1344/matFiles/thermo_data.mat'   ,... % provide path to thermodynamic data (including the file name)
                        'output_PATH',                        'YourPathToOutput/');                              % provide output path where the files will be saved
redGEM_salmonella(paramSalmonella);
end
