%% sensitivity analysis
clear
clc
load('./models/iNTS_SL1344_v1_0.mat')
model = iNTS_SL1344_v1_0;
load('/Users/eva/Documents/MATLAB/GIT/iNTS_SL1344/matFiles/exchange.mat')
f = find(ismember(model.rxns,exchange));
model.lb(f) = -25;

% Original coefficients for the biomass equation
obj = find(model.c);
idx = find(model.S(:,obj));
original_coefficients = model.S(idx,obj)'; 
% Define the number of perturbations
num_perturbations = 100;

% Preallocate a matrix to store the perturbed coefficients
perturbed_coefficients = zeros(num_perturbations, length(original_coefficients));

% Generate the perturbed coefficients
% Generate the perturbed coefficients
for i = 1:num_perturbations
    % Generate a random perturbation factor between 0.8 and 1.2 for each coefficient
    perturbation_factors = 0.8 + (1.2 - 0.8) * rand(1, length(original_coefficients));
    
    % Apply the perturbation to each coefficient while preserving its sign
    perturbed_coefficients(i, :) = abs(original_coefficients) .* perturbation_factors .* sign(original_coefficients);
end
    
% Original lower bound of the ATPM reaction
atpm = find(ismember(model.rxns,'ATPM'));
original_lower_bound = model.lb(atpm); % Replace with the actual lower bound value

% Define the number of perturbations
num_perturbations = 100;

% Preallocate a vector to store the perturbed lower bounds
perturbed_lower_bounds = zeros(num_perturbations, 1);

% Generate the perturbed lower bounds
for i = 1:num_perturbations
    % Generate a random perturbation factor between 0.8 and 1.2
    perturbation_factor = 0.8 + (1.2 - 0.8) * rand();
    
    % Apply the perturbation to the lower bound
    perturbed_lower_bounds(i) = original_lower_bound * perturbation_factor;
end

% for each new biomass composition, GAM and NGAM convert the model to
% thermo
load('./matFiles/thermo_data.mat')
mkdir('sensitivity_analysis')
for i = 1:length(perturbed_coefficients)
modeli = model;
modeli.S(idx,obj) = perturbed_coefficients(i,:)';
modeli.lb(atpm) =  perturbed_lower_bounds(i);
modeli.ub(atpm) =  perturbed_lower_bounds(i);
sol = solveFBAmodelCplex(modeli);
soli(i,1) = sol.f;

tmodel = prepModelforTFA(modeli,thermo_data,model.CompartmentData);
ttmodel = convToTFA(tmodel, thermo_data, {'AIRC3';'GARFT';'UDCPPtppi'}, 'DGo', {}, 0.001, 1, 1);
%save the new models
save(strcat('model',num2str(i),'.mat'),'ttmodel')
end

