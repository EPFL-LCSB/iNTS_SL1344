%% script to test the model performance against Biolog Data
% go to the iNTS_SL1344 folder
cd('/PathTo/iNTS_SL1344/')

% add paths
addpath(genpath('./'));
addpath(genpath('/PathTo/mattfa/'));

%load files
load('./matFiles/Biolog_IDs.mat')
load('./models/iNTS_SL1344_v1_0.mat')

%% test carbon sources (PM1) - FBA

model = iNTS_SL1344_v1_0;

aerobicTag = 1; % aerobic conditions 1; anaerobic 0

% set media
model.lb(find(ismember(model.rxns,'EX_glc__D_e'))) = 0;
model.lb(find(ismember(model.rxns,'EX_nh4_e'))) = -10;


if aerobicTag
    model.lb(find(ismember(model.rxns,'EX_o2_e'))) = -20;
else
    model.lb(find(ismember(model.rxns,'EX_o2_e'))) = 0;
end


for i = 1:length(pm1)
    modeli = model;
    % for each compound in pm1 identify the corresponding exchange
    f = find(ismember(model.mets,strcat(pm1{i},'_e')));
    if f

        f1 = find(model.S(f,:));
        for k = 1:length(f1)
            s = find(model.S(:,f1(k)));
            if length(s)==1
                lb_ex = model.rxns(f1(k));
                break
            end
        end


        fi = find(ismember(modeli.rxns,lb_ex));
        modeli.lb(fi) = -25;
        % make sure the compound is uptaken
        modeli.ub(fi) = -10^-6;

        sol = solveFBAmodelCplex(modeli);

        if not(isempty(sol.f))
            growth_fba(i,1) = sol.f;
        else
            growth_fba(i,1) = 0;
        end
    else
        growth_fba(i,1) = 0;
    end

end

%% test carbon sources (PM1) - TFA

ttmodel = iNTS_SL1344_v1_0;

aerobicTag = 1; % aerobic conditions 1; anaerobic 0

% set media
ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_glc__D_e'))) = 0;
ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_nh4_e'))) = 10;

if aerobicTag
    ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_o2_e'))) = 20;
else
    ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_o2_e'))) = 0;
end

for i = 1:length(pm1)

    modeli = ttmodel;

    % impose lower bound on growth rate to solve faster (optional)
    fj = find(ismember(modeli.varNames,strcat('F_',ttmodel.rxns(find(ttmodel.c)))));
    modeli.var_lb(fj) = 10^-6;

    % for each compound in pm1 identify the corresponding exchange
    f = find(ismember(model.mets,strcat(pm1{i},'_e')));
    if f

        f1 = find(model.S(f,:));
        for k = 1:length(f1)
            s = find(model.S(:,f1(k)));
            if length(s)==1
                lb_ex = model.rxns(f1(k));
                break
            end
        end

        fi = find(ismember(modeli.varNames,strcat('R_',lb_ex)));
        modeli.var_ub(fi) = 25;
        % make sure the compound is uptaken
        modeli.var_lb(fi) = 10^-6;
        % "open" net flux varibale
        fi = find(ismember(modeli.varNames,strcat('NF_',lb_ex)));
        modeli.var_lb(fi) = - 25;

        sol = solveTFAmodelCplex(modeli);

        if not(isempty(sol.val))
            growth_tfa(i,1) = sol.val;
        else
            growth_tfa(i,1) = 0;
        end
    else
        growth_tfa(i,1) = 0;
    end
end


%% test nitrogen sources (PM2) - FBA
model = iNTS_SL1344_v1_0;

aerobicTag = 1; % aerobic conditions 1; anaerobic 0

% set media
model.lb(find(ismember(model.rxns,'EX_glc__D_e'))) = 0;
model.lb(find(ismember(model.rxns,'EX_pyr_e'))) = -25;
model.lb(find(ismember(model.rxns,'EX_nh4_e'))) = 0;


if aerobicTag
    model.lb(find(ismember(model.rxns,'EX_o2_e'))) = -20;
else
    model.lb(find(ismember(model.rxns,'EX_o2_e'))) = 0;
end

for i = 1:length(pm2)
    modeli = model;
    if not(ismember(pm2{i},dipeptide))
        f = find(ismember(model.mets,strcat(pm2{i},'_e')));

        if f
            f1 = find(model.S(f,:));
            for k = 1:length(f1)
                s = find(model.S(:,f1(k)));
                if length(s)==1
                    lb_ex = model.rxns(f1(k));
                    break
                end
            end
            fi = find(ismember(modeli.rxns,lb_ex));
            modeli.lb(fi) = -5;
            % make sure the compound is uptaken
            modeli.ub(fi) = -10^-6;
            sol = solveFBAmodelCplex(modeli);
            if not(isempty(sol.f))
                growth_n_fba(i,1) = sol.f;
            else
                growth_n_fba(i,1) = 0;
            end
        else
            growth_n_fba(i,1) = 0;
        end
    else
        fi = find(ismember(modeli.rxns,strcat('EX_',pm2{i},'_e')));
        modeli.lb(fi) = -5;
        sol = solveFBAmodelCplex(modeli);
        if not(isempty(sol.f))
            growth_n_fba(i,1) = sol.f;
        else
            growth_n_fba(i,1) = 0;
        end
    end
end

%% test nitrogen sources (PM2) - TFA
ttmodel = iNTS_SL1344_v1_0;

aerobicTag = 1; % aerobic conditions 1; anaerobic 0

% set media
ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_glc__D_e'))) = 0;
ttmodel.var_lb(find(ismember(ttmodel.varNames,'NF_pyr_e'))) = -25;
ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_pyr_e'))) = 25;
ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_nh4_e'))) = 0;

if aerobicTag
    ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_o2_e'))) = 20;
else
    ttmodel.var_ub(find(ismember(ttmodel.varNames,'R_EX_o2_e'))) = 0;
end

for i = 1:length(pm2)
    modeli = ttmodel;

    if not(ismember(pm2{i},dipeptide))
        f = find(ismember(model.mets,strcat(pm2{i},'_e')));

        if f
            f1 = find(model.S(f,:));
            for k = 1:length(f1)
                s = find(model.S(:,f1(k)));
                if length(s)==1
                    lb_ex = model.rxns(f1(k));
                    break
                end
            end
            fi = find(ismember(modeli.varNames,strcat('R_',lb_ex)));
            modeli.var_ub(fi) = 5;
            modeli.var_lb(fi) = 10^-6;
            fi = find(ismember(modeli.varNames,strcat('NF_',lb_ex)));
            modeli.var_lb(fi) = -5;
            sol = solveTFAmodelCplex(modeli);
            if not(isempty(sol.val))&& sol.val>0
                growth_n_tfa(i,1) = sol.val;
            else
                growth_n_tfa(i,1) = 0;
            end
        else
            growth_n_tfa(i,1) = 0;
        end
    else
        fi = find(ismember(modeli.varNames,strcat('R_EX_',pm2{i},'_e')));
        modeli.var_ub(fi) = 5;
        modeli.var_lb(fi) = 10^-6;
        fi = find(ismember(modeli.varNames,strcat('NF_EX_',pm2{i},'_e')));
        modeli.var_lb(fi) = -5;
        sol = solveTFAmodelCplex(modeli);
        if not(isempty(sol.val))&& sol.val>0

            growth_n_tfa(i,1) = sol.val;
        else
            growth_n_tfa(i,1) = 0;
        end
    end
end

