function [r2_fm_vm_ss,r2_fm_rm_ss,r2_vm_rm_ss,beta] = ptg_glm_model(runs,ptg_values,ptg_variables,ks,on,sm,spk,sav,dist)
% Ptg Paper - adamhsugi@gmail.com
% Da Cunha's Lab - Dpt of Pharmacology - Universidade Federal do Parana 
% Curitiba, Brazil
%
%   GLM model for variables of ptg project
%   with k-fold cross validation
%   Inputs:
%       - runs = vector with the same length of ptg_values tagging each
%       value to a specific run (trial/reward approach)
%       - ptg_values = time series array with all the neural data of individual neurons + behavioral data for the glm
%       regression (e.g. spk/bin - spd - acc...)
%       - ptg_variables = name of all the variables included in the model
%       - ks = number of "folds" of the cross validation (default = 5)
%       - on = plots are made if on is equal to 1
%       - sm = smoothing window (moving average smoothing window)
%       - spk = index for figures generating
%       - sav = save on
%       - dist = glm distribution type (e.g. poisson or normal)
%   Outputs:
%       - r2_fm_vm_ss = r-squared (variance explained) by the comparing
%       full model with var model (fullmodel - var)
%       - r2_fm_rm_ss = r-squared (variance explained) by the comparing
%       full model with the real data
%       - r2_vm_rm_ss = r-squared (variance explained) by the comparing
%       var model with the real data
%
%   ptg_values generated by regression_data_ptg_arms.m script or just
%   generate an table with columns as variables (1st column as firing rate)

% variables to store the calculated r-squared for each models
r2_fm_vm_ss = []; 
r2_fm_rm_ss = [];
r2_vm_rm_ss = [];
beta = [];
% ptg_variables = {'spks_bin', 'speeds', 'accs', 'fracs_run',...
%                  'arms2','arms3','run_durations', 'trials',...
%                  'angles','h_dirs', 'r_l_mov','random_uni',...
%                  'random_norm'};

% Separating data for fivefold cross-validation
if nargin <2
   ks = 5;
   on = 1;
elseif nargin <3
   on = 1;
elseif nargin <4
    sm = 3;
elseif nargin <5
    sav = 0;
elseif nargin <6
    dist = 'poisson';
end

if strcmp('normal',dist)
    lk = 'identity';
else
    lk = 'log';
end

%%%%%%%%%%%%% Five-fold data sets %%%%%%%%%%%%%%%%%%%
% Settings for random selection of trials for training and test datasets
nruns = unique(runs); % simplify the runs vector and check trials included in the analysis
num_runs = size(nruns,1); % number of trials
runs_perdataset = floor(num_runs/ks); % maximum number of trials in each data set to balance analysis
perm = randperm(num_runs); % generating random tags for trials
combinations = []; % to store combinations of trials
for k = 1:ks
    combinations = [combinations; perm((k-1)*runs_perdataset+1:(k)*runs_perdataset)]; % storing combinations of trials for each dataset
end

% Separating runs data
runs_data = [];
runs_data{num_runs} = [];
for r = 1:num_runs
    runs_data{nruns(r)} = ptg_values(find(runs == nruns(r)),:); % using combinations to separate neural and behavioral data of each dataset
end
runs_data = runs_data(~cellfun('isempty',runs_data));

% Creating data sets
datasets = [];
datasets{ks} = []; %storing datasets
for k = 1:ks    
    dats = runs_data(combinations(k,:));
    data = [];
    for d = 1:size(dats,2)
        data = [data; dats{d}];
    end
    datasets{k} = data; % separating datasets
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating models - one model for each "fold" of the cross validation AND
% one for each variable
var = 1;
ks_comb = [1, 2, 3, 4, 5;...
           2, 3, 4, 5, 1;...
           3, 4, 5, 1, 2;...
           4, 5, 1, 2, 3;...
           5, 1, 2, 3, 4]; % combination of training and test data sets  (e.g. [tr tr tr tr test])

var_nums = 2:size(ptg_variables,2); %number of the variables
for v = 2:size(ptg_variables,2)
    r2_fm_vm_s = []; % storing R-squared for full model vs var model
    r2_fm_rm_s = []; % storing R-squared for full model vs real data 
    r2_vm_rm_s = []; % storing R-squared for var model vs real data
    
    v_ex = v; % var evaluated
%     v_ex = [10];
    var_nums = 2:size(ptg_variables,2); %number of the variables
    for e = 1:size(v_ex,2)
        var_nums(find(var_nums == v_ex(e))) = []; % evaluated variable (excluding from the full model)
    end
    
    contrs = []; % for storing contributions for each variable       
    for k = 1:ks % one R-squared for each validation
        % full model training and test dataset
        fmdata_tr = [datasets{ks_comb(k,1)}; datasets{ks_comb(k,2)}; datasets{ks_comb(k,3)}; datasets{ks_comb(k,4)}];
        fmdata_resp = fmdata_tr(:,1); % responses - spks/bin
        fmdata_tr = fmdata_tr(:,2:end); % full model variables data
        fmdata_te = datasets{ks_comb(k,5)};
        fmdata_te_data = fmdata_te(:,2:end);
                        
        % var model training and test dataset (var model = full model - variable been evaluated)
        vardata_tr = [datasets{ks_comb(k,1)}; datasets{ks_comb(k,2)}; datasets{ks_comb(k,3)}; datasets{ks_comb(k,4)}];
        vardata_tr = vardata_tr(:,var_nums); % var model variables data       
        vardata_resp = fmdata_resp; % responses - spks/bin        
        vardata_te = datasets{ks_comb(k,5)}; 
        vardata_te_data = vardata_te(:,var_nums); % test data        
         
        real_resp = fmdata_te(:,1); % real data
        
        % fitting the models
        [fm_bs,~,fmstats] = glmfit(fmdata_tr, fmdata_resp,dist,'constant','on');
        fm_resp = glmval(fm_bs,fmdata_te_data,lk);
        beta = [beta fm_bs];
        
        [var_bs,~,varstats] = glmfit(vardata_tr, vardata_resp,dist,'constant','on');
        var_resp = glmval(var_bs,vardata_te_data,lk);
        
        % Plotting data
        % showing full model and var model results
        if on == 1
            if v < ((size(ptg_variables,2)-1)/2)+2
                fig = figure(100+spk);
                sp = subplot(((size(ptg_variables,2)-1)/2),5,(((v-2)*5)+k));
            else
                fig = figure(1+100+spk);
                sp = subplot(((size(ptg_variables,2)-1)/2),5,(((v-2)*5)+k)-30);
            end                
            
            % title for each run
            
            if ((v-2)*5)+k < 6
                title(sprintf('Run %d',((v-2)*5)+k))
            elseif ((v-2)*5)+k < 36 & ((v-2)*5)+k > 30
                 title(sprintf('Run %d',((v-2)*5)+k-30))
            end
            fig.Position = [200 100 700 600];
            hold on
            plot(fm_resp(1:60),'Color',[.8 .2 .2],'LineStyle','-','LineWidth',2)
            plot(smooth(real_resp(1:60),sm),'Color','k','LineWidth',2)
            plot(var_resp(1:60),'Color',[.4 .6 .8],'LineWidth',3)
            xlim([0 60])
            if (((v-2)*5)+k) == 60 
                a = gca;
                yyaxis left
                ya = gca;
                ya.YAxis(1).Visible = 'off';                                
                ya.YAxis(2).Visible = 'on';
                yyaxis right
                ya.YColor = [0.1500    0.1500    0.1500];
                ylabel('Spks/bin')                
                xlabel('Time (s)')
                xticks([0 60])
                xticklabels([0 6])
                set(a,'FontSize',15)
            else
                a = gca;                
                axis off
%                 ylabel(sprintf('%s',ptg_variables{v}))
            end            
        end        
        % calculating R-squared for both models
        r2_fm_vm = corr(fm_resp,var_resp).^2;
        r2_fm_rm = corr(fm_resp,real_resp).^2;
        r2_vm_rm = corr(var_resp,real_resp).^2;
        
        r2_fm_vm_s = [r2_fm_vm_s r2_fm_vm];
        r2_fm_rm_s = [r2_fm_rm_s r2_fm_rm];
        r2_vm_rm_s = [r2_vm_rm_s r2_vm_rm];
        
    end
    
    % storing R-squared for all models and dataset combinations
    r2_fm_vm_ss = [r2_fm_vm_ss; r2_fm_vm_s];
    r2_fm_rm_ss = [r2_fm_rm_ss; r2_fm_rm_s];
    r2_vm_rm_ss = [r2_vm_rm_ss; r2_vm_rm_s];
end
% saving figures
if sav == 1
    set(gcf,'PaperPosition', [0 0 15 10],'Renderer','Painters'); %Position plot at left hand corner with width 5 and height 5.
    set(gcf,'PaperSize', [15 10]);
    print(gcf,sprintf('GLM_model %d',spk),'-dpdf','-r300')
end