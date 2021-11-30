% Counterbalanced fivefold cross validation
%
% Ptg Paper - gabriell.baltazar99@gmail.com
% Da Cunha's Lab - Dpt of Pharmacology - Universidade Federal do ParanÃ¡
% Curitiba, Brazil
%
% Script to compute five (for fivefold cross-validation of GLM model) or
% two (for twofold cross-validation of 2 Gaussians model) data pools that
% are used in further validation methods
%

[filename, pathname] = uigetfile('*.mat', 'Select mat data files', 'Multiselect', 'on');
if ~iscell(filename)
    filename={filename};
end

y=1;
fit=[];
fit{5}=[];
reg=[];
reg{5}=[];
infos=[];
infos{66,5}=[];

for a=1:size(filename,2)  %carrega arquivos  
    
    data=load(sprintf('%s%s', pathname, filename{a}));     %load the files selected in the beggining of the script one at a time    
    
    fields=fieldnames(data);
    fields1=strncmp(fields,'SPK',3);
    fields2=find(fields1);        
    
    for b=1:length(fields2)  %separa quantos neurônios tem no arquivo pra trabalhar com eles separados
        
        e5=strcat('data.',(fields{fields2(b)}));                
        e6=eval(e5);                
        
        animal=filename{a}(1:end-4);
        neuron=fields{fields2(b)};
        name=[animal ' ' neuron];
        name1={codes(sprintf('%s',name))};                                                  
        
        runs=[];
        
        for c=1:length(data.cutteddata)
            
            if c~=1 && c~=length(data.cutteddata) && data.cutteddata(c).tag(1)==1
                if data.cutteddata(c).tsst(1)-1.050>data.cutteddata(c-1).tsend(1) && data.cutteddata(c).tsend(1)+1.050<data.cutteddata(c+1).tsst(1)
                    if data.cutteddata(c).efc(1)>0.8 && data.cutteddata(c).efc<1.2
                        runs=[runs;c];
                    end
                end                
            end
            
        end
        
        for d=1:length(runs)
            
            d1=runs(d);
            bin=0.150;
            
            bef=find(e6>=data.cutteddata(d1).tsst(1)-1.050 & e6<data.cutteddata(d1).tsst(1));                                
            beflast=data.cutteddata(d1).tsst(1)-1.050;
            bef1=e6(bef);

            acc=find(e6>=data.cutteddata(d1).tsst(1) & e6<=data.cutteddata(d1).tspk(1));
            accb=(data.cutteddata(d1).tspk(1)-data.cutteddata(d1).tsst(1))/10;
            acclast=data.cutteddata(d1).tsst(1);
            acc1=e6(acc);

            dea=find(e6>data.cutteddata(d1).tspk(1) & e6<=data.cutteddata(d1).tsend(1));
            deab=(data.cutteddata(d1).tsend(1)-data.cutteddata(d1).tspk(1))/10;
            dealast=data.cutteddata(d1).tspk(1);
            dea1=e6(dea);

            aft=find(e6>data.cutteddata(d1).tsend(1) & e6<=data.cutteddata(d1).tsend(1)+1.050);
            aftlast=data.cutteddata(d1).tsend(1);
            aft1=e6(aft);
            
            for e=1:7                
                beflast1=beflast+bin;                                                                                    
                befbin=bef1>beflast & bef1<=beflast1;                
                if ~isempty(befbin)
                    befbin1=find(befbin);
                    befbinsfr(d,e)=length(befbin1)/bin;                    
                else
                    befbinsfr(d,e)=0;
                end                                              
                beflast=beflast1;                                                                
                
                aftlast1=aftlast+bin;                        
                aftbin=aft1>aftlast & aft1<=aftlast1;
                if ~isempty(aftbin)
                    aftbin1=find(aftbin);
                    aftbinsfr(d,e)=length(aftbin1)/bin;
                else
                    aftbinsfr(d,e)=0;
                end                                    
                aftlast=aftlast1;
                
                clear('befbin','aftbin');                
            end
            
            for e1=1:10                
                acclast1=acclast+accb;                                                                                   
                accbin=acc1>acclast & acc1<=acclast1;                
                if ~isempty(accbin)
                    accbin1=find(accbin);
                    accbinsfr(d,e1)=length(accbin1)/accb;                        
                else
                    accbinsfr(d,e1)=0;
                end                                                        
                acclast=acclast1;
                
                dealast1=dealast+deab;                                                                                    
                deabin=dea1>dealast & dea1<=dealast1;              
                if ~isempty(deabin)
                    deabin1=find(deabin);
                    deabinsfr(d,e1)=length(deabin1)/deab;                                        
                else
                    deabinsfr(d,e1)=0;
                end                                    
                dealast=dealast1;
                
                clear('accbin','deabin');                
            end            
            clear('bef','acc','dea','aft','bef1','acc1','dea1','aft1','beflast','acclast','dealast','aftlast','bin','accb','deab');                        
        end
        
        totalbinsfr(1:length(runs),1:7)=befbinsfr;
        totalbinsfr(1:length(runs),8:17)=accbinsfr;
        totalbinsfr(1:length(runs),18:27)=deabinsfr;
        totalbinsfr(1:length(runs),28:34)=aftbinsfr;
        
        totalz=zscore(totalbinsfr,0,'all');                
        
        rwd_arms=rwds(name(1:3));                
        
        ks=5;
        nruns=unique(runs);
        num_runs = size(nruns,1);
        runs_perdataset = floor(num_runs/ks);
        n(y)=runs_perdataset;
        combinations = [];
        arml(1:5,1)=0;
        armm(1:5,1)=0;
        armh(1:5,1)=0;
        j1=1;
        
        for i=1:runs_perdataset
            point(1:5,1)=1;
            for j=1:ks
                
                test_arm=data.cutteddata(runs(j1)).arm(1);                
                
                if i==1
                    if test_arm==rwd_arms(1)
                        arml(j,1)=arml(j,1)+1;
                    elseif test_arm==rwd_arms(2)
                        armm(j,1)=armm(j,1)+1;
                    elseif test_arm==rwd_arms(3)
                        armh(j,1)=armh(j,1)+1;
                    end
                    
                    combinations(j,i)=j1;
                else
                                                                                
                    if test_arm==rwd_arms(1)
                        arm_test=arml;
                    elseif test_arm==rwd_arms(2)
                        arm_test=armm;
                    elseif test_arm==rwd_arms(3)
                        arm_test=armh;
                    end
                    
                    [mini,ind]=min(arm_test);
                    
                    if point(ind)==1
                        combinations(ind,i)=j1;
                        point(ind)=0;
                        if test_arm==rwd_arms(1)
                            arml(ind,1)=arml(ind,1)+1;
                        elseif test_arm==rwd_arms(2)
                            armm(ind,1)=armm(ind,1)+1;
                        elseif test_arm==rwd_arms(3)
                            armh(ind,1)=armh(ind,1)+1;
                        end
                    else
                        arm_test=1;
                        ind1=ind;
                        while arm_test==1
                            ind1=ind1+1;
                            [mini,ind]=min(arm_test(ind1:end));
                            
                            if point(ind)==1
                                arm_test=0;
                                combinations(ind,i)=j1;
                                point(ind)=0;
                                if test_arm==rwd_arms(1)
                                    arml(ind,1)=arml(ind,1)+1;
                                elseif test_arm==rwd_arms(2)
                                    armm(ind,1)=armm(ind,1)+1;
                                elseif test_arm==rwd_arms(3)
                                    armh(ind,1)=armh(ind,1)+1;
                                end                                                                                           
                            end
                            
                            if ind1>=5
                                arm_test=0;
                                test_arm2=find(point);
                                combinations(test_arm2(1),i)=j1;
                                point(test_arm2(1))=0;
                                
                                if test_arm==rwd_arms(1)
                                    arml(test_arm2(1),1)=arml(test_arm2(1),1)+1;
                                elseif test_arm==rwd_arms(2)
                                    armm(test_arm2(1),1)=armm(test_arm2(1),1)+1;
                                elseif test_arm==rwd_arms(3)
                                    armh(test_arm2(1),1)=armh(test_arm2(1),1)+1;
                                end
                                clear('test_arm2')
                            end
                            
                            
                        end
                    end
                    
                end
                
                j1=j1+1;
            end
        end
        clear('arml','armm','armh','point');
        
        for t=1:size(combinations,1)
            infos{y,t}(1:15)=0;
            for t1=1:size(combinations,2)
                
                test_arm=data.cutteddata(runs(combinations(t,t1))).arm(1);
                if test_arm==rwd_arms(1)
                    infos{y,t}(1)=infos{y,t}(1)+1;
                elseif test_arm==rwd_arms(2)
                    infos{y,t}(2)=infos{y,t}(2)+1;
                elseif test_arm==rwd_arms(3)
                    infos{y,t}(3)=infos{y,t}(3)+1;                    
                end
                
                test_trial=data.cutteddata(runs(combinations(t,t1))).trial(1);
                infos{y,t}(test_trial+3)=infos{y,t}(test_trial+3)+1;
                
                clear('test_arm','test_trial');
                
            end
        end
        clear('rwd_arms');
        
        % Separating runs data
        runs_data = [];
        runs_data{num_runs} = [];
        for r = 1:num_runs
            runs_data{nruns(r)} = totalz(r,:);
        end
        runs_data = runs_data(~cellfun('isempty',runs_data));
        
        % Creating data sets
        datasets = [];
        datasets{ks} = []; %storing datasets
        for k = 1:ks    
            dats = runs_data(combinations(k,:));
            data1 = [];
            for d = 1:size(dats,2)
                data1 = [data1; dats{d}];
            end
            datasets{k} = data1; % separating datasets
            clear('dats','data1');
        end
        
        combs=[1, 2, 3, 4, 5;...
               2, 3, 4, 5, 1;...
               3, 4, 5, 1, 2;...
               4, 5, 1, 2, 3;...
               5, 1, 2, 3, 4];
                
        for f1=1:size(combs,1)                       
            for f2=1:size(combs,2)                    
                if f2~=length(combs)
                    test((f2-1)*runs_perdataset+1:(f2)*runs_perdataset,1:34)=datasets{combs(f1,f2)};
                else
                    test1=mean(test,1);
                    test2=mean(datasets{combs(f1,f2)},1);
                end
            end
            fit{f1}(y,:)=test1;
            reg{f1}(y,:)=test2;
            clear('test','test1','test2')
        end
        y=y+1;
        clear('nruns','num_runs','runs_perdataset','perm','combinations','runs_data','dataset','combs');
        clear('befbinsfr','accbinsfr','deabinsfr','aftbinsfr','totalbinsfr','totalz','e5','e6','runs');
    end    
end

data_info=[];
data_info{5}=[];

for w=1:size(infos,2)
    for w1=1:size(infos,1)
        data_info{w}(w1,1:15)=infos{w1,w};
    end
end

info1=data_info{1};
info2=data_info{2};
info3=data_info{3};
info4=data_info{4};
info5=data_info{5};

poss1=fit{1};
poss1(1:66,35:37)=NaN;
poss1(1:66,38:71)=reg{1};

poss2=fit{2};
poss2(1:66,35:37)=NaN;
poss2(1:66,38:71)=reg{2};

poss3=fit{3};
poss3(1:66,35:37)=NaN;
poss3(1:66,38:71)=reg{3};

poss4=fit{4};
poss4(1:66,35:37)=NaN;
poss4(1:66,38:71)=reg{4};

poss5=fit{5};
poss5(1:66,35:37)=NaN;
poss5(1:66,38:71)=reg{5};

arc_name='Fivefold(3).xlsx';

xlswrite(arc_name,poss1,'Combination1');
xlswrite(arc_name,poss2,'Combination2');
xlswrite(arc_name,poss3,'Combination3');
xlswrite(arc_name,poss4,'Combination4');
xlswrite(arc_name,poss5,'Combination5');
xlswrite(arc_name,info1,'Informations dataset 1');
xlswrite(arc_name,info2,'Informations dataset 2');
xlswrite(arc_name,info3,'Informations dataset 3');
xlswrite(arc_name,info4,'Informations dataset 4');
xlswrite(arc_name,info5,'Informations dataset 5');