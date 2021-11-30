% New z-score plots for new alignments
%
% Ptg Paper - gabriell.baltazar99@gmail.com
% Da Cunha's Lab - Dpt of Pharmacology - Universidade Federal do Paraná
% Curitiba, Brazil
%
% Scritp to compute mean z-scores and plot it in heatmaps, aligned to five 
% events of the rats runs: locomotion onset, peak of acceleration, peak of 
% speed, peak of deceleration and locomotion end.
%
% Z-scores are computed using bins of 100 miliseconds, in windows raging
% from -3 seconds to +3 seconds according to the current event.
%

% Selecting files
[filename, pathname] = uigetfile('*.mat', 'Select mat data files', 'MultiSelect','on');
if ~iscell(filename)
    filename = {filename};
end

% setting parameters
window = [-3 3]; % window (xlim) of zscore plot
binsz = .1; % bins length
numbins = 60;   %number of bins
normon = 1; % normalization of time for runs
smoothon = 1;   % smoothing (1 applies, 0 not applies)
rgb = 'parula';% colormap_ptg(); %# colormap used
bootson = 0;
scron = 0; % scrambling data on

% excel file table writing
writet = 0; % writing excel table (1 writes, 0 not writes)
excelfilename = 'test_locoff';  % excel file name in case of writet=1

% figure parameters
classon = 0;    % prepares the heatmaps for specified neuron class
class = 'pv_test1'; % class analysed if classon=1 (possibilities: pv_test1; vp_test1; locoff_test1;
codeon = 1; % printing neuron codes on the figure (1 for coding, 0 for not coding)
axisfont = 10;  % axes font (10 for pdf, 15 for png figures)
titlefont = 8;  % titles font (8 for pdf, 15 for png figures)

% variables for storing neurons data
zs1 = [];
zs2 = [];
zs3 = [];
zss = [];
z_ss = [];
spkss = [];
runs_durations = [];
cruns_durations = [];
codess=[];

% setting neurons classes
pv={'0401','0402','0407','0408','0412','0413','0414','0415','0417','0418','0902','0904','0921','0925','1106','1109','1110'};
vp={'0410','0416','0804','0908','0911','0922','0926','1101','1107','1111'};
locoff={'0409','0801','0808','0906','0910','0912','0914','0915','0917','0918','1103'};

% Files loop
for al=1:5
    alig=al;
    for f = 1:size(filename,2)
        
        % Loading .mat files
        load(sprintf('%s%s',pathname,filename{f}));
        list = who('-file',sprintf('%s%s',pathname,filename{f}));
        spknames = [];
        
        for l = 1:size(list,1)
            if findstr('SPK',list{l})
                       spknames = [spknames {list{l}}];
            end
        end
        
        filen = [filename{f}(1:end-4)];
        filen = {filen};
        rat = filen{1}(2:3);
        rwdarms = rwds(filen{1});
                                                                                               
        %collecting the coordinates of the leds in just two variables: 
        %LEDX and LEDY

        clear('test_ledx','test_ledy','test_leds','max_val','max_valx','max_valy','LEDX','LEDY');

        test_ledx=size(LED1_X,1)-size(LED2_X,1);
        test_ledy=size(LED1_Y,1)-size(LED2_Y,1);
        test_leds=size(LED1_X,1)-size(LED1_Y,1);

        if test_ledx==0 && test_ledy==0 && test_leds==0
            max_val=size(LED1_X,1);
        else            
            if test_ledx<0
                max_valx=size(LED2_X,1)+test_ledx;
            elseif test_ledx>0
                max_valx=size(LED1_X,1)-test_ledx;
            else
                max_valx=size(LED1_X,1);
            end

            if test_ledy<0
                max_valy=size(LED2_Y,1)+test_ledy;
            elseif test_ledy>0
                max_valy=size(LED1_Y,1)-test_ledy;
            else
                max_valy=size(LED1_Y,1);
            end

            if max_valx-max_valy>0
                max_val=max_valy;
            else
                max_val=max_valx;
            end
        end

        for b2=1:max_val        

            if LED1_X(b2,1)==0
                LEDX(b2,1)=LED2_X(b2,1);
            elseif LED2_X(b2,1)==0
                LEDX(b2,1)=LED1_X(b2,1);
            else
                LEDX(b2,1)=mean(LED1_X(b2,1),LED2_X(b2,1));
            end

            if LED1_X(b2,2)==0
                LEDX(b2,2)=LED2_X(b2,2);
            elseif LED2_X(b2,2)==0
                LEDX(b2,2)=LED1_X(b2,2);
            else
                LEDX(b2,2)=LED1_X(b2,2);
            end

            if LED1_Y(b2,1)==0
                LEDY(b2,1)=LED2_Y(b2,1);
            elseif LED2_Y(b2,1)==0
                LEDY(b2,1)=LED1_Y(b2,1);
            else
                LEDY(b2,1)=mean(LED1_Y(b2,1),LED2_Y(b2,1));
            end

            if LED1_Y(b2,2)==0
                LEDY(b2,2)=LED2_Y(b2,2);
            elseif LED2_Y(b2,2)==0
                LEDY(b2,2)=LED1_Y(b2,2);
            else
                LEDY(b2,2)=LED1_Y(b2,2);
            end

        end                                                                     

        % Neurons loop
        for s = 1:size(spknames,2)
                       
            % variables storing data for each neuron
            h_spks1 = [];
            h_spks2 = [];
            h_spks3 = [];
            h_spkss = [];
            h_data = [];
            s_spkss1 = [];
            s_data = [];

            counts_spks = []; % for CI neuron
            code = codes([filen{1} ' ' spknames{s}]);
            codess = [codess {code}];
            spk = eval(sprintf('%s',spknames{s}));
            h_s = histcounts(spk,'BinWidth',binsz);
            spk_count = size(spk,2);
            spk_counts = 0;
            count = 0;

            mean_s = mean(h_s);
            std_s = sqrt(mean_s);

            if scron == 1
                spks = scrambled_spks(spk);
            else 
                spks = spk;
            end

            for c = 1:size(cutteddata,2)                 
                if cutteddata(c).tag == 1 
                    if cutteddata(c).efc >= 0.8 & cutteddata(c).efc <= 1.2 | isempty(cutteddata(c).efc)
                                                                                                                                     
                        %calculating velocity and acceleration for each
                        %neuron                                                                        

                        clear('f1','f2','vel','acc','maci','mici','max_acc','min_acc','mean_acc');

                        first_bin=cutteddata(c).tsst(1);
                        v=1;

                        f1=LEDX(:,2)-first_bin;                        
                        f2=f1(f1>0);
                        index=length(LEDX)-length(f2);

                        while first_bin<cutteddata(c).tsend(1)

                            delta_x=LEDX(index+3)-LEDX(index);
                            delta_y=LEDY(index+3)-LEDY(index);
                            square_dist=(delta_x^2)+(delta_y^2);
                            dist=sqrt(square_dist);                                                
                            delta_time=LEDX(index+3,2)-LEDX(index,2);
                            vel(v,1)=dist/delta_time;
                            vel(v,2)=LEDX(index+3,2);

                            if v>1
                                delta_vel=vel(v,1)-vel(v-1,1);
                                delta_time1=vel(v,2)-vel(v-1,2);
                                acc(v-1,1)=delta_vel/delta_time1;
                                acc(v-1,2)=vel(v,2);
                            end

                            v=v+1;
                            index=index+3;                            
                            first_bin=first_bin+delta_time;

                            clear('delta_x','delta_y','square_dist','dist','delta_time','delta_time1')
                        end

                        for c1=2:size(acc,1)-1                            
                            mean_acc(c1-1)=mean(acc(c1-1:c1+1,1));                            
                        end

                        [~,maci]=max(mean_acc(1,:));
                        max_acc=acc(maci,2);
                        [~,mici]=min(mean_acc(1,:));
                        min_acc=acc(mici,2);

                        %alignment loop
                        if alig == 1 % locomotion onset
                             ts = cutteddata(c).tsst(1);                                
                        elseif alig == 2 % peak of speed
                             ts = max_acc;
                        elseif alig == 3 % peak of acceleration                            
                             ts = cutteddata(c).tspk(1);
                        elseif alig == 4 % peak of deceleration
                             ts = min_acc;
                        elseif alig == 5 % locomotion end
                             ts = cutteddata(c).tsend(1);
                        end 
                            
                        ts1=cutteddata(c).tspk(1);    
                                                                                                                                                                                                  
                        if normon == 1
                            
                            c_spks = spks(find(spks > ts-3 & spks < ts+3)); %                             
                            runs_durations = [runs_durations cutteddata(c).tsend(1)-cutteddata(c).tsst(1)];
                            h_spks = histcounts(c_spks,'NumBins', numbins,'BinLimits',[ts-3 ts+3]);
                            
                        else                                

                            if isempty(c_spks)
                                c_spks = -10;
                            end
                            h_spks = histcounts(c_spks,'BinWidth',binsz,'BinLimits',window);
                            
                            if isempty(s_spks)
                                s_spks = -10;
                            end
                            s_spks1 = histcounts(s_spks,'BinWidth',binsz,'BinLimits',window);
                            
                        end
                        
                        % storing all the data
                        h_spkss = [h_spkss; h_spks]; 
                        h_data = [h_data h_spks]; % for std and mean calculations 

                        % storing data separately for each target arm
                        if cutteddata(c).arm == rwdarms(1)
                            h_spks1 = [h_spks1; h_spks];
                        elseif cutteddata(c).arm == rwdarms(2)
                            h_spks2 = [h_spks2; h_spks];
                        elseif cutteddata(c).arm == rwdarms(3)
                            h_spks3 = [h_spks3; h_spks];
                        end 
                        
                        spk_counts = spk_counts + size(c_spks,2);
                        count = count +1;
                        
                    end
                end
            end               

            % calculating mean and std for zscore (based on whole data)
            mean_z = mean(h_data);
            std_z = std(h_data);
            std_z = sqrt(mean_z);
            
            % calculating zscore of the mean data of each arm (and whole data)
            if isempty(h_spks1)                
                meana1 = mean(zeros(5,60));
            else
                meana1 = mean(h_spks1,1);
            end
            meana2 = mean(h_spks2,1);
            meana3 = mean(h_spks3,1);
            means = mean(h_spkss,1);
            
            mean_s = mean(s_spkss1,1);

            z1 = zscore(meana1);
            z2 = zscore(meana2);
            z3 = zscore(meana3);
            zs = zscore(means);
            z1 = (meana1-mean_z)/std_z;
            z2 = (meana2-mean_z)/std_z;
            z3 = (meana3-mean_z)/std_z;
            
            z_s = zscore(mean_s);
            
            if bootson == 1
                meana1 = sum(h_spks1,1);
                meana2 = sum(h_spks2,1);
                meana3 = sum(h_spks3,1);
                means = sum(h_spkss,1);
                mean_z = mean(mean(counts_spks,1));
                std_z = mean(std(counts_spks,1));
                z1 = (meana1-mean_z)/std_z;
                z2 = (meana2-mean_z)/std_z;
                z3 = (meana3-mean_z)/std_z;
                zs = (means-mean_z)/std_z;
            end

            % compiling all the arms and whole data
            if smoothon == 1
                zs1 = [zs1; smooth(z1)'];
                zs2 = [zs2; smooth(z2)'];
                zs3 = [zs3; smooth(z3)'];
                zss = [zss; smooth(zs)'];
                z_ss = [z_ss; smooth(z_s)'];
            else
                zs1 = [zs1; z1];
                zs2 = [zs2; z2];
                zs3 = [zs3; z3];
                zss = [zss; zs];
                z_ss = [z_ss; z_s];
            end

            % Done msgs
            spkss = [spkss {code}];
            fprintf('%s Done!\n',code)                                
             
        end                  
    end
                            
    %sorting data           
    
    for x=1:size(zss,1)
        
        for y=2:size(zss,2)-1
            test_z(y-1)=mean(zss(x,y-1:y+1));
        end
        [~,peak_act(x,1)]=max(test_z);
        
        clear('test_z');
        
    end                

     [~,inda] = sort(peak_act,'descend');

     zsss = zss(inda,:);
     codesss = codess(inda);    
     
%     %writing data to excell file
    
    if writet==1
        
        tabela(:,1)=codesss;

        for p=1:size(zsss,1)
            for p1=1:size(zsss,2)
                tabela(p,p1+1)={zsss(p,p1)};
            end
        end
        
    end
      
    %plotting figure        
    
    if classon==1
        r1=1;    

        for r=1:size(zsss,1)

            pv_test=strcmp(pv,codesss(r));
            pv_test1=find(pv_test);        
            vp_test=strcmp(vp,codesss(r));
            vp_test1=find(vp_test);        
            locoff_test=strcmp(locoff,codesss(r));
            locoff_test1=find(locoff_test);

            if ~isempty(class) %|| ~isempty(
                zsssc(r1,:)=zsss(r,:);
                codesssc(r1)=codesss(r);
                test_r(r1)=r;
                r1=r1+1;            
            end

            clear('pv_test','pv_test1','vp_test','vp_test1','locoff_test','locoff_test1');

        end 
    else
        zsssc=zsss;
        codesssc=codesss;
    end
    
    %clear('zsss','codesss')   
    
    spkss = ones(1,size(zsssc,1));
    cmax = 2;
    cmin = -2;    

    if al==1
        fig = figure();
    end            
    
    subplot(1,5,al);        
    hold on
    imagesc('CData',zsssc,'XData',window(1):binsz:window(2))                
    plot([0 0],[-0.5 size(spkss,2)+1],'k--')
    xlim(window); 
    ylim([.5 size(spkss,2)+0.5 ]);
    colormap(rgb);
    caxis([cmin cmax]);
    xlabel('Time (s)','FontSize',axisfont);
    yticks([1 5 10]);

    if codeon==1
        
        yticks(1:1:length(codesssc));
        ax=gca;
        if classon==0
            ax.FontSize=7.5;
        end
        yticklabels(codesssc);

        for o=1:length(codesssc)

            pv_test=strcmp(pv,codesssc(o));
            pv_test1=find(pv_test);

            vp_test=strcmp(vp,codesssc(o));
            vp_test1=find(vp_test);

            locoff_test=strcmp(locoff,codesssc(o));
            locoff_test1=find(locoff_test);

            if ~isempty(pv_test1)
                ax.YTickLabel{o} = ['\color{green}' ax.YTickLabel{o}];
            elseif ~isempty(vp_test1)
                ax.YTickLabel{o} = ['\color{blue}' ax.YTickLabel{o}];
            elseif ~isempty(locoff_test1)
                ax.YTickLabel{o} = ['\color{black}' ax.YTickLabel{o}];
            else
                ax.YTickLabel{o} = ['\color{gray}' ax.YTickLabel{o}];
            end

            clear('pv_test','pv_test1','vp_test','vp_test1','locoff_test','locoff_test1');
        
        end
        
    end
    
    if al==1
        
        title('Locomotion Onset','FontSize',titlefont);
        ylabel('Neurons','FontSize',axisfont);
        if writet==1
            sheet='Locomotion Onset';
        end
        
    elseif al==2
        
        title('Peak of Acceleration','FontSize',titlefont);
        if writet==1
            sheet='Peak of Acceleration';
        end
        
    elseif al==3        
        
        title('Peak of Speed','FontSize',titlefont); 
        if writet==1
            sheet='Peak of Speed';
        end
        
    elseif al==4
        
        title('Peak of Deceleration','FontSize',titlefont); 
        if writet==1
            sheet='Peak of Deceleration';
        end
        
    elseif al==5
        
        title('Locomotion End','FontSize',titlefont); 
        if writet==1
            sheet='Locomotion End';
        end
        
        cb=colorbar;
        cb.Position = cb.Position + 1e-10;
        set(cb,'Position',[0.93 0.20 0.01 0.5]);
    end
    
    hold off

    if writet==1
        xlswrite(excelfilename,zsss,sheet);
    end

    clear('peak_act','zs','zss','zsss','zsssc','std_z','mean_z','means','codess','codesss','codesssc','h_spkss','h_data','s_spkss1','s_data','peak_act_s');
    codess=[];
    zss = [];
    %set(gcf,'position',[0,0,2000,2000])
                                                                       
end