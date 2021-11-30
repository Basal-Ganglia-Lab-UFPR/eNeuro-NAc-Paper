% Behavioral analysis
%
% Ptg Paper - gabriell.baltazar99@gmail.com
% Da Cunha's Lab - Dpt of Pharmacology - Universidade Federal do Paraná
% Curitiba, Brazil
%
% Script to compute mean latencies of rats runs to each arm (L, M or H) and
% choice preference over the same arms
%

%selecting locomotion data to work with
[filename, pathname] = uigetfile('*.mat', 'Select mat data files', 'Multiselect', 'on');
if ~iscell(filename)
    filename={filename};
end

y=1;
w=1;
choices(1:66,1:12)=0;
rat2='r04';
index_l=1;
index_m=1;
index_h=1;

for a=1:size(filename,2)  %carrega arquivos  
    
    clear('data','fields','fields1','fields2');
    
    data=load(sprintf('%s%s', pathname, filename{a}));     %load the files selected in the beggining of the script one at a time    
    
    fields=fieldnames(data);
    fields1=strncmp(fields,'SPK',3);
    fields2=find(fields1);        
        
    clear('e5','e6','animal','neuron','name','name1','rat','rewards');

    e5=strcat('data.',(fields{fields2(1)}));                
    e6=eval(e5);                

    animal=filename{a}(1:end-4);
    neuron=fields{fields2(1)};
    name=[animal ' ' neuron];
    name1={codes(sprintf('%s',name))};
    neurons(a,1)=name1;

    rat=filename{a}(1:3);
    rewards=rwds(rat); 

    l=1;
    m=1;
    h=1; 
    z=1;
    counter=0;                        

    for c=1:length(data.cutteddata)                      
        if data.cutteddata(c).tag(1)==1                                

            if data.cutteddata(c).arm(1)==rewards(1)
                rew=1;
            elseif data.cutteddata(c).arm(1)==rewards(2)
                rew=2;
            elseif data.cutteddata(c).arm(1)==rewards(3)
                rew=3;
            else
                rew=9;
            end

            if rew==1
                latencies_l(l)=data.cutteddata(c).tsend(1)-data.cutteddata(c).tsst(1);
                l=l+1;
            elseif rew==2
                latencies_m(m)=data.cutteddata(c).tsend(1)-data.cutteddata(c).tsst(1);
                m=m+1;
            elseif rew==3
                latencies_h(h)=data.cutteddata(c).tsend(1)-data.cutteddata(c).tsst(1);
                h=h+1;
            end

            if c~=1 && data.cutteddata(c).trial(1)~=data.cutteddata(c-1).trial(1)
                d=1;
                e=1;
                counter=0;
                
                if data.cutteddata(c).arm(1)==rewards(1)
                    arm_preference(z)='L';
                elseif data.cutteddata(c).arm(1)==rewards(2)
                    arm_preference(z)='M';
                elseif data.cutteddata(c).arm(1)==rewards(3)
                    arm_preference(z)='H';
                else
                    arm_preference(z)=NaN;
                end
                
                z=z+1;                
            end

            if c==1
                choice=0;
                counter=0;
                d=1;
                e=1;
                
                if data.cutteddata(c).arm(1)==rewards(1)
                    arm_preference(z)='L';
                elseif data.cutteddata(c).arm(1)==rewards(2)
                    arm_preference(z)='M';
                elseif data.cutteddata(c).arm(1)==rewards(3)
                    arm_preference(z)='H';
                else
                    arm_preference(z)=NaN;
                end
                
                z=z+1;                
            else                    
                if d==1
                    if counter==0
                        choice=0;                            
                    elseif counter==1
                        choice=3;
                    elseif counter==2
                        choice=6;
                        d=0;
                    end

                end

            end                                               

            if e==1
                choices(y,choice+rew)=choices(y,choice+rew)+1;
                counter=counter+1;
                if d==0
                    e=0;
                end
            end
        end                        
    end                

    rat1=rat;

    if strcmp(rat1,rat2)==0                        
        
        w=1;
        
%        xlswrite('behavioral_analysis(2).xlsx',latencies_l_rat,sprintf('%s_l',rat2));
%        xlswrite('behavioral_analysis(2).xlsx',latencies_m_rat,sprintf('%s_m',rat2));
%        xlswrite('behavioral_analysis(2).xlsx',latencies_h_rat,sprintf('%s_h',rat2));
        xlswrite('behavioral_analysis(3).xlsx',arm_preference_rat,rat2);

        clear('latencies_l_rat','latencies_m_rat','latencies_h_rat','arm_preference_rat');

    end

    if length(latencies_l)>12
        latencies_l_1=latencies_l(1:12);
    elseif length(latencies_l)<12
        latencies_l_1=latencies_l;
        latencies_l_1(length(latencies_l)+1:12)=NaN;
    else
        latencies_l_1=latencies_l;
    end
    
    if length(latencies_m)>12
        latencies_m_1=latencies_m(1:12);
    elseif length(latencies_m)<12
        latencies_m_1=latencies_m;
        latencies_m_1(length(latencies_m)+1:12)=NaN;
    else
        latencies_m_1=latencies_m;
    end
    
    if length(latencies_h)>12
        latencies_h_1=latencies_h(1:12);
    elseif length(latencies_h)<12
        latencies_h_1=latencies_h;
        latencies_h_1(length(latencies_h)+1:12)=NaN;
    else
        latencies_h_1=latencies_h;
    end    
    
    latencies_l_rat(w,1:12)=latencies_l_1;
    latencies_m_rat(w,1:12)=latencies_m_1;
    latencies_h_rat(w,1:12)=latencies_h_1;
    
    if length(arm_preference)>12
        arm_preference_1=arm_preference(1:12);
    elseif length(arm_preference)<12
        arm_preference_1=arm_preference;
        arm_preference_1(length(arm_preference)+1:12)=NaN;
    else
        arm_preference_1=arm_preference;
    end
    
    arm_preference_rat(w,1:12)=arm_preference_1;
    
    w=w+1;

    y=y+1;

    rat2=rat;

    clear('latencies_l','latencies_m','latencies_h','arm_preference','arm_preference_1');    
end