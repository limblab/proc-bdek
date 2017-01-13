
%
%Set up list of all days for inclusion
FileName = {'Mihili','07112013',  2        ;...
            'Mihili','07152013', [2,3]     ;...
            'Mihili','07192013', [2,3,4]   ;...
            'Mihili','08062013', [2,3]     ;...
            'Mihili','08122013', [2,3]     ;...
            'Mihili','08152013',  2        ;...
            'MrT'   ,'05042013',  2        ;...
            'MrT'   ,'05052013',  2        ;...
            'MrT'   ,'05062013',  2        };
        
        
%  FileName = {'Mihili','07192013',  2        };

[BDAY,GDAY,WDAY,PS] = deal(cell(size(FileName,1),1));
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2});
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5);
    
    alldays(2).tt = [alldays(1).tt ; alldays(2).tt];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tune_ranges = [600 800];
    tune_align = 'target';
    brain_area = 'PMd';
    MODthresh = 0;
    if exist('pdi','var')
        prediction_day_indices = pdi;
    else
        prediction_day_indices = 2;
    end

    [neurons, PD,PD1,nears,fars,nearneurs,farneurs,DFPD,VM_P] = deal(cell(length(tune_ranges)-1,1));
    CO_index = 1;
    for ranges = 1:length(tune_ranges)-1
        tune_range = tune_ranges(ranges:ranges+1);
        % Initialize
        co_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
        %neurons = cell(length(co_units),length(tune_range)-1);

        for i = 1:length(co_units) % Loop through units
             %clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

            % set time bin edges
            t1 = tune_range(1);
            t2 = tune_range(2);

            [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,10),...
                alldays(CO_index).tt,t1,t2,tune_align,i,1);

            neurons{ranges}{i}.tuning = wrapped_tune;
            neurons{ranges}{i}.tuning_low = wrapped_low;
            neurons{ranges}{i}.tuning_high = wrapped_high;

            VM_P{ranges}{i} = p;

            pd_purp = wrapped_cents(abs(circ_dist(wrapped_cents,p(end)))==min(abs(circ_dist(wrapped_cents,p(end)))));
            if p(2)<0
                PD{ranges}{i} = mod(pd_purp+pi,2*pi);
            else
                PD{ranges}{i} = pd_purp;
            end
            PD1{ranges}(i,:) = PD{ranges}{i};

        end

        reppd = repmat(PD1{ranges}',size(alldays(prediction_day_indices).tt,1),1);
        reprch = repmat(alldays(prediction_day_indices).tt(:,10),1,size(reppd,2));

        dfrompd = circ_dist(reppd,reprch);

        DFPD{ranges} = dfrompd;

    end
    centers = round(wrapped_cents*1000)./1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    firing_script;
    clc; fprintf('%d/%d\n',daynum,size(FileName,1)); 
    Cisek_plot_loop;

    BDAY{daynum} = B;
    GDAY{daynum} = G;
    WDAY{daynum} = W;
    
    PS{daynum} = PARAMS;

end
%% Monkey M
[all_bases,all_bases_n,all_mods,all_mods_n,all_width,all_width_n] = deal(cell(5,1));
cfp = {'k','b','r','g'};
figure; hold on; 
subplot(2,3,1); hold on; title('Baseline');
for i = 1:6%length(PS)
    for j = 1:length(PS{i})
        
        all_bases{j}(i,:) = PS{i}{j}(1,:);
        all_bases{j}(i,:) = (PS{i}{j}(1,:) - min(PS{i}{j}(1,:))) ./ ...
            (max(PS{i}{j}(1,:)) - min(PS{i}{j}(1,:)));
        
        plot(all_bases{j}(i,:),cfp{j},'LineWidth',0.3);   
    end
end
for i = 1:3; plot(mean(all_bases{i}(1:6,:)),cfp{i},'LineWidth',3); end

subplot(2,3,2); hold on; title('Modulation');
for i = 1:6%length(PS)
    for j = 1:length(PS{i})
        all_mods{j}(i,:) = PS{i}{j}(2,:);
        all_mods{j}(i,:) = (PS{i}{j}(2,:) - min(PS{i}{j}(2,:))) ./ ...
            (max(PS{i}{j}(2,:)) - min(PS{i}{j}(2,:)));
 
        plot(all_mods{j}(i,:),cfp{j}); 
    end
end
for i = 1:3; plot(mean(all_mods{i}(1:6,:)),cfp{i},'LineWidth',3); end
subplot(2,3,3); hold on; title('Width');
for i = 1:6%length(PS)
    for j = 1:length(PS{i})
        all_width{j}(i,:) = PS{i}{j}(3,:);
        all_width{j}(i,:) = (PS{i}{j}(3,:) - min(PS{i}{j}(3,:))) ./ ...
            (max(PS{i}{j}(3,:)) - min(PS{i}{j}(3,:)));
  
        plot(all_width{j}(i,:),cfp{j});   
    end
end
for i = 1:3; plot(mean(all_width{i}(1:6,:)),cfp{i},'LineWidth',3); end

%% Monkey T
subplot(2,3,4); hold on; title('Baseline');
for i = 7:9%length(PS)
    for j = 1:length(PS{i})
        
        all_bases{j}(i,:) = PS{i}{j}(1,:);
        all_bases{j}(i,:) = (PS{i}{j}(1,:) - min(PS{i}{j}(1,:))) ./ ...
            (max(PS{i}{j}(1,:)) - min(PS{i}{j}(1,:)));
        
        plot(all_bases{j}(i,:),cfp{j},'LineWidth',0.3);   
    end
end
for i = 1:3; plot(mean(all_bases{i}(7:9,:)),cfp{i},'LineWidth',3); end

subplot(2,3,5); hold on; title('Modulation');
for i = 7:9%length(PS)
    for j = 1:length(PS{i})
        all_mods{j}(i,:) = PS{i}{j}(2,:);
        all_mods{j}(i,:) = (PS{i}{j}(2,:) - min(PS{i}{j}(2,:))) ./ ...
            (max(PS{i}{j}(2,:)) - min(PS{i}{j}(2,:)));
 
        plot(all_mods{j}(i,:),cfp{j}); 
    end
end
for i = 1:3; plot(mean(all_mods{i}(7:9,:)),cfp{i},'LineWidth',3); end
subplot(2,3,6); hold on; title('Width');
for i = 7:9%length(PS)
    for j = 1:length(PS{i})
        all_width{j}(i,:) = PS{i}{j}(3,:);
        all_width{j}(i,:) = (PS{i}{j}(3,:) - min(PS{i}{j}(3,:))) ./ ...
            (max(PS{i}{j}(3,:)) - min(PS{i}{j}(3,:)));
  
        plot(all_width{j}(i,:),cfp{j});   
    end
end
for i = 1:3; plot(mean(all_width{i}(7:9,:)),cfp{i},'LineWidth',3); end


%% 
Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
MrT_inds = find(strcmp(FileName(:,1),'MrT'));

figure; hold on; 
subplot(1,3,1); hold on; title('Baseline','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot([daynum + x_off, daynum+x_off+0.25], BDAY{daynum}{i}(1,:),'.-','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off,      [BDAY{daynum}{i}(2,1) BDAY{daynum}{i}(3,1)],'Color',cfp{i},'LineWidth',1);
        plot(daynum*[1 1]+x_off+0.25, [BDAY{daynum}{i}(2,2) BDAY{daynum}{i}(3,2)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);

subplot(1,3,2); hold on; title('Modulation','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot([daynum + x_off, daynum+x_off+0.25], GDAY{daynum}{i}(1,:),'.-','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off,      [GDAY{daynum}{i}(2,1) GDAY{daynum}{i}(3,1)],'Color',cfp{i},'LineWidth',1);
        plot(daynum*[1 1]+x_off+0.25, [GDAY{daynum}{i}(2,2) GDAY{daynum}{i}(3,2)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);
 
subplot(1,3,3); hold on; title('Width','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot([daynum + x_off, daynum+x_off+0.25], WDAY{daynum}{i}(1,:),'.-','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off,      [WDAY{daynum}{i}(2,1) WDAY{daynum}{i}(3,1)],'Color',cfp{i},'LineWidth',1);
        plot(daynum*[1 1]+x_off+0.25, [WDAY{daynum}{i}(2,2) WDAY{daynum}{i}(3,2)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);

%% 
Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
MrT_inds = find(strcmp(FileName(:,1),'MrT'));

figure; hold on; 
subplot(1,3,1); hold on; title('Baseline','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot(BDAY{daynum}{i}(1,1), BDAY{daynum}{i}(1,2),'.','Color',cfp{i},'MarkerSize',10);
        plot(BDAY{daynum}{i}(2:3,1),[1 1]*BDAY{daynum}{i}(1,2),'Color',cfp{i},'LineWidth',1);
        plot([1 1]*BDAY{daynum}{i}(1,1),BDAY{daynum}{i}(2:3,2),'Color',cfp{i},'LineWidth',1);
    end
    plot(cellfun(@(x) x(1,1),BDAY{daynum}),cellfun(@(x) x(1,2),BDAY{daynum}),'k');
end

%%

subplot(1,3,2); hold on; title('Modulation','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot([daynum + x_off, daynum+x_off+0.25], GDAY{daynum}{i}(1,:),'.-','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off,      [GDAY{daynum}{i}(2,1) GDAY{daynum}{i}(3,1)],'Color',cfp{i},'LineWidth',1);
        plot(daynum*[1 1]+x_off+0.25, [GDAY{daynum}{i}(2,2) GDAY{daynum}{i}(3,2)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);
 
subplot(1,3,3); hold on; title('Width','FontSize',16);
for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if length(BDAY{daynum}) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','r'}; end
    for i = 1:length(BDAY{daynum})
        plot([daynum + x_off, daynum+x_off+0.25], WDAY{daynum}{i}(1,:),'.-','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off,      [WDAY{daynum}{i}(2,1) WDAY{daynum}{i}(3,1)],'Color',cfp{i},'LineWidth',1);
        plot(daynum*[1 1]+x_off+0.25, [WDAY{daynum}{i}(2,2) WDAY{daynum}{i}(3,2)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);
    