%%
session_nums = [2];
tbin = 3;
symmetric_plot = 1;
boot_count = 100;
align_index = 10;
cutoffFR = 2;
binwindow = 1000;
smooth_type = 'bin';
num_bins = 5;

plotcolors = {'r','g','b'};

%%
%%%%%%%%%%%%%%%%%%%%%%%
[FM,PM,EPM] = deal(cell(length(session_nums)));

for neurnum = 1:191

    clc; fprintf('Neuron: %d/%d\n',neurnum,191);
    tuning_array = neurboot{tbin}{neurnum}.tuning;

    modulations = (max(tuning_array,[],2) - min(tuning_array,[],2))';
    good_neurs = find(modulations > 0);

    ecountsneur = cell(1,length(session_nums));
    for z = 1:length(session_nums)

        day_pred = session_nums(z);
        tt_pred = alldays(day_pred).tt;

        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        for i = 1:size(alldays(session_nums(z)).tt,1)    
            for q = 1:size(tuning_array,1)       
                ecountsneur{z}(i,q) = tuning_array(q,round(100*alldays(session_nums(z)).tt(i,align_index)));
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%

    for q = 1:length(session_nums) % loop through prior sessions

        %fprintf('Prior %d/%d\n',q,length(session_nums));

        % current prior session
        prisess = session_nums(q);

        cur_tt = alldays(prisess).tt; % trial table
        dub_tt = [cur_tt; cur_tt]; % duplicated trial table

        ecounts = ecountsneur{prisess-1};
        acounts = repmat(trial_counts{prisess-1}{tbin}(:,neurnum),1,size(tuning_array,1));

        minfr = min(neurboot{tbin}{neurnum}.tuning,[],2);
        maxfr = max(neurboot{tbin}{neurnum}.tuning,[],2);

        %%% Get preferred directions from tuning descriptions %%%
        pds = cell(size(neurboot{tbin}{neurnum}.tuning,1),1);
        for i = 1:length(minfr)
           [heights,temp_pd] = findpeaks(neurboot{tbin}{neurnum}.tuning(i,:));
           peak_pd = temp_pd(heights==max(heights));
           
           pds{i} = peak_pd;
        end

        %%% Only use single-peaked neurons %%%
        numpeaks = cellfun(@length,pds);
        well_tuned = find(numpeaks==1); % single peaked tuning curve
        %well_tuned = find(numpeaks>0);
        prefd = round(cellfun(@mean,pds));

        %%% Extract movement directions
        move_dirs = alldays(prisess).tt(:,align_index);

        %%% Of single-peaked neurons, eliminate those with too low FR %%%
        maxtc = maxfr(well_tuned);
        mintc = minfr(well_tuned);
        badn = find(maxfr - minfr <= cutoffFR);
        vwell_tuned = well_tuned(~ismember(well_tuned,badn));

        %%% Compile distances from PD --> Movement direction
        perc_max = zeros(length(vwell_tuned),length(move_dirs));
        from_move = zeros(length(vwell_tuned),length(move_dirs));
        eperc_max = zeros(length(vwell_tuned),length(move_dirs));
        
        if ~isempty(vwell_tuned)
        
            for i = 1:length(move_dirs)

                % distances
                from_move(:,i) = wrapped_cents(prefd(vwell_tuned)) - move_dirs(i);

                % constrain to -pi:pi
                from_move(from_move > pi) = from_move(from_move > pi) - 2*pi;

                % normalize all neurons to calculate activation in % modulation
                perc_max(:,i) = (acounts(i,vwell_tuned)' - minfr(vwell_tuned))./...
                    (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Actual Spikes
                eperc_max(:,i) = (ecounts(i,vwell_tuned)' - minfr(vwell_tuned))./...
                    (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Predicted (CO) spikes

            end

            FM{q}{neurnum} = from_move;
            PM{q}{neurnum} = perc_max;
            EPM{q}{neurnum} = eperc_max;
            
        else
            
        FM{q}{neurnum} = [];
        PM{q}{neurnum} = [];
        EPM{q}{neurnum} = [];
        
        end
     
    end
end

%%
[ALL_FM,ALL_PM,ALL_EPM] = deal(cell(length(session_nums),1));
for q = 1:length(session_nums) 
    ALL_FM{q} = vertcat(FM{q}{:}); fprintf('Concatenating (1/3)...');
    ALL_PM{q} = vertcat(PM{q}{:}); clc; fprintf('Concatenating (2/3)...');
    ALL_EPM{q} = vertcat(EPM{q}{:}); clc; fprintf('Concatenating (3/3)...');

    prisess = session_nums(q);
    cur_tt = alldays(prisess).tt; % trial table
    dub_tt = [cur_tt; cur_tt]; % duplicated trial table
%%% Get likelihood trials %%%
   
    likes = unique(cur_tt(:,3));
    likind = cell(length(likes),1);
    if symmetric_plot == 1
        for z = 1:length(likes)
            likind{z} = find(dub_tt(:,3)==likes(z));
        end

        ALL_FM{q} = [ALL_FM{q}, -ALL_FM{q}];
        ALL_PM{q} = [ALL_PM{q}, ALL_PM{q}];
        ALL_EPM{q} = [ALL_EPM{q}, ALL_EPM{q}]; 

    else
        for z = 1:length(likes)
            likind{z} = find(cur_tt(:,3)==likes(z));
        end
    end

    %%% Initialize %%%

    [allfm,allpm,allepm,fm_increasing,orderfm,pm_increasing,epm_increasing,...
        binndfminc,binndpminc,binndepminc,averaged_fm,averaged_pm,...
        boot_pm_low,boot_pm_high,averaged_epm,legent] = deal(cell(length(likes),1));
    
    %%% Reshape to combine all trials and all neurons
    allfm_all = reshape(ALL_FM{q},numel(ALL_FM{q}),1); clc; fprintf('Reshaping (1/3)...');
    allpm_all = reshape(ALL_PM{q},numel(ALL_PM{q}),1); clc; fprintf('Reshaping (2/3)...');
    allepm_all = reshape(ALL_EPM{q},numel(ALL_EPM{q}),1); clc; fprintf('Reshaping (3/3)...');

    %%% Create mask for convolution window averaging
    mask = ones(1,binwindow)/binwindow;
    
    %%% For all likelihoods, compute smoothed and binned metrics
    for z = 1:length(likes)
        
        fprintf('Likelihood %d/%d\n',z,length(likes));
        
        allfm{z} = reshape(ALL_FM{q}(:,likind{z}),numel(ALL_FM{q}(:,likind{z})),1);
        allpm{z} = reshape(ALL_PM{q}(:,likind{z}),numel(ALL_PM{q}(:,likind{z})),1);
        allepm{z} = reshape(ALL_EPM{q}(:,likind{z}),numel(ALL_EPM{q}(:,likind{z})),1);
    
        [fm_increasing{z},orderfm{z}] = sortrows(allfm{z});
        
        pm_increasing{z} = allpm{z}(orderfm{z});
        epm_increasing{z} = allepm{z}(orderfm{z});
        
        switch smooth_type
            
            case 'bin'
        
                binedges = -pi:(pi/num_bins):pi;
                for pbin = 1:length(binedges)-2
                    
                    fprintf('TC bin %d/%d\n',pbin,length(binedges)-2);

                    binbounds = [binedges(pbin) binedges(pbin+2)];
                    bininds = fm_increasing{z} >= binbounds(1) & fm_increasing{z} <= binbounds(2);

                    binndfminc{z}{q}(1,pbin) = nanmean(fm_increasing{z}(bininds));
                    binndpminc{z}{q}(1,pbin) = nanmean(pm_increasing{z}(bininds));
                    binndepminc{z}{q}(1,pbin) = nanmean(epm_increasing{z}(bininds));

                    [binndpminc{z}{q}(2,pbin), binndpminc{z}{q}(3,pbin)] = ...
                        boot_bounds(boot_count,@nanmean,pm_increasing{z}(bininds),2.5,97.5);
                end
            
            case 'conv'

                averaged_fm{z}{q} = conv(fm_increasing{z},mask,'valid');

                averaged_pm{z}{q} = conv(pm_increasing{z},mask,'valid');
                boot_wrap = @(x) conv(pm_increasing{z},hist(x,1:binwindow)/binwindow,'valid');
                [boot_pm_low{z}{q}, boot_pm_high{z}{q}] = boot_bounds(boot_count,boot_wrap,1:binwindow,2.5,97.5);

                averaged_epm{z}{q} = conv(epm_increasing{z},mask,'valid');
                
        end
        legent{z} = num2str(likes(z));
        
    end
end
%%
figure; hold on;
for q = 1:length(session_nums)
    %%% PLOTTING %%%
    switch smooth_type
        case 'conv'
   
        subplot(1,length(session_nums),q); hold on;
        for z = 1:length(likes)
            plot(averaged_fm{z}{q},averaged_pm{z}{q},'Color',plotcolors{z},'LineWidth',3);    
        end

        h = legend(legent);
        v = get(h,'title');
        set(v,'string','Likelihoods');

        for z = 1:length(likes)    
            patch([averaged_fm{z}{q}' fliplr(averaged_fm{z}{q}')],...
                  [boot_pm_low{z}{q}' fliplr(boot_pm_high{z}{q}')],...
                  plotcolors{z},'FaceAlpha',.5,'EdgeAlpha',0);

            plot(averaged_fm{z}{q},averaged_epm{z}{q},'k--','LineWidth',2);
        end 
        ylim([0 1.5]);
        title(sprintf('Prior: %d',priors{session_nums(q)}.val),'FontSize',18);
        
        case 'bin'
            
        subplot(1,length(session_nums),q); hold on;
        for z = 1:length(likes)
            plot(binndfminc{z}{q},binndpminc{z}{q}(1,:),'Color',plotcolors{z},'LineWidth',3);    
        end

        h = legend(legent);
        v = get(h,'title');
        set(v,'string','Likelihoods');

        for z = 1:length(likes)    
            patch([binndfminc{z}{q} fliplr(binndfminc{z}{q})],...
                  [binndpminc{z}{q}(2,:) fliplr(binndpminc{z}{q}(3,:))],...
                  plotcolors{z},'FaceAlpha',.5,'EdgeAlpha',0);

            plot(binndfminc{z}{q},binndepminc{z}{q},'k--','LineWidth',2);
        end 
        ylim([0 1.5]);
        title(sprintf('Prior: %d',priors{session_nums(q)}.val),'FontSize',18);    
    end

end

