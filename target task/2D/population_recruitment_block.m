%%
session_nums = 2;
tbin = 2;
symmetric_plot = 0;
boot_count = 1000;
align_index = 10;
cutoffFR = 0;
binwindow = 1000;
smooth_type = 'bin';
num_bins = 6;
num_blocks = 5;

blocks = round(linspace(1,size(alldays(session_nums).tt,1),num_blocks+1));

plotcolors = {'r','g','b'};
figure;
for bs = 1:length(blocks)-1
    
    trial_block = blocks(bs):blocks(bs+1);
    tt_pred = alldays(session_nums).tt(trial_block,:); 
    
    fprintf('time bin: %d/%d\n',bs,length(blocks)-1);
    subplot(1,length(blocks)-1,bs); hold on;

    cur_tt = tt_pred; % trial table
    dub_tt = [cur_tt; cur_tt]; % duplicated trial table

    ecounts = expect_counts{1}{tbin};
    acounts = trial_counts{1}{tbin};

    minfr = zeros(length(neurons{tbin}),1);
    maxfr = zeros(length(neurons{tbin}),1);

    %%% Get preferred directions from tuning descriptions %%%
    pds = cell(length(neurons{tbin}),1);
    for i = 1:length(neurons{tbin})

        minfr(i) = min(neurons{tbin}{i}.tuning);
        maxfr(i) = max(neurons{tbin}{i}.tuning);

       [heights,pds{i}] = findpeaks(neurons{tbin}{i}.tuning);

    end

    %%% Only use single-peaked neurons %%%
    numpeaks = cellfun(@length,pds);
    well_tuned = find(numpeaks==1); % single peaked tuning curve
    %well_tuned = find(numpeaks>0);
    prefd = round(cellfun(@mean,pds));

    %%% Extract movement directions
    move_dirs = cur_tt(:,align_index);

    %%% Of single-peaked neurons, eliminate those with too low FR %%%
    maxtc = maxfr(well_tuned);
    mintc = minfr(well_tuned);
    badn = find(maxfr - minfr <= cutoffFR);
    vwell_tuned = well_tuned(~ismember(well_tuned,badn));

    %%% Compile distances from PD --> Movement direction
    perc_max = zeros(length(vwell_tuned),length(move_dirs));
    from_move = zeros(length(vwell_tuned),length(move_dirs));
    eperc_max = zeros(length(vwell_tuned),length(move_dirs));
    for i = 1:length(move_dirs)

        % distances
        from_move(:,i) = wrapped_cents(prefd(vwell_tuned)) - move_dirs(i);

        % constrain to -pi:pi
        from_move(from_move > pi) = from_move(from_move > pi) - 2*pi;

        % normalize all neurons to calculate activation in % modulation
        perc_max(:,i) = (acounts(trial_block(i),vwell_tuned)' - minfr(vwell_tuned))./...
            (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Actual Spikes
        eperc_max(:,i) = (ecounts(trial_block(i),vwell_tuned)' - minfr(vwell_tuned))./...
            (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Predicted (CO) spikes

    end

    %%% Get likelihood trials %%%
    likes = unique(cur_tt(:,3));
    likind = cell(length(likes),1);
    if symmetric_plot == 1
        for z = 1:length(likes)
            likind{z} = find(dub_tt(:,3)==likes(z));
        end

        from_move = [from_move, -from_move]; %#ok<AGROW>
        perc_max = [perc_max, perc_max]; %#ok<AGROW>
        eperc_max = [eperc_max, eperc_max]; %#ok<AGROW>

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
    allfm_all = reshape(from_move,numel(from_move),1);
    allpm_all = reshape(perc_max,numel(perc_max),1);
    allepm_all = reshape(eperc_max,numel(eperc_max),1);

    %%% Create mask for convolution window averaging
    mask = ones(1,binwindow)/binwindow;

    %%% For all likelihoods, compute smoothed and binned metrics
    for z = 1:length(likes)

        fprintf('Likelihood %d/%d\n',z,length(likes));

        allfm{z} = reshape(from_move(:,likind{z}),numel(from_move(:,likind{z})),1);
        allpm{z} = reshape(perc_max(:,likind{z}),numel(perc_max(:,likind{z})),1);
        allepm{z} = reshape(eperc_max(:,likind{z}),numel(eperc_max(:,likind{z})),1);

        [fm_increasing{z},orderfm{z}] = sortrows(allfm{z});

        pm_increasing{z} = allpm{z}(orderfm{z});
        epm_increasing{z} = allepm{z}(orderfm{z});

        switch smooth_type

            case 'bin'

            binedges = -pi:(pi/num_bins):pi;
            for pbin = 1:length(binedges)-2

                binbounds = [binedges(pbin) binedges(pbin+2)];
                bininds = fm_increasing{z} >= binbounds(1) & fm_increasing{z} <= binbounds(2);

                binndfminc{z}(1,pbin) = nanmean(fm_increasing{z}(bininds));
                binndpminc{z}(1,pbin) = nanmean(pm_increasing{z}(bininds));
                binndepminc{z}(1,pbin) = nanmean(epm_increasing{z}(bininds));

                [binndpminc{z}(2,pbin), binndpminc{z}(3,pbin)] = ...
                    boot_bounds(boot_count,@nanmean,pm_increasing{z}(bininds),2.5,97.5);
            end

            case 'conv'

                averaged_fm{z} = conv(fm_increasing{z},mask,'valid');

                averaged_pm{z} = conv(pm_increasing{z},mask,'valid');
                boot_wrap = @(x) conv(pm_increasing{z},hist(x,1:binwindow)/binwindow,'valid');
                [boot_pm_low{z}, boot_pm_high{z}] = boot_bounds(boot_count,boot_wrap,1:binwindow,2.5,97.5);

                averaged_epm{z} = conv(epm_increasing{z},mask,'valid');

        end
        legent{z} = num2str(likes(z));

    end

    for q = 1:length(session_nums)

        %%% PLOTTING %%%
        switch smooth_type
            case 'conv'

            %subplot(1,length(session_nums),q); hold on;
            for z = 1:length(likes)
                plot(averaged_fm{z},averaged_pm{z},'Color',plotcolors{z},'LineWidth',3);    
            end

            h = legend(legent);
            v = get(h,'title');
            set(v,'string','Likelihoods');

            for z = 1:length(likes)    
                patch([averaged_fm{z}' fliplr(averaged_fm{z}')],...
                      [boot_pm_low{z}' fliplr(boot_pm_high{z}')],...
                      plotcolors{z},'FaceAlpha',.5,'EdgeAlpha',0);

                plot(averaged_fm{z},averaged_epm{z},'k--','LineWidth',2);
            end 
            ylim([0 1.5]);
            title(sprintf('Prior: %d',priors{session_nums}.val),'FontSize',18);

            case 'bin'

            %subplot(1,length(session_nums),q); hold on;
            for z = 1:length(likes)
                plot(binndfminc{z},binndpminc{z}(1,:),'Color',plotcolors{z},'LineWidth',3);    
            end

            h = legend(legent);
            v = get(h,'title');
            set(v,'string','Likelihoods');

            for z = 1:length(likes)    
                patch([binndfminc{z} fliplr(binndfminc{z})],...
                      [binndpminc{z}(2,:) fliplr(binndpminc{z}(3,:))],...
                      plotcolors{z},'FaceAlpha',.5,'EdgeAlpha',0);

                plot(binndfminc{z},binndepminc{z},'k--','LineWidth',2);
            end 
            ylim([0 1.5]);
            title(sprintf('Prior: %d',priors{session_nums}.val),'FontSize',18);    
        end
    end
end

