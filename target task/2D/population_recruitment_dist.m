%%
session_nums = 2;
tbins = [1];
symmetric_plot = 0;
boot_count = 1000;
align_index = 10;
cutoffFR = 0;
binwindow = 1000;
smooth_type = 'bin_nonoverlap';
num_bins = 5;

plotcolors = {'r','g','b','c'};
figure;
[vwt,tbpds,FM,PM] = deal(cell(length(tbins),1));
for time_bins = 1:length(tbins)
    fprintf('time bin: %d/%d\n',time_bins,length(tbins));
    tbin = tbins(time_bins);
    subplot(1,length(tbins),time_bins); hold on
    for q = 1:length(session_nums) % loop through prior sessions

        fprintf('Prior %d/%d\n',q,length(session_nums));

        % current prior session
        prisess = session_nums(q);

        cur_tt = alldays(prisess).tt; % trial table
        dub_tt = [cur_tt; cur_tt]; % duplicated trial table

        ecounts = expect_counts{prisess-1}{tbin};
        acounts = trial_counts{prisess-1}{tbin};

        minfr = zeros(length(neurons{tbin}),1);
        maxfr = zeros(length(neurons{tbin}),1);

        %%% Get preferred directions from tuning descriptions %%%
        [heights,pds] = deal(cell(length(neurons{tbin}),1));
        for i = 1:length(neurons{tbin})

            minfr(i) = min(neurons{tbin}{i}.tuning);
            maxfr(i) = max(neurons{tbin}{i}.tuning);

           [heights{i},pds{i}] = findpeaks(neurons{tbin}{i}.tuning);

        end

        %%% Only use single-peaked neurons %%%
        numpeaks = cellfun(@length,pds);
        well_tuned = find(numpeaks>=1); % single peaked tuning curve
        prefd = round(cellfun(@mean,pds));
%         well_tuned = find(numpeaks>0);
%         anonfunc = @(x) find(x==max(x),1,'first');
%         pdinds = cellfun(anonfunc,heights,'UniformOutput',0);
%         prefd = NaN(length(pdinds),1);
%         for i = 1:length(pdinds)
%             if ~isempty(pdinds{i} > 0)
%                 prefd(i,:) = pds{i}(pdinds{i});
%             end
%         end

        %%% Extract movement directions
        move_dirs = cur_tt(:,align_index);
%         move_dirs = zeros(size(cur_tt,1),1);

        %%% Of single-peaked neurons, eliminate those with too low FR %%% 
        maxtc = maxfr(well_tuned);
        mintc = minfr(well_tuned);
        badn = find(maxfr - minfr <= cutoffFR);
        vwell_tuned = well_tuned(~ismember(well_tuned,badn));
        
        vwt{time_bins} = vwell_tuned;
        tbpds{time_bins} = prefd(vwell_tuned);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %vwell_tuned = vwell_tuned(ismember(vwell_tuned,top100));
        %vwell_tuned = vwell_tuned(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        %%% Compile distances from PD --> Movement direction
        perc_max = zeros(length(vwell_tuned),length(move_dirs));
        from_move = zeros(length(vwell_tuned),length(move_dirs));
        eperc_max = zeros(length(vwell_tuned),length(move_dirs));
        perc_exp = zeros(length(vwell_tuned),length(move_dirs));
        for i = 1:length(move_dirs)

            % distances
            %from_move(:,i) = circ_dist(wrapped_cents(prefd(vwell_tuned)),move_dirs(i));

            % constrain to -pi:pi
            %from_move(from_move > pi) = from_move(from_move > pi) - 2*pi;

            % normalize all neurons to calculate activation in % modulation
            perc_max(:,i) = (acounts(i,vwell_tuned)' - minfr(vwell_tuned))./...
                (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Actual Spikes
            eperc_max(:,i) = (ecounts(i,vwell_tuned)' - minfr(vwell_tuned))./...
                (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Predicted (CO) spikes
            perc_exp(:,i) = (acounts(i,vwell_tuned)' - ecounts(i,vwell_tuned)')./...
                (maxfr(vwell_tuned)-minfr(vwell_tuned)); % Actual Spikes
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            perc_max = perc_exp;
            from_move = eperc_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%                 binedges = -pi:(pi/num_bins):pi;
                binedges = 0:(1/num_bins):1;
                bininds = cell(1,length(binedges)-2);
                nbininds = cell(1,length(binedges)-2);
                for pbin = 1:length(binedges)-2

                    binbounds = [binedges(pbin) binedges(pbin+2)];
                    bininds{pbin} = fm_increasing{z} >= binbounds(1) & fm_increasing{z} <= binbounds(2);
%                 end
%                 minbinnum = min(cellfun(@sum,bininds));
% 
%                 for pbin = 1:length(binedges)-2
%                     
%                     actinds = find(bininds{pbin}==1);
%                     shuffchoices = randperm(length(actinds));
%                     nbininds{pbin} = actinds(shuffchoices(1:minbinnum));

                    binndfminc{z}{q}(1,pbin) = nanmean(fm_increasing{z}(bininds{pbin}));
                    binndpminc{z}{q}(1,pbin) = nanmean(pm_increasing{z}(bininds{pbin}));
                    binndepminc{z}{q}(1,pbin) = nanmean(epm_increasing{z}(bininds{pbin}));

                    [binndpminc{z}{q}(2,pbin), binndpminc{z}{q}(3,pbin)] = ...
                        boot_bounds(boot_count,@nanmean,pm_increasing{z}(bininds{pbin}),2.5,97.5);
                end
                
                case 'conv'

                    averaged_fm{z}{q} = conv(fm_increasing{z},mask,'valid');

                    averaged_pm{z}{q} = conv(pm_increasing{z},mask,'valid');
                    boot_wrap = @(x) conv(pm_increasing{z},hist(x,1:binwindow)/binwindow,'valid');
                    [boot_pm_low{z}{q}, boot_pm_high{z}{q}] = boot_bounds(boot_count,boot_wrap,1:binwindow,2.5,97.5);

                    averaged_epm{z}{q} = conv(epm_increasing{z},mask,'valid');
                    
                case 'bin_equal'

                binedges = 1:num_bins:length(fm_increasing{z});
                bininds = cell(1,length(binedges)-1);
                nbininds = cell(1,length(binedges)-1);
                for pbin = 1:length(binedges)-1

                    binbounds = [binedges(pbin) binedges(pbin+1)];

                    binndfminc{z}{q}(1,pbin) = nanmean(fm_increasing{z}(binbounds(1):binbounds(2)));
                    binndpminc{z}{q}(1,pbin) = nanmean(pm_increasing{z}(binbounds(1):binbounds(2)));
                    binndepminc{z}{q}(1,pbin) = nanmean(epm_increasing{z}(binbounds(1):binbounds(2)));

                    [binndpminc{z}{q}(2,pbin), binndpminc{z}{q}(3,pbin)] = ...
                        boot_bounds(boot_count,@nanmean,pm_increasing{z}(binbounds(1):binbounds(2)),2.5,97.5);
                end
                
                case 'bin_nonoverlap'

%                 binedges = -pi:(pi/num_bins):pi;
                binedges = 0:(1/num_bins):1;
                bininds = cell(1,length(binedges)-1);
                nbininds = cell(1,length(binedges)-1);
                for pbin = 1:length(binedges)-1

                    binbounds = [binedges(pbin) binedges(pbin+1)];
                    bininds{pbin} = fm_increasing{z} >= binbounds(1) & fm_increasing{z} <= binbounds(2);


                    binndfminc{z}{q}(1,pbin) = nanmean(fm_increasing{z}(bininds{pbin}));
                    binndpminc{z}{q}(1,pbin) = nanmean(pm_increasing{z}(bininds{pbin}));
                    binndepminc{z}{q}(1,pbin) = nanmean(epm_increasing{z}(bininds{pbin}));

                    [binndpminc{z}{q}(2,pbin), binndpminc{z}{q}(3,pbin)] = ...
                        boot_bounds(boot_count,@nanmean,pm_increasing{z}(bininds{pbin}),2.5,97.5);

                end

            end
            legent{z} = num2str(likes(z));

        end
    end

    for q = 1:length(session_nums)

        %%% PLOTTING %%%
        switch smooth_type
            case 'conv'

            %subplot(1,length(session_nums),q); hold on;
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

            %subplot(1,length(session_nums),q); hold on;
            for z = 1:length(likes)
                plot(binndfminc{z}{q},binndpminc{z}{q}(1,:),'Color',plotcolors{z},'LineWidth',3);    
                %plot(fm_increasing{z},pm_increasing{z},'.','Color',plotcolors{z});
            end

            h = legend(legent);
            v = get(h,'title');
            set(v,'string','Likelihoods');

            for z = 1:length(likes)    
                patch([binndfminc{z}{q} fliplr(binndfminc{z}{q})],...
                      [binndpminc{z}{q}(2,:) fliplr(binndpminc{z}{q}(3,:))],...
                      plotcolors{z},'FaceAlpha',.25,'EdgeAlpha',0);

                plot(binndfminc{z}{q},zeros(size(binndfminc{z}{q})),'--','Color','k','LineWidth',2);
            end 
            
            catall = horzcat(binndpminc{:});
            catall2 = horzcat(catall{:});
            %ylim([0 ceil(max(catall2(:)))]);
            title(sprintf('Prior: %d',priors{session_nums(q)}.val),'FontSize',18);    
            
            case 'bin_equal'

            %subplot(1,length(session_nums),q); hold on;
            for z = 1:length(likes)
                plot(binndfminc{z}{q},binndpminc{z}{q}(1,:),'Color',plotcolors{z},'LineWidth',3);    
            end

            h = legend(legent);
            v = get(h,'title');
            set(v,'string','Likelihoods');

            for z = 1:length(likes)    
                patch([binndfminc{z}{q} fliplr(binndfminc{z}{q})],...
                      [binndpminc{z}{q}(2,:) fliplr(binndpminc{z}{q}(3,:))],...
                      plotcolors{z},'FaceAlpha',.25,'EdgeAlpha',0);

                %plot(binndfminc{z}{q},binndepminc{z}{q},'--','Color','k','LineWidth',2);
            end 
            
            catall = horzcat(binndpminc{:});
            catall2 = horzcat(catall{:});
            %ylim([0 ceil(max(catall2(:)))]);
            title(sprintf('Prior: %d',priors{session_nums(q)}.val),'FontSize',18);    
            
             case 'bin_nonoverlap'

            %subplot(1,length(session_nums),q); hold on;
            for z = 1:length(likes)
                plot(binndfminc{z}{q},binndpminc{z}{q}(1,:),'Color',plotcolors{z},'LineWidth',3);    
                %plot(fm_increasing{z},pm_increasing{z},'.','Color',plotcolors{z});
            end

            h = legend(legent);
            v = get(h,'title');
            set(v,'string','Likelihoods');

            for z = 1:length(likes)    
                patch([binndfminc{z}{q} fliplr(binndfminc{z}{q})],...
                      [binndpminc{z}{q}(2,:) fliplr(binndpminc{z}{q}(3,:))],...
                      plotcolors{z},'FaceAlpha',.25,'EdgeAlpha',0);

                plot(binndfminc{z}{q},zeros(size(binndfminc{z}{q})),'--','Color','k','LineWidth',2);
            end 
            
            catall = horzcat(binndpminc{:});
            catall2 = horzcat(catall{:});
            %ylim([0 ceil(max(catall2(:)))]);
            title(sprintf('Prior: %d',priors{session_nums(q)}.val),'FontSize',18);    
        end

    end
    
    FM{time_bins} = binndfminc;
    PM{time_bins} = binndpminc;
end


A = cellfun(@(x) x{1}(1,:),binndpminc,'UniformOutput',0);
ALo = cellfun(@(x) x{1}(2,:),binndpminc,'UniformOutput',0);
AHi = cellfun(@(x) x{1}(3,:),binndpminc,'UniformOutput',0);

average_width = vertcat(A{:}); 
average_widthH = vertcat(AHi{:});
average_widthL = vertcat(ALo{:});

A_width = average_width(:,[end 1]);
AH_width = average_widthH(:,[end 1]);
AL_width = average_widthL(:,[end 1]);
