brain_area = 'PMd';
loop_alignments = {'target','go'};
 loop_ranges = {-200:100:1000 , -100:100:500};
% loop_ranges = {-200:200:1000 , -200:200:600};

prediction_day_indices = 2;

comb_tt = [alldays(1).tt; alldays(prediction_day_indices).tt];

figure; hold on;
for loopthrough = 1:length(loop_alignments)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tune_ranges = loop_ranges{loopthrough};
    tune_align = loop_alignments{loopthrough};

    %
    time_bins = tune_ranges;
    time_align = tune_align;

    [unit_counts,mintrial,maxtrial] = deal(cell(length(prediction_day_indices),1));
    for z = 1:length(prediction_day_indices)

        prediction_day = prediction_day_indices(z);
        day_units = eval(sprintf('alldays(1).%s_units',brain_area));
        if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

        day_pred = prediction_day;
        tt_pred = comb_tt;

        % find rasters for each unit
        unit_pred = cell(length(day_units),1);
        for q = 1:length(day_units)

            clc;
            fprintf('Day: %d/%d\nUnit: %d/%d\n',z,length(prediction_day_indices),q,length(day_units));

            [rast_out,rast_inds] = raster_plot(day_units{q},tt_pred,[time_bins(1)./1000 time_bins(end)./1000],...
                time_align,0,'none');

            out = nan(size(vertcat(rast_out{:})));
            for m = 1:length(rast_out)
                out(rast_inds{m},:) = rast_out{m};
            end
            out(isnan(out))=0;

            rast = zeros(size(out,1),length(time_bins)-1);
            binsizes = time_bins - time_bins(1); binsizes(1) = 1;
            for v = 1:length(time_bins)-1
                rast(:,v) = sum(out(:,binsizes(v):binsizes(v+1)),2);
            end
            unit_counts{z}{q} = rast./repmat(diff(time_bins)/1000,size(rast,1),1);
        end
    end
    %%
    [trial_counts,firing_absolute,spike_count_raw] = deal(cell(length(prediction_day_indices),1));
    for bin = 1:length(time_bins) - 1
        clc; fprintf('bin %d/%d\n',bin,length(time_bins)-1);
        for z = 1:length(prediction_day_indices)

            day_pred = prediction_day_indices(z);
            tt_pred = comb_tt;  

            day_units = eval(sprintf('alldays(1).%s_units',brain_area));
            if strcmp(brain_area,'M1'); day_units(end-1:end) =[]; end

            for i = 1:size(tt_pred,1)   
                for q = 1:length(day_units)       
                    trial_counts{z}{bin}(i,q) = unit_counts{z}{q}(i,bin);
                end
            end
            firing_absolute{z}{bin} = trial_counts{z}{bin};
            spike_count_raw{z}{bin} = trial_counts{z}{bin}*(diff(time_bins(1:2))/1000);
        end
    end

    xsforplot = .5*(time_bins(1:end-1)+time_bins(2:end));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    fanofunc = @(x) var(x)./mean(x);
    crossfunc = @(x,y) cos(x).*sin(y) - sin(x).*cos(y);                     
    inbinfunc = @(x,b1,b2) crossfunc(x,b1) < 0 & crossfunc(x,b2) > 0;
    cobins = (0:pi/4:2*pi) + pi/8;

    FRthresh = 0;
    NEL = 10;

    TT = comb_tt;
    priord = circ_mean(TT(:,2));

    LHvals = flipud(unique(TT(:,3)));
    movedirs = TT(:,10);
    [FANOS,FF] = deal(cell(length(LHvals),1));
    [BF] = deal(zeros(length(LHvals),3));

    if sum(LHvals > 1000)
        colp = {'k','b','g','r'};
    else
        colp = {'b','g','r'};
    end

    for TIMES = 1:length(spike_count_raw{1})
        timebinindex = TIMES;
        COUNTS = spike_count_raw{1}{timebinindex};%.*ODinds;

        [LBIN_counts,FANOS,ALLFANOS] = deal(cell(length(LHvals),1));
        for likei = 1:length(LHvals)

            % Find indices for current likelihood condition
            likecondinds = find(TT(:,3)==LHvals(likei));
            COUNTS_like = COUNTS(likecondinds,:);
            TT_like = TT(likecondinds,:);

            if LHvals(likei) > 1000 % If we're in center out
                for cobin = 1:(length(cobins)-1) % loop through CO targets

                    % Find CO trials within current target
                    incobin = find(inbinfunc(TT_like(:,10),cobins(cobin),cobins(cobin+1)));
                    TT_cobin = TT_like(incobin,:);
                    COUNTS_cobin = COUNTS_like(incobin,:);

                    % Sort the reaches by distance to center of the target
                    [comoves_ordered,coorderind] = ...
                        sortrows(circ_dist(TT_cobin(:,10),cobins(cobin)+pi/8));

                    % Separate into groups of desired trial length
                    coorder_pad = [coorderind; nan(NEL-rem(length(coorderind),NEL),1)];
                    cogroups = reshape(coorder_pad,NEL,[]);
                    cogroups(:,isnan(sum(cogroups))) = [];
                    for ci = 1:size(cogroups,2) % Loop through n-trial groups
                        cogroupmoves = TT_cobin(cogroups(:,ci),10);
                        cogroupspan = abs(circ_dist(cogroupmoves(1),cogroupmoves(end)));

                        if cogroupspan < pi/4
                            % Assign spike counts
                            LBIN_counts{likei}{cobin,ci} = COUNTS_cobin(cogroups(:,ci),:);
                        else
                            LBIN_counts{likei}{cobin,ci} = NaN;
                        end
                    end

                end
            else

                [moves_ordered,orderind] = sortrows(circ_dist(TT_like(:,10),priord));
                order_pad = [orderind; nan(NEL-rem(length(orderind),NEL),1)];
                groups = reshape(order_pad,NEL,[]);
                groups(:,isnan(sum(groups))) = [];

                for i = 1:size(groups,2)
                    groupmoves = TT_like(groups(:,i),10); 
                    groupspan = abs(circ_dist(groupmoves(1),groupmoves(end)));

                    if groupspan < pi/4
                        LBIN_counts{likei}{i} = COUNTS_like(groups(:,i),:);
                    else
                        LBIN_counts{likei}{i} = NaN;
                    end
                end
            end

            FANOS{likei} = cellfun(fanofunc,LBIN_counts{likei},'UniformOutput',0);
            ALLFANOS{likei} = horzcat(FANOS{likei}{:});

            [BF(likei,1),BF(likei,2)] = boot_bounds(1000,@nanmean,ALLFANOS{likei},2.5,97.5);
            BF(likei,3) = nanmean(ALLFANOS{likei});

            FF{likei}(:,TIMES) = BF(likei,[1 2 3])';

        end

    end
    %%
 
    if loopthrough == 1

        xsforplotT = .5*(time_bins(1:end-1)+time_bins(2:end));
        max_targ_x = max(xsforplotT);
        for i = 1:length(LHvals)
            plot(xsforplotT,FF{i}(3,:),colp{i});
            patch([xsforplotT fliplr(xsforplotT)],[FF{i}(1,:) fliplr(FF{i}(2,:))],...
                colp{i},'FaceAlpha',0.25,'EdgeAlpha',0.25);
        end

    else

        xsatgo = .5*(time_bins(1:end-1)+time_bins(2:end));
        xsforplotG = xsatgo + max_targ_x + 100;
        min_go_x = min(xsatgo);

        for i = 1:length(LHvals)
            plot(xsforplotG,FF{i}(3,:),colp{i});
            patch([xsforplotG fliplr(xsforplotG)],[FF{i}(1,:) fliplr(FF{i}(2,:))],...
                colp{i},'FaceAlpha',0.25,'EdgeAlpha',0.25);
        end 
    
    end
end
