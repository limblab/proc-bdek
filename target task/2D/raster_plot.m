function[raster,allinds,fig_hand,SpikeTime] = raster_plot(train,trial_table,t1_t2,alignment,do_plot,position_rank,TitleString,varargin)
%% [raster,allinds,fig_hand,SpikeTime] = raster_plot(train,trial_table,t1_t2,alignment,do_plot,position_rank,TitleString,varargin)
if nargin < 7
    TitleString = '';
end

t1 = t1_t2(1);
t2 = t1_t2(2);

if strcmp(alignment,'center')
    aligntype = 4;
elseif strcmp(alignment,'target')
    aligntype = 5; % 5 for target, 6 for go cue
elseif strcmp(alignment,'go')
    aligntype = 6;
elseif strcmp(alignment,'end')
    aligntype = 7;
elseif isnumeric(alignment)
    aligntype = alignment;
else
    fprintf('choose ''target'' or ''go''\n');
end

if do_plot==1; Plotraster = 1; else Plotraster = 0; end% 1 to plot raster

if strcmp(position_rank,'endpoint')
    rank_anchor = 1; 
elseif strcmp(position_rank,'target')
    rank_anchor = 2;
elseif strcmp(position_rank,'trials')
    rank_anchor = 3;
elseif strcmp(position_rank,'none')
    rank_anchor = 0;
elseif strcmp(position_rank,'uncertainty');
    rank_anchor = 4;
end

% 1: end position
% 2: target position

% Pick out rows from trial table for low and high uncertainty

all_likes = unique(trial_table(:,3));
allinds = cell(length(all_likes),1);
allvartrials = cell(length(all_likes),1);
[LA] = deal(cell(length(all_likes),1));
for i = 1:length(all_likes)
    allinds{i} = find(trial_table(:,3)==all_likes(i));
    allvartrials{i} = [allinds{i} trial_table(allinds{i},:)];
    LA{i} = size(allvartrials{i},1);
end
aligntype = aligntype + 1; % Increment to account for appended column to allvartrials

PlotColors = {'r.','g.','b.'};

if length(all_likes) > 1
    LA_prev = cellfun(@(x) x,LA);
    LA_prev = [0; LA_prev(1:end-1)];
else
    LA_prev = 0;
end
%% Trials
if Plotraster == 1
	fig_hand = figure; hold on;
else
    fig_hand = 0;
end

[raster, final_pos, SpikeTime] = deal(cell(length(all_likes),1));
for q = 1:length(all_likes)
    
    raster{q} = nan(LA{q},round(1000*(t2-t1))); % Raster with NaN
    final_pos{q} = NaN(LA{q},1);
    for i = 1:LA{q}

        % Get Time of start
        timestart = allvartrials{q}(i,aligntype)+t1;

        % Get Endpoint/target pos
        if rank_anchor == 1
            final_pos{q}(i) = allvartrials{q}(i,11);
        elseif rank_anchor == 2
            final_pos{q}(i) = allvartrials{q}(i,3);
        elseif rank_anchor == 0
            final_pos{q}(i) = 0;
        end

        % Align spikes to start and get rid of those not in region
        aligned_ts = round(1000*(train - timestart));
        aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
        
        SpikeTime{q}(i,:) = mean(aligned_ts); 

        % Fill out rasters
        if strcmp(position_rank,'trials')
            raster{q}(i,aligned_ts) = i;
        elseif strcmp(position_rank,'none')
            raster{q}(i,aligned_ts) = 1;
        elseif strcmp(position_rank,'uncertainty')
            raster{q}(i,aligned_ts) = i + sum(LA_prev(1:q));
        else
            raster{q}(i,aligned_ts) = final_pos{q}(i);
        end

    end

%% Do Plots
    if do_plot ~= 0
        xplotter = repmat((1000*t1):(1000*t2-1),LA{q},1);
        plot(xplotter,raster{q},PlotColors{q},'MarkerSize',15);
    end
end

allrast = vertcat(raster{:});

if do_plot ~= 0

    xlabel(sprintf('Time from %s (ms)',alignment),'FontSize',14);
    if strcmp(position_rank,'none')
      ylabel('Trial Number','FontSize',14);
    else 
      ylabel(sprintf('%s location (rad)',position_rank),'FontSize',14);
    end
    title(TitleString,'FontSize',16);

    if (min(allrast(:)) ~= max(allrast(:)) && ~isnan(min(allrast(:))) && ~isnan(min(allrast(:))))
        ylim([min(allrast(:)) max(allrast(:))]);
    end

end

