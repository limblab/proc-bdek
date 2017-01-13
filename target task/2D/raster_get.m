function[raster] = raster_get(train,trial_table,t1_t2,alignment,upper_time_lims,varargin)

t1 = t1_t2(1); % in seconds
t2 = t1_t2(2); % in seconds

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

raster = zeros(size(trial_table,1),round(1000*(t2-t1))); % Raster with NaN
for i = 1:size(trial_table,1)

    % Get Time of start
    timestart = trial_table(i,aligntype)+t1;

    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
    
    raster(i,aligned_ts) = 1; 
    if nargin > 4
        upperlimind = round(1000*(upper_time_lims(i)-timestart));
        indnums = 1:size(raster,2);
        raster(i,indnums > upperlimind) = NaN;
    end
end




end

