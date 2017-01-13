function rasterplot(bdf,train,trial_table,t1,t2)
%%%
% OUTPUTS = [meanFR, final_pos, target_pos, maxspeed, maxspeed_time, direction, UNC];

aligntype = 5; %5=target, 6=go;
colorson = 1;


raster = nan(length(trial_table),round(1000*(t2-t1)));

plottypes = cell(2,1);
plottypes{1} = 'b';
plottypes{2} = 'r';
plottypes{3} = 'k';

figure; hold on;
for i = 1:length(trial_table)
    
    if colorson == 1
        if trial_table(i,end)== min(trial_table(:,end))
            trialtype = 1;
        else 
            trialtype = 2;
        end
    else
        trialtype = 3;
    end

    % Get Time of start
    timestart = trial_table(i,aligntype)+t1;
    bdfstart = find(bdf.pos(:,1)<timestart,1,'last');
    bdfgo = find(bdf.pos(:,1)<trial_table(i,6),1,'last');
    bdfend = bdfstart + 1000*(t2-t1)-1;
    bdfrew = find(bdf.pos(:,1)>trial_table(i,7),1,'first');
    
    % Align spikes to start and get rid of those not in region
    aligned_ts = round(1000*(train - timestart));
    aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];
    
    % Fill out rasters
    raster(i,aligned_ts) = 1;
        
    num_spikes(i,:) = [trialtype nansum(raster(i,:))];

    plot((1000*t1):1000*t2-1,i*raster(i,:),[plottypes{trialtype} '.']);

end
% figure; hold on;
% for i =1:length(trial_table)
%     bar(i,num_spikes(i,2),plottypes{num_spikes(i,1)});
% end

