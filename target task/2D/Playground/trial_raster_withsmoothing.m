function[trial_raster,dat] = trial_raster_withsmoothing(units,trial_table,t1_pm,t2_pm,maxISI,causal,subsamplesize)
% [trial_raster,dat] = trial_raster(units,trial_table,t1_pm,t2_pm) is a
% function that creates a raster (trial_raster) in cell array form. In
% addition, the second output 'dat' is a struct containing trial number 
% and spike information for use with the GPFA library from Byron Yu
% (users.ece.cmu.edu/~byronyu/software.shtml)
% INPUTS
%   units: This is a cell array in which each cell contains the spike times
%          for a single neuron
%   trial_table: Trial table, in which at least one column contains
%          timestamps of an event on which you want to align the raster
%   t1_pm: this is a vector with size 1x2. First element is the column for
%          which you want to align the start of the raster. The second
%          element is an offset (in milliseconds) from that column.
%              example: t1_pm = [5 200] means start the raster 200 
%                       milliseconds after the time in column 5. 
%   t2_pm: Same as t2_pm but for the raster end time. 
%   
%   NOTE: the raster for each trial does NOT need to be of constant length.
%         (i.e. trial_table(i,t2_pm(1)) - trial_table(i,t1_pm(1)) can be
%         different for every trial i)
%   
% maxISI: represents two standard deviations of the Gaussian kernel (in ms) 
%        Example: maxISI = 50 means that the amplitude of the continuous
%        signal will only significantly increase when an ISI is less than
%        50 ms. 
% causal: If (1), then maxISI is doubled and only the causal half of the
% gaussian is used. 
sd = maxISI/2;
N = 6*sd+1;
if causal
    filtprof = fspecial('gaussian',[2*N 1],2*sd)*2;
    filtprof(1:floor(end/2)) = 0; 
else
    filtprof = fspecial('gaussian',[2*N 1],2*sd);
end

align1 = t1_pm(1); pm1 = t1_pm(2); % Break out the aligment values for clarity
align2 = t2_pm(1); pm2 = t2_pm(2);

% Create a cell array in which each cell contains information from one trial
trial_raster = cell(size(trial_table,1),1); 
dat = struct;
% Loop through trials
for i = 1:size(trial_table,1);
    
    if ~isnan(sum(trial_table(i,[align1,align2])))
        tlen = round(1000*((trial_table(i,align2)+(pm2./1000))-(trial_table(i,align1)+(pm1./1000))));

        t1 = trial_table(i,align1)+(pm1./1000)-3*maxISI/1000; %start time
        t2 = trial_table(i,align2)+(pm2./1000)+3*maxISI/1000; %end time

        % Initialize raster (NxT) where N is number of neurons and T is time
        raster = zeros(length(units),round(tlen/subsamplesize)); 

        for neur = 1:length(units)
            temp_rast = zeros(1,round(1000*(t2-t1)));
            % Align spikes to start and get rid of those not in region
            aligned_ts = round(1000*(units{neur}(2:end) - t1));
            aligned_ts(aligned_ts<=0 | aligned_ts>=((t2-t1)*1000)) = [];

            % Fill out rasters
            temp_rast(aligned_ts) = 1;
            cont_signal = conv(temp_rast,filtprof,'same')*1000;
            cont_signal = cont_signal((3*maxISI+1):(end-3*maxISI));

            raster(neur,:) = bin_array(cont_signal,1,round(tlen/subsamplesize),'mean');
        end

        trial_raster{i} = raster; 
        dat(i).spikes = raster;
    else
        trial_raster{i} = NaN(length(units),1);
        dat(i).spikes = NaN(length(units),1);
    end
    dat(i).trialId = i;
%     fprintf('trial: %d/%d\n',i,size(trial_table,1));
end