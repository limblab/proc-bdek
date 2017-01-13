function [corrected_ts, lags] = photo_correct(tt,bdf,columns,fail_type)

% [corrected_ts,photo] = photo_correct(tt,bdf,columns) 
% This function uses photodetector timestamps to replace system commands to
% correct for monitor lag.


%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%          tt: Trial Table
%         bdf: bdf with photodector onset pulses in AI15 of NSP (id[143 1])
%     columns: Indices of desired columns in trial table 
%   fail_type: Dictates action when no appropriate photodetector timestamp
%              is found
%             
%             'nan' - Return NaN
%             'mean' - Use average lag computed from all other trials    
%             'command' - Leave command signal timestamp 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >> [corrected_ts , lags] = photo_correct(tt,bdf,[5 6 7],'mean');
% >> new_tt = tt;
% >> new_tt(:,[5 6 7]) = corrected_ts;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
photo_ts = sortrows(vertcat(bdf.units(vertcat(bdf.units.id)==143).ts));

% Initialize
corrected_ts = zeros(size(tt,1),length(columns)); 
photo = zeros(size(tt,1),1);

for j = 1:length(columns) % Loop through trial table columns
    for i = 1:size(tt,1) % Loop through trial table rows
        
        % Find system command ts
        comm_ts = tt(i,columns(j)); 
        if isnan(comm_ts)
            comm_ts = max(tt(:,columns(j))); % Guarantee fail state for NaN
        end
        
        % Find photodetector ts immediately following system command ts
        nextind = find(photo_ts > comm_ts,1,'first');
        if isempty(nextind); 
            photo(i) = comm_ts + inf; 
        else
            photo(i) = photo_ts(nextind);
        end
    
        if comm_ts < 0
            corrected_ts(i,j) = comm_ts;
        % Check for failure state (lag greater than 200 ms)
        elseif (photo(i) - comm_ts) > 0.2   
            
            % Do appropriate failure state action
            if strcmp(fail_type,'nan')
                corrected_ts(i,j) = NaN;   
            elseif strcmp(fail_type,'mean')
                corrected_ts(i,j) = -1;
            elseif strcmp(fail_type,'command')
                corrected_ts(i,j) = tt(i,j); 
            end
            
        else   % Use photodetector timestamp 
            corrected_ts(i,j) = photo(i);
        end
    end
end

% Find average lag
lags = corrected_ts - tt(:,columns);
mean_lag = nanmean(lags(lags > 0));

tt_in_use = tt(:,columns);

if strcmp (fail_type,'mean')
    % In case of fail type 'mean', replace ts with mean lag-calculated ts. 
    corrected_ts(corrected_ts == -1) = tt_in_use(corrected_ts == -1) + mean_lag;
    lags(lags < 0) = nan; 
end

end

