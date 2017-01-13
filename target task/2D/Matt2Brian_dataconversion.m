% Important columns in trial table
%    3: trial type - use as marker to split trials (e.g. early/late)
%    5: target on time 
%    6: go cue time
%    7: movement onset time
%   10: trial direction of interest [rad] (target location or movement dir)

tt = zeros(size(BL.trial_table,1),11);
tt(:,[10 5 6 7]) = BL.trial_table(:,[1 2 3 4]);
tt(:,3) = 1;
timeoffset = tt(end,7)+1000; % For creating 1 timescale

tt2 = zeros(size(AD.trial_table,1),11);
tt2(:,10) = AD.trial_table(:,1);
tt2(:,[5 6 7]) = AD.trial_table(:,[2 3 4])+timeoffset;
tt2(:,3) = 2;

alldays(1).tt = tt;
alldays(2).tt = [tt; tt2];

% Re-format M1
for i = 1:length(BL.M1.units); 
    alldays(1).M1_units{i} = [BL.M1.units(i).ts'; AD.M1.units(i).ts'+timeoffset];
end

% Re-format PMd
for i = 1:length(BL.PMd.units); 
    alldays(1).PMd_units{i} = [BL.PMd.units(i).ts'; AD.PMd.units(i).ts'+timeoffset];
end
