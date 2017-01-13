function [tt,DB_inds,labels] = getTT_CK(bdf)
% getTT - returns a table containing the key timestamps 
%                   for all of the 1D Uncertainty trials in BDF
%
% Each row of the table corresponds to a single trial.  Columns are as
% follows:
%    1: Trial Databurst Timestamp
%    2: Cue 1 Location
%    3: Cue 2 Location
%    4: Trial Center Timestamp
%    5: Trial Center Hold Timestamp
%    6: Outer Target Timestamp
%    7: Trial Memory Timestamp
%    8: Trial Cue Delay Timestamp
%    9: Trial Go     Timestamp (i.e. Movement Start)  NaN if not R or F
%    10: Trial End    Timestamp
%    11: Color 1
%    12: Color 2
%    13: True location
%    14: Rate 1
%    15: Rate 2
%    16: Trial End  Result   -- R (32), A (33), F (34),
%                               I (35), or NaN

%%

%cd C:\Users\limblab\Desktop\s1_analysis\proc\bdek\Uncert\;

fn = bdf.meta.filename;
while ~isempty(fn)
    [tok fn]=strtok(fn,'\\');
end
fn=strtok(tok,'.');
fn = ['tt_' fn '.mat'];

%% Get the Databurst timestamps and pull the first timestamp
dtb = bdf.databursts;
burstlengths = cellfun(@(x) length(x),dtb(:,2));
dtb(burstlengths < 10,:) = [];
for dtbi=1:length(dtb)
    dtb_ts(dtbi) = dtb{dtbi,1};
    dtb_cue1(dtbi) = bytes2float(dtb{dtbi,2}(10:13));
    dtb_cue2(dtbi) = bytes2float(dtb{dtbi,2}(14:17));
    dtb_col1(dtbi) = bytes2float(dtb{dtbi,2}(18:21));
    dtb_col2(dtbi) = bytes2float(dtb{dtbi,2}(22:25));
    dtb_true(dtbi) = bytes2float(dtb{dtbi,2}(26:29));
    dtb_rate1(dtbi) = bytes2float(dtb{dtbi,2}(30:33));
    dtb_rate2(dtbi) = bytes2float(dtb{dtbi,2}(34:37));
    % Set NAN for erroneous perturbation databurst (dropped)
%     if abs(dtb_perts(dtbi))>50
%         dtb_perts(dtbi)=NaN;
%     end
end
% Only keep databursts that occur from 1 second and on
dtb_start_ind = find(dtb_ts>=1.0,1,'first');

%% Load and Process Words

% Find the first databurst time and truncate any words preceding it
words_first_ind = find(bdf.words(:,1) > dtb_ts(dtb_start_ind),1,'first');

%Find last trial end code index
words_end_ind=find(bdf.words(:,2)==hex2dec('20')|bdf.words(:,2)==hex2dec('21')|bdf.words(:,2)==hex2dec('22')|bdf.words(:,2)==hex2dec('23'),1,'last');
words_end_ts = bdf.words(words_end_ind,1);% last trial end time

% Truncate Words
wds = bdf.words(words_first_ind:words_end_ind,:);
wds_ts    = wds(:,1);
wds_codes = wds(:,2);

%% Process Full Databurst
% Truncate Databursts following the last end code time (if necessary)
dtb_last_ind = find(dtb_ts < words_end_ts,1,'last');
dtb_range = [dtb_start_ind:dtb_last_ind];

dtb_ts = dtb_ts(dtb_range);
dtb_cue1 = dtb_cue1(dtb_range);
dtb_cue2 = dtb_cue2(dtb_range);
dtb_col1 = dtb_col1(dtb_range);
dtb_col2 = dtb_col2(dtb_range);
dtb_true = dtb_true(dtb_range);
dtb_rate1 = dtb_rate1(dtb_range);
dtb_rate2 = dtb_rate2(dtb_range);

%% Process each trial
all_trial_word_inds = find(wds_codes>=32 & wds_codes<=34);
all_trial_word_ts = wds_ts(all_trial_word_inds);

% CENTER LOC   x30
% OUTER TARGET x40
% GO x31
% Reward       x20
% Abort        x21
% Failure      x22
% Incomplete   x23

center_code_inds     = find(wds_codes==hex2dec('30'));
center_hold_inds     = find(wds_codes==hex2dec('A0'));
outer_code_inds      = find(wds_codes==hex2dec('40'));
mem_delay_inds       = find(wds_codes==hex2dec('81'));
targcue_inds         = find(wds_codes==hex2dec('82'));
go_code_inds         = find(wds_codes==hex2dec('31'));
outer_hold_inds      = find(wds_codes==hex2dec('A1'));

reward_code_inds     = find(wds_codes==hex2dec('20'));
abort_code_inds      = find(wds_codes==hex2dec('21'));
failure_code_inds    = find(wds_codes==hex2dec('22'));
incomplete_code_inds = find(wds_codes==hex2dec('23'));

center_ts       = wds_ts(center_code_inds);
cthold_ts       = wds_ts(center_hold_inds);
outer_ts        = wds_ts(outer_code_inds);
memdelay_ts     = wds_ts(mem_delay_inds);
targcue_ts      = wds_ts(targcue_inds);
go_ts           = wds_ts(go_code_inds);
othold_ts       = wds_ts(outer_hold_inds);
reward_ts       = wds_ts(reward_code_inds);
failure_ts      = wds_ts(failure_code_inds);
abort_ts        = wds_ts(abort_code_inds);
incomplete_ts   = wds_ts(incomplete_code_inds);

% GETUNTRIALTABLE1D - returns a table containing the key timestamps 
%                   for all of the 1D Uncertainty trials in BDF
%
% Each row of the table corresponds to a single trial.  Columns are as
% follows:
%    1: Trial Databurst Timestamp
%    2: Cue 1 Location
%    3: Cue 2 Location
%    4: Trial Center Timestamp
%    5: Trial Center Hold Timestamp
%    6: Outer Target Timestamp
%    7: Trial Memory Timestamp
%    8: Trial Cue Delay Timestamp
%    9: Trial Go     Timestamp (i.e. Movement Start)  NaN if not R or F
%    10: Trial End    Timestamp
%    11: Color 1
%    12: Color 2
%    13: True location
%    14: Rate 1
%    15: Rate 2
%    16: Trial End  Result   -- R (32), A (33), F (34),
%                               I (35), or NaN

dropTrial = false;
trial_counter=1;
DB_inds = nan(length(dtb_range),1);
% For each trial complete code
for trial_i=1:length(all_trial_word_inds)
    
    %Output 1 and 2
    dtb_ind = find(dtb_ts < all_trial_word_ts(trial_i),1,'last');
    
    if ~isempty(dtb_ind) && ~ismember(dtb_ind,DB_inds)
        out_dtb_ts = dtb_ts(dtb_ind);
        out_dtb_cue1 = dtb_cue1(dtb_ind);
        out_dtb_cue2 = dtb_cue2(dtb_ind);
        out_dtb_col1 = dtb_col1(dtb_ind);
        out_dtb_col2 = dtb_col2(dtb_ind);
        out_dtb_true = dtb_true(dtb_ind);
        out_dtb_rate1 = dtb_rate1(dtb_ind);
        out_dtb_rate2 = dtb_rate2(dtb_ind);
        
        out_dtb_ind = dtb_ind;
  
    else
        out_dtb_ts = NaN;
        out_dtb_cue1 = NaN;
        out_dtb_cue2 = NaN;
        out_dtb_col1 = NaN;
        out_dtb_col2 = NaN;
        out_dtb_true = NaN;
        out_dtb_rate1 = NaN;
        out_dtb_rate2 = NaN;
        out_dtb_ind = NaN;
    end
    
    out_end_ts = floor(1000*wds_ts(all_trial_word_inds(trial_i)))/1000;
    result_code = wds_codes(all_trial_word_inds(trial_i));
    if result_code==hex2dec('20')
        out_result = result_code;
        out_go_ts = go_ts(find(go_ts<out_end_ts,1,'last'));
        out_targ_ts = outer_ts(find(outer_ts<out_end_ts,1,'last'));
        
    elseif result_code==hex2dec('22')
        out_result = result_code;     
        out_go_ts = go_ts(find(go_ts<out_end_ts,1,'last'));        
        out_targ_ts = outer_ts(find(outer_ts<out_end_ts,1,'last'));    
    elseif result_code==hex2dec('21')
        out_result = result_code;
        out_go_ts = NaN;
        out_targ_ts = NaN;
    elseif result_code==hex2dec('23')
        out_result = result_code;
        out_go_ts = NaN;
        out_targ_ts = NaN;
    else
        out_result = NaN;
        out_go_ts = NaN;
        out_targ_ts = NaN; 
    end
    if isempty(out_go_ts)
        out_go_ts=NaN;
    end
    if isempty(out_targ_ts)
    	out_targ_ts=NaN;
    end
    if out_targ_ts < out_dtb_ts || out_targ_ts < 1
        out_targ_ts = NaN;
    end
    if out_go_ts < out_dtb_ts || out_go_ts < 1
        out_go_ts = NaN;
    end
    if out_end_ts < out_dtb_ts || out_end_ts < 1
        out_end_ts = NaN;
    end
    %Output 3
    cn_ind = find(center_ts<out_end_ts,1,'last');
    try
        cn_ind2 = find(center_ts>dtb_ts(out_dtb_ind),1,'first');
    catch err
        cn_ind2 = [];
    end
    if ~isempty(cn_ind) && ~isempty(cn_ind2)
      
        if cn_ind==cn_ind2
            out_center_ts = center_ts(cn_ind);
        elseif (center_ts(cn_ind) < dtb_ts(dtb_ind)) && (center_ts(cn_ind2) < out_end_ts)
            out_center_ts = center_ts(cn_ind2);
        elseif (center_ts(cn_ind) > dtb_ts(dtb_ind)) && (center_ts(cn_ind2) > out_end_ts)
            out_center_ts = center_ts(cn_ind);
        else
            out_center_ts = NaN;
        end
    else
        out_center_ts = NaN;
    end
    
    % Later outputs
    if ~isnan(out_center_ts)
        out_cthold_ts = cthold_ts(find(cthold_ts>out_center_ts,1,'first'));
        out_memdelay_ts = memdelay_ts(find(memdelay_ts>out_center_ts,1,'first'));
        out_targcue_ts = targcue_ts(find(targcue_ts>out_center_ts,1,'first'));
        
        if out_targcue_ts > out_end_ts
            out_targcue_ts = -1;
        end
        
    else
        out_cthold_ts = NaN;
        out_memdelay_ts = NaN;
        out_targcue_ts = NaN;
    end
    
    if isempty(out_memdelay_ts); out_memdelay_ts = NaN; end
    if isempty(out_targcue_ts); out_targcue_ts = NaN; end
 

    if 1%~dropTrial
        DB_inds(trial_counter) = out_dtb_ind;
        tt(trial_counter,:) = [...
            out_dtb_ts, ...
            out_dtb_cue1, ...
            out_dtb_cue2, ...
            out_center_ts, ...
            out_cthold_ts, ...
            out_targ_ts, ...
            out_memdelay_ts,...
            out_targcue_ts,...
            out_go_ts, ...
            out_end_ts, ...
            out_dtb_col1,... 
            out_dtb_col2,... 
            out_dtb_true,... 
            out_dtb_rate1,... 
            out_dtb_rate2,...
            out_result];  % Result of trial ('R', 'A', 'I', or 'N')
       
        trial_counter=trial_counter+1; 
    end
    dropTrial=false;
end

bad_trials = sum(isnan(tt(:,4:7)),2)>=1;
tt(bad_trials,:) = [];

labels = {'1: DB time';'2: Cue 1 Loc';'3: Cue 2 Loc';'4: CT on';'5: CT hold';'6: OTs on';'7: Memory';...
    '8: CT cue';'9: GO';'10: Trial End';'11: Color 1';'12: Color 2';'13: True target Loc';'14: Rate 1';'15: Rate 2';...
    '16: End Result'};

return