function [tt,priors,slices,bump,fulls,labels,posoff] = getTT_CC(bdf)
% getTT - returns a table containing the key timestamps 
%                   for all of the 1D Uncertainty trials in BDF
%
% Each row of the table corresponds to a single trial.  Columns are as
% follows:
labels = {'(1)  DB timestamp';...       
          '(2)  Target Location';...   
          '(3)  Visual Uncertainty';...
          '(4)  Bump Magnitude';...
          '(5)  Stimulate';...
          '(6)  Center TS';...
          '(7)  Bump TS';...
          '(8)  Go TS';...
          '(9)  End TS';...
          '(10) End Result';...
          '(11) Centroid';...
          '(12) Endpoint'};

%%
% Truncate Words
wds = bdf.words;
wds_ts    = wds(:,1);
wds_codes = wds(:,2);

% CENTER LOC   x30
% OUTER TARGET x40
% GO           x31
% Reward       x20
% Abort        x21
% Failure      x22
% Incomplete   x23

center_ts       = wds_ts(wds_codes==hex2dec('30'));
%outer_ts        = wds_ts(wds_codes==hex2dec('40'));
bump_ts         = wds_ts(wds_codes==hex2dec('50'));
go_ts           = wds_ts(wds_codes==hex2dec('31'));
%reward_ts       = wds_ts(wds_codes==hex2dec('20'));
%abort_ts        = wds_ts(wds_codes==hex2dec('21'));
%failure_ts      = wds_ts(wds_codes==hex2dec('22'));
%incomplete_ts   = wds_ts(wds_codes==hex2dec('23'));
complete_ts     = wds_ts(ismember(wds_codes,hex2dec({'20','21','22','23'})));
results         = wds_codes(ismember(wds_codes,hex2dec({'20','21','22','23'})));

%% Get the Databurst timestamps and pull the first timestamp
dtb = bdf.databursts;
burstlengths = cellfun(@(x) length(x),dtb(:,2));
dtb(burstlengths < 10,:) = [];
dtb = [dtb; {inf,[]}];

dtb_ts = vertcat(dtb{:,1});
[dtb_targ,prior_list,dtb_bumpmag,dtb_bumpdur,dtb_backdur,dtb_stim,...
    dtb_visual_kap,dtb_slicenum,dtb_center_ts,dtb_bump_ts,dtb_go_ts,...
    dtb_complete_ts,dtb_result,dtb_centroid] = deal(nan(length(dtb)-1,1));
dtb_slices = nan(length(dtb),10);
for dtbi=1:length(dtb)-1
    dtb_targ(dtbi)      = bytes2float(dtb{dtbi,2}(10:13));
    prior_list(dtbi)    = bytes2float(dtb{dtbi,2}(14:17));
    dtb_bumpmag(dtbi)   = bytes2float(dtb{dtbi,2}(18:21));
    dtb_bumpdur(dtbi)   = bytes2float(dtb{dtbi,2}(22:25));
    dtb_backdur(dtbi)   = bytes2float(dtb{dtbi,2}(26:29));
    dtb_stim(dtbi)      = dtb{dtbi,2}(30);
    dtb_visual_kap(dtbi)= bytes2float(dtb{dtbi,2}(31:34));
    dtb_slicenum(dtbi)  = bytes2float(dtb{dtbi,2}(35:38));
    dtb_slices(dtbi,:)  = bytes2float(dtb{dtbi,2}(39:end))';
    dtb_centroid(dtbi)  = circ_mean(dtb_slices(dtbi,:),[],2);
    
    % Timestamps
    center_tsi = find(center_ts > dtb_ts(dtbi) & center_ts < dtb_ts(dtbi+1));
    bump_tsi = find(bump_ts > dtb_ts(dtbi) & bump_ts < dtb_ts(dtbi+1));
    go_tsi = find(go_ts > dtb_ts(dtbi) & go_ts < dtb_ts(dtbi+1));
    complete_tsi = find(complete_ts > dtb_ts(dtbi) & complete_ts < dtb_ts(dtbi+1));
    
    if isempty(center_tsi); dtb_center_ts(dtbi)=NaN; else dtb_center_ts(dtbi) = center_ts(center_tsi(1)); end
    if isempty(bump_tsi); dtb_bump_ts(dtbi)=NaN; else dtb_bump_ts(dtbi) = bump_ts(bump_tsi(1)); end
    if isempty(go_tsi); dtb_go_ts(dtbi)=NaN; else dtb_go_ts(dtbi) = go_ts(go_tsi(1)); end
    if isempty(complete_tsi) 
        dtb_complete_ts(dtbi)=NaN; 
        dtb_result(dtbi)=NaN;    
    else
        dtb_complete_ts(dtbi) = complete_ts(complete_tsi(1));
        dtb_result(dtbi) = results(complete_tsi(1));
    end
end
slicenum = mode(dtb_slicenum);
bumpdur = mode(dtb_bumpdur);
backdur = mode(dtb_backdur);

% Ditch trials for which the databurst occured before 1s
badfirst = find(dtb_ts < 1);
%% TT
tt_all=[...
    dtb_ts(1:end-1), ...
    dtb_targ, ...
    dtb_visual_kap, ...
    dtb_bumpmag,...
    dtb_stim,...
    dtb_center_ts, ...
    dtb_bump_ts, ...
    dtb_go_ts, ...
    dtb_complete_ts, ...
    dtb_result, ... % Result of trial ('R', 'A', 'I', or 'N')
    dtb_centroid]; 
% Only keep trials that occur from 1 second and on
tt_all(badfirst,:) = [];

%% Slices
slices_all = dtb_slices(:,1:slicenum);
slices_all(badfirst,:) = [];

%% Prior 
prior_list(badfirst) = [];
% Deal with weird prior databursts
check_bb = @(burst) burst < 10e-5 | burst > 1e5+1;
bad_bursts = find(check_bb(prior_list));
good_bursts = find(~check_bb(prior_list));
for i = 1:length(bad_bursts)
    bb = bad_bursts(i);
    ind_dists = abs(good_bursts - bb);
    replacer_ind = good_bursts(find(ind_dists==min(ind_dists),1,'first'));
    replacer = prior_list(replacer_ind);
    prior_list(bb)= replacer;
end
%%
BDF = bdf;
BDF.pos(:,2) = BDF.pos(:,2) - 3;
BDF.pos(:,3) = BDF.pos(:,3) + 33;

[endangle,endangle_curs] = deal(NaN(size(tt_all,1),1));
for i = 1:size(tt_all,1)
    locatbump_i = find(bdf.pos(:,1) > (tt_all(i,7)+bumpdur),1,'first');
    if ~isnan(locatbump_i)
        posoff(i,:) = BDF.pos(locatbump_i,2:3);

        pos_o = BDF.pos;
        pos_o(:,2:3) = pos_o(:,2:3)-repmat(posoff(i,:),size(BDF.pos,1),1);
    else
        posoff(i,:) = nan;
        pos_o = NaN;
    end
    
    %bdf_end_i = find(bdf.pos(:,1) > tt_all(i,9),1,'first');
    bdf_end_i = find(pos_o(:,1) > tt_all(i,9),1,'first');
    if ~isempty(bdf_end_i)
        %endpoint = BDF.pos(bdf_end_i,2:3);
        endpoint = pos_o(bdf_end_i,2:3);
        endpoint_hand = BDF.pos(bdf_end_i,2:3);
        endangle(i) = atan2(endpoint(2),endpoint(1));
        endangle_curs(i) = atan2(endpoint_hand(2),endpoint_hand(1));
    else
        endangle(i) = NaN;
        endangle_curs(i) = NaN;
    end
end

tt_all = [tt_all, endangle, endangle_curs];

% Eliminate junk
bad_trials = [find(isnan(sum(tt_all,2)));...
              find(sum(slices_all,2)>10000 | sum(slices_all,2)<-10000)];

tt_all(bad_trials,:) = [];
slices_all(bad_trials,:) = [];
prior_list(bad_trials) = [];

fulls.tt = tt_all;
fulls.slices = slices_all;
fulls.priors = prior_list;

completed = find(ismember(tt_all(:,10),[32 34]));

tt_all = tt_all(completed,:);
slices_all = slices_all(completed,:);
prior_list = prior_list(completed);

%% Separate into prior blocks 
prior_switches = [0 find(diff(prior_list)~=0) length(prior_list)];
[tt,priors,slices] = deal(cell(length(prior_switches)-1,1));
for i = 1:length(prior_switches)-1
    tt{i} = tt_all((prior_switches(i)+1):prior_switches(i+1),:);
    priors{i}.inds = (prior_switches(i)+1):prior_switches(i+1);
    priors{i}.val = prior_list(prior_switches(i)+1);
    slices{i} = slices_all((prior_switches(i)+1):prior_switches(i+1),:);
end

if length(prior_switches)==2
    tt = tt{1};
    priors = priors{1};
    slices = slices{1};
end

bump.out = bumpdur;
bump.in = backdur;
return