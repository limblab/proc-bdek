brain_area = 'first';

loop_alignmentsPD = {'target','target','go'};
% loop_alignmentsPD = {'target','go'};
loop_rangesPD = {[50 250],[600 800],[50 250]};
% loop_rangesPD = {[600 800],[50 250]};

FR_thresh = 0;

CO_index = 1;
reachdir_col = 10;

[PDS,PDS_bounds,PDS_dist,MAGS_dist,curve_b] = deal(cell(1,length(loop_alignmentsPD)));
% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
for loopthrough = 1:length(loop_alignmentsPD)
    
    tune_align = loop_alignmentsPD{loopthrough};
    tune_range = loop_rangesPD{loopthrough};
    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',loopthrough,length(loop_alignmentsPD),i,length(co_units));
        [PD,PD_bounds,PD_dist,MAG_dist,curve_params] = PD_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),...
                                      alldays(CO_index).tt,tune_range(1),tune_range(2),tune_align,i,[],1000,[2.5 97.5]);
        PDS{loopthrough}(i,:) = PD;
        PDS_bounds{loopthrough}(i,:) = PD_bounds;
        PDS_dist{loopthrough}(i,:) = PD_dist;
        MAGS_dist{loopthrough}(i,:) = MAG_dist;
        curve_b{loopthrough}(i,:) = curve_params;
    end        
end
%%
PDS_mat = cell2mat(PDS);
% PDS_intervals = cell2mat(cellfun(@(x) (x(:,2)-x(:,1))./pi.*180,PDS_bounds,'UniformOutput',0));
% PDS_intervals = cellfun(@(x) prctile(x,[2.5,97.5],2)./pi.*180,PDS_dist,'UniformOutput',0);
interval_widths = cell2mat(cellfun(@(x) circ_dist(prctile(x,97.5,2),prctile(x,2.5,2))./pi.*180,PDS_dist,'UniformOutput',0));
intervals_ok = interval_widths>0&interval_widths<90;

% create matrix of intervals/magnitudes that pass test
MAGS = cell2mat(cellfun(@(x) mean(x,2),MAGS_dist,'UniformOutput',0));
MAGS_ok = cell2mat(cellfun(@(x) ttest(x')',MAGS_dist,'UniformOutput',0));
neurs_ok = double(intervals_ok&MAGS_ok);
neurs_ok(neurs_ok==0) = NaN;

smallest_intervals = (interval_widths.*neurs_ok)==repmat(min((interval_widths.*neurs_ok),[],2),1,size(interval_widths,2));
biggest_mags = (MAGS.*neurs_ok)==repmat(max((MAGS.*neurs_ok),[],2),1,size(MAGS,2));

passed_test = double(smallest_intervals);
passed_test(~passed_test) = nan;

best_PDS = nansum(passed_test.*PDS_mat,2);
best_PDS(sum(isnan(passed_test),2)==length(loop_rangesPD)) = NaN;
