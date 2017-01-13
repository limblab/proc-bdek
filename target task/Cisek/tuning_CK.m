tune_ranges = {[500 1000],[0 350],[0 1000],[-100 200]};%,[0 200]};
tune_align = {6,7,8,20};%,20};
brain_area = 'M1';
MODthresh = 0;

tuning_dir_col = 17;

prediction_day_indices = 2;
CO_index = 1;


co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
% [PDS,reppd,reprch,dfrompd] = deal(cell(length(tune_align),1));
% for aligns = 1:length(tune_align)

[~,~,~,fulldat] = ...
         tuning_types_PD_func_COS(alldays,brain_area,CO_index,tuning_dir_col,tune_align,tune_ranges);

    PDS = get_PDS_fulldat(fulldat,pi/2);

%     for i = 1:size(PDS{aligns},2)
%         reppd{aligns}{i} = repmat(PDS{aligns}(:,i)',size(alldays(prediction_day_indices).tt,1),1);
%         reprch{aligns}{i} = repmat(alldays(prediction_day_indices).tt(:,tuning_dir_col),1,size(reppd{aligns}{i},2));%pref;
% 
%         dfrompd{aligns}{i} = circ_dist(reppd{aligns}{i},reprch{aligns}{i});
%     end
    
% end



