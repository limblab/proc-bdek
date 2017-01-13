brain_area = 'PMd';

loop_alignments = {'target','go'};
loop_ranges = {[600 800],[50 250]};
FR_thresh = 0;

CO_index = 1;
reachdir_col = 10;

[percdiff, xsfp, diffb, tds, totd,neurons,FD,PDneurs,ODneurs,ORTHneurs,PD,PD1,TA,...
    Firings,dfrompd,mbound,excited,suppressed] = deal(cell(1,length(loop_alignments)));
% Initialize
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
for loopthrough = 1:length(loop_alignments)
    
tune_align = loop_alignments{loopthrough};
tune_ranges = loop_ranges{loopthrough};

    for ranges = 1:length(tune_ranges)-1
        tune_range = tune_ranges(ranges:ranges+1);
        for i = 1:length(co_units) % Loop through units
             clc; fprintf('%d/%d Tuning: (%d/%d)\n',loopthrough,length(loop_alignments),i,length(co_units));
% 
%             % set time bin edges
%             t1 = tune_range(1);
%             t2 = tune_range(2);
% 
%             if  t1 <= 0 && t2 <= 0 && strcmp(tune_align,'target')
%                 [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,[nan; alldays(CO_index).tt(1:end-1,reachdir_col)],...
%                 alldays(CO_index).tt,t1,t2,tune_align,i,0);
%             else
%                 [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),...
%                 alldays(CO_index).tt,t1,t2,tune_align,i,0);
%             end
%             if length(wrapped_low)~=length(wrapped_tune); wrapped_low=wrapped_tune; wrapped_high=wrapped_tune; end
%             
%             %PD = wrapped_cents(wrapped_tune==max(wrapped_tune));
%             PD = atan2(p(3),p(2));
%             
%             if length(PD)>1 || isempty(PD)
%                 PD1{loopthrough}{ranges}(i,:) = NaN;
%             else
%                 PD1{loopthrough}{ranges}(i,:) = PD;
%             end
%             
            
            [PD,PD_bounds] = PD_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),alldays(CO_index).tt,t1,t2,tune_align,i,1000,[2.5 97.5]);


            neurons{loopthrough}{ranges}{i}.tuning = wrapped_tune';
        end        
    end
end

allneurs = horzcat(neurons{:});
[tunearrays] = deal(cell(1,length(allneurs)));
for i = 1:length(allneurs)
    tunearrays{i} = cell2mat(cellfun(@(x) x.tuning,allneurs{i},'UniformOutput',0))';
end

tuned = cell2mat(cellfun(@(x) ~isnan(sum(x,2)),tunearrays,'UniformOutput',0));

tune_counts = [sum(tuned) sum(sum(tuned,2)==0)];

tuning_depths = cell2mat(cellfun(@(x) max(x,[],2)-min(x,[],2),tunearrays,'UniformOutput',0));
tuning_prefs = double(tuning_depths == repmat(max(tuning_depths,[],2),1,length(allneurs)));
tuning_prefs(tuning_prefs==0) = NaN;

pref_counts = nansum(tuning_prefs);

%%
PDS_ranges = cell2mat(cellfun(@(x) x{1},PD1,'UniformOutput',0));
T = PDS_ranges.*tuning_prefs;
no_tuning = sum(isnan(T),2)==size(T,2);

best_PDS = nansum(T,2);
best_PDS(no_tuning) = NaN;


%% Get Counts
% tt_co = alldays(CO_index).tt;
% tbin = 50;
% 
% maxtriallength = max(tt_co(:,7)-tt_co(:,5));
% triallength = round(maxtriallength.*1000/tbin).*tbin/1000;
% allrast = zeros(length(co_units),size(tt_co,1).*triallength*1000/tbin);
% for q = 1:length(co_units)
%     clc; fprintf('Unit: %d/%d\n',q,length(co_units)); 
%     raster = raster_plot(co_units{q},tt_co,[0 triallength],'target',0,'none');
%     binrast = bin_array(raster{1},size(raster{1},1),length(raster{1})/tbin,'sum')./(tbin/1000);
%     
%     allrast(q,:) = reshape(binrast,[],1);
% end
% op_ranges = prctile(allrast,99,2);
% op_ranges = var(allrast,[],2);

