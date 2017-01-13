% for i = 1:length(alldays)
%     alldays(i).tt(floor(10*(alldays(i).tt(:,3)-floor(alldays(i).tt(:,3))))==9,:) = [];
%     alldays(i).tt(isnan(alldays(i).tt(:,3)),:) = [];
% end
tune_ranges = [600 800];
tune_align = 'target';
brain_area = 'PMd_narrow';
MODthresh = 0;
if exist('pdi','var')
    prediction_day_indices = pdi;
else
    prediction_day_indices = 2;
end
CO_index = 1;
group_col = 3; 
reachdir_col = 10;

[neurons, PD,PD1,nears,fars,nearneurs,farneurs,DFPD,VM_P] = deal(cell(length(tune_ranges)-1,1));
for ranges = 1:length(tune_ranges)-1
    tune_range = tune_ranges(ranges:ranges+1);
    % Initialize
    co_units = eval(sprintf('alldays(1).%s_units',brain_area));
    if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
    %neurons = cell(length(co_units),length(tune_range)-1);

    for i = 1:length(co_units) % Loop through units
         clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));

        % set time bin edges
        t1 = tune_range(1);
        t2 = tune_range(2);

        [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,reachdir_col),...
            alldays(CO_index).tt,t1,t2,tune_align,i);
       
        %%Smooth and fix edges on tuning curves
        front_cent = cents;
        back_cent = cents;

        front_tune = tune; front_low = tune_low; front_high = tune_high;
        back_tune = tune; back_low = tune_low; back_high = tune_high;

        tune_cents = [back_cent-2*pi cents front_cent+2*pi];
        tune_pad = [back_tune;tune;front_tune];
        low_pad = [back_low;tune_low;front_low]; 
        high_pad = [back_high;tune_high;front_high];

        interp_cents = tune_cents(1):.01:tune_cents(end);
        tune_interp = interp1(tune_cents,tune_pad,interp_cents);
        low_interp = interp1(tune_cents,low_pad,interp_cents);
        high_interp = interp1(tune_cents,high_pad,interp_cents);

        smooth_tune = smooth(tune_interp,100);
        smooth_low = smooth(low_interp,100);
        smooth_high= smooth(high_interp,100);

        wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
        wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
        wrapped_low = smooth_low(interp_cents >= 0 & interp_cents <2*pi); 
        wrapped_high = smooth_high(interp_cents >= 0 & interp_cents <2*pi);

        neurons{ranges}{i}.tuning = wrapped_tune;
        neurons{ranges}{i}.tuning_low = wrapped_low;
        neurons{ranges}{i}.tuning_high = wrapped_high;
      
        PD{ranges}{i} = wrapped_cents(wrapped_tune==max(wrapped_tune));
        if length(PD{ranges}{i}) < 2; PD1{ranges}(i,:) = PD{ranges}{i}; else
            PD1{ranges}(i,:) = NaN; end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,reachdir_col),...
%             alldays(CO_index).tt,t1,t2,tune_align,i,1);
%         
%         VM_P{ranges}{i} = p;
%         
%         PD{ranges}{i} = wrapped_cents(wrapped_tune==max(wrapped_tune));
%         PD{ranges}{i} = wrapped_cents(abs(circ_dist(wrapped_cents,p(end)))==min(abs(circ_dist(wrapped_cents,p(end)))));
%         PD1{ranges}(i,:) = PD{ranges}{i};
%         if length(PD{ranges}{i})>1
%             PD1{ranges}(i,:) = NaN;%PD{ranges}{1};
%         else
%             PD1{ranges}(i,:) = PD{ranges}{i};
%         end
%         
%         neurons{ranges}{i}.tuning = wrapped_tune;
%         neurons{ranges}{i}.tuning_low = wrapped_low;
%         neurons{ranges}{i}.tuning_high = wrapped_high;
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    end
%     nears{ranges} = find(abs(circ_dist(PD1{ranges},circ_mean(alldays(prediction_day_indices).tt(:,reachdir_col))))< 2*pi/4);
%     fars{ranges} = find(abs(circ_dist(PD1{ranges},circ_mean(alldays(prediction_day_indices).tt(:,reachdir_col))))> 2*pi/4);
%     
    reppd = repmat(PD1{ranges}',size(alldays(prediction_day_indices).tt,1),1);
    %reprch = repmat(alldays(prediction_day_indices).tt(:,reachdir_col),1,size(reppd,2));
    reprch = pi/2*ones(size(reppd));
    
    dfrompd = circ_dist(reppd,reprch);
    
%     nearinds = abs(dfrompd) < pi/2; 
%     farinds  = abs(dfrompd) > pi/2; 
%      
%     nearneurs{ranges} = nan(size(nearinds)); nearneurs{ranges}(nearinds)=1;
%     farneurs{ranges} = nan(size(farinds)); farneurs{ranges}(farinds)=1;
    
    DFPD{ranges} = dfrompd;
    
end
centers = round(wrapped_cents*1000)./1000;
