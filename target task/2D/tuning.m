%TUNE_RANGES = {[0:25:800],[0:25:400]};
TUNE_RANGES = {[0:50:1000],[0:50:500],[0:50:1000],[0:50:400]};
%TUNE_ALIGNS = {'target','go'};
TUNE_ALIGNS = {6,7,8,9};
brain_area = 'PMd';
CO_index = 1;
prediction_day = 2;
MODthresh = 0;

%
co_units = eval(sprintf('alldays(1).%s_units',brain_area));
if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end

[TUNES, PDS,MODS] = deal(cell(length(co_units),1));
bestPD = NaN(length(co_units),3);

for i = 1:length(co_units)
    
    for loopthrough = 1:length(TUNE_ALIGNS)

        tune_ranges = TUNE_RANGES{loopthrough};
        tune_align = TUNE_ALIGNS{loopthrough};
        indcount = 1;
        for ranges = 1:length(tune_ranges)-1
        tune_range = tune_ranges(ranges:ranges+1);

             clc; fprintf('Neuron: %d/%d\n',...
                 i,length(co_units));

            % set time bin edges
            t1 = tune_range(1);
            t2 = tune_range(2);

            [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(CO_index).tt(:,17),...
                alldays(CO_index).tt,t1,t2,tune_align,i);

            % Smooth and fix edges on tuning curves
            front_cent = cents;
            back_cent = cents;

            front_tune = tune; front_low = tune_low; front_high = tune_high;
            back_tune = tune; back_low = tune_low; back_high = tune_high;

            tune_cents = [back_cent-2*pi cents front_cent+2*pi];
            tune_pad = [back_tune;tune;front_tune];
 
            interp_cents = tune_cents(1):.01:tune_cents(end);
            tune_interp = interp1(tune_cents,tune_pad,interp_cents);
           
            smooth_tune = smooth(tune_interp,100);
            
            wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
            wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 

            TUNES{i}{loopthrough}(indcount,:) = wrapped_tune;
            pds = wrapped_cents(wrapped_tune==max(wrapped_tune));
            if length(pds)==1
            	PDS{i}(indcount) = pds;
            else
                PDS{i}(indcount) = NaN;
            end
            MODS{i}(indcount) = max(wrapped_tune)-min(wrapped_tune);
         
            indcount = indcount +1; 
        end
    end
    
    mod_compare = find(MODS{i}==max(MODS{i}));
    if length(mod_compare) == 1
        bestPD(i,:) = [PDS{i}(mod_compare), mod_compare, MODS{i}(mod_compare)];
        %bestPD(i,:) = PDS{i}(4);
    else
        bestPD(i,:) = [PDS{i}(mod_compare(1)) NaN MODS{i}(mod_compare(1))];
        %bestPD(i,:) = NaN;
    end
    
end


reppd = repmat(bestPD(:,1)',size(alldays(prediction_day).tt,1),1);
reprch = repmat(alldays(prediction_day).tt(:,10),1,size(reppd,2));

dfrompd = circ_dist(reppd,reprch);