brain_area = 'PMd';

loop_alignments = {'target','go'};
loop_ranges = {[-200 0:50:1000],[-50:50:400]};
boot_type = 'trials';
FR_thresh = 0;

CO_index = 1;
if exist('pdi','var')
    prediction_day_indices = pdi;
else
    prediction_day_indices = 2;
end

[EorI,tune_widths,Av_rsp] = deal(cell(1,2));
baselines = zeros(length(eval(sprintf('alldays(1).%s_units',brain_area))),1);
for loopthrough = 1:2
    
tune_align = loop_alignments{loopthrough};
tune_ranges = loop_ranges{loopthrough};

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

            if  t1 <= 0 && t2 <= 0 && loopthrough == 1

            [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,[nan; alldays(CO_index).tt(1:end-1,10)],...
            alldays(CO_index).tt,t1,t2,tune_align,i,1);
            
            else
                
            [wrapped_cents,wrapped_tune,p,wrapped_low,wrapped_high] = co_tuning_VM(co_units,alldays(CO_index).tt(:,10),...
            alldays(CO_index).tt,t1,t2,tune_align,i,1);

            end
            
            tune_widths{loopthrough}(i,ranges) = p(3);

            if t1 <=0 && t2<=0 && loopthrough == 1
                baselines(i) = mean(wrapped_tune);
            else
                EorI{loopthrough}(i,ranges) = 2*(mean(wrapped_tune) > baselines(i))-1;
            end
            
            Av_rsp{loopthrough}(i,ranges) = mean(wrapped_tune);
            
        end
    end
end
