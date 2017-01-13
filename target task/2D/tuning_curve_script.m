neurs = cell(3,1);
for j = 1:3 %DAYS
    for i = 1:length(co_units) % Loop through units
        
        cur_tt = alldays(j).tt;
        linds = find(cur_tt(:,3)==max(cur_tt(:,3)));
        hinds = find(cur_tt(:,3)==min(cur_tt(:,3)));
        
        indices = {linds,hinds};
        
        for q = 1:2
            clc; fprintf('Tuning: (%d/%d)\n',i,length(co_units));

            % set time bin edges
            t1 = tune_range(1);
            t2 = tune_range(2);

            [cents,tune,tune_low,tune_high] = co_tuning(alldays(j).PMd_units,alldays(j).tt(indices{q},10),alldays(j).tt(indices{q},:),t1,t2,tune_align,i);

            % Smooth and fix edges on tuning curves
            front_cent = cents;
            back_cent = cents;

            front_tune = tune;
            front_low = tune_low;
            front_high = tune_high;

            back_tune = tune;
            back_low = tune_low;
            back_high = tune_high;

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

            neurs{j}{q}{i}.tuning = wrapped_tune;
            neurs{j}{q}{i}.tuning_low = wrapped_low;
            neurs{j}{q}{i}.tuning_high = wrapped_high; 
            %    {session portion}{trial portion}{unit}
            %neurons{i}.raster = tune_rast;

        end
    end
end
%%
dif1 = cell(2,1);
dif2 = cell(2,1);
for i =1:length(neurs{1}{1})
    for q = 1:2
        base = nanmean(neurs{1}{q}{i}.tuning);
        modul = max(neurs{1}{q}{i}.tuning) - min(neurs{1}{q}{i}.tuning);
        d1 = neurs{2}{q}{i}.tuning - neurs{1}{q}{i}.tuning;
        d2 = neurs{3}{q}{i}.tuning - neurs{1}{q}{i}.tuning;

        dif1{q}(i) = nanmean(d1(~isinf(d1)))./modul; 
        dif1{q}(isnan(dif1{q}))=1; dif1{q}(isinf(dif1{q}))=1;

        dif2{q}(i) = nanmean(d2(~isinf(d2)))./modul;
        dif2{q}(isnan(dif2{q}))=1; dif2{q}(isinf(dif2{q}))=1;
    end
end
   
[ld1.l, ud1.l] = boot_bounds(5000,@mean,dif1{1},2.5,97.5);   
[ld1.h, ud1.h] = boot_bounds(5000,@mean,dif1{2},2.5,97.5);   

[ld2.l, ud2.l] = boot_bounds(5000,@mean,dif2{1},2.5,97.5);   
[ld2.h, ud2.h] = boot_bounds(5000,@mean,dif2{2},2.5,97.5);   

figure; hold on;
errorbar([1,2],[mean(dif1{1}),mean(dif1{2})],[ld1.l,ld1.h],[ud1.l,ud1.h]);
title('Diff 1');
figure; hold on;
errorbar([1,2],[mean(dif2{1}),mean(dif2{2})],[ld2.l,ld2.h],[ud2.l,ud2.h]);
title('Diff 2');
%%
figure; hold on;
col_plots = {'b','g','r'};
for i = 1:length(neurs{1})
    for j = 1:3
        plot(wrapped_cents,neurs{j}{i}.tuning,col_plots{j});
    end
    pause; cla
end

