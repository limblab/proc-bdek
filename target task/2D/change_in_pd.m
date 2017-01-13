%%
% session_nums = [1 3];
% tune_ranges = [400 800];
% 
% tune_align = 'target';
% brain_area = 'PMd';
% 
% plotcolors = {'r','g','b'};
% tbin = 1;
% 
% neurs = cell(length(tune_ranges)-1,length(session_nums));
% for priorsess = 1:length(session_nums)
%     
%     ps = session_nums(priorsess);
%     
%     clc; fprintf('trial block %d/%d\n',priorsess,length(session_nums));
%     
%     for ranges = 1:length(tune_ranges)-1
%         tune_range = tune_ranges(ranges:ranges+1);
%         % Initialize
%         co_units = eval(sprintf('alldays(1).%s_units',brain_area));
%         if strcmp(brain_area,'M1'); co_units(end-1:end) =[]; end
% 
%         for i = 1:length(co_units) % Loop through units
%             clc; fprintf('%d/%d Tuning: (%d/%d)\n',ranges,length(tune_ranges)-1,i,length(co_units));
% 
%             % set time bin edges
%             t1 = tune_range(1);
%             t2 = tune_range(2);
% 
%             [cents,tune,tune_low,tune_high] = co_tuning(co_units,alldays(ps).tt(:,10),alldays(ps).tt,t1,t2,tune_align,i);
% 
%             % Smooth and fix edges on tuning curves
%             front_cent = cents;
%             back_cent = cents;
% 
%             front_tune = tune; front_low = tune_low; front_high = tune_high;
%             back_tune = tune; back_low = tune_low; back_high = tune_high;
% 
%             tune_cents = [back_cent-2*pi cents front_cent+2*pi];
%             tune_pad = [back_tune;tune;front_tune];
%             low_pad = [back_low;tune_low;front_low]; 
%             high_pad = [back_high;tune_high;front_high];
% 
%             interp_cents = tune_cents(1):.01:tune_cents(end);
%             tune_interp = interp1(tune_cents,tune_pad,interp_cents);
%             low_interp = interp1(tune_cents,low_pad,interp_cents);
%             high_interp = interp1(tune_cents,high_pad,interp_cents);
% 
%             smooth_tune = smooth(tune_interp,100);
%             smooth_low = smooth(low_interp,100);
%             smooth_high= smooth(high_interp,100);
% 
%             wrapped_cents = interp_cents(interp_cents >=0 & interp_cents< 2*pi);
%             wrapped_tune = smooth_tune(interp_cents >= 0 & interp_cents <2*pi); 
%             wrapped_low = smooth_low(interp_cents >= 0 & interp_cents <2*pi); 
%             wrapped_high = smooth_high(interp_cents >= 0 & interp_cents <2*pi);
% 
%             neurs{ranges,priorsess}{i}.tuning = wrapped_tune;
%             neurs{ranges,priorsess}{i}.tuning_low = wrapped_low;
%             neurs{ranges,priorsess}{i}.tuning_high = wrapped_high;
%             %neurs{i}.raster = tune_rast;
% 
%         end
%     end
% end

all_pds = cell(length(session_nums),1);

for q = 1:length(session_nums) % loop through prior sessions
   
    fprintf('Prior %d/%d\n',q,length(session_nums));
    
    % current prior session
    prisess = session_nums(q);
    
    cur_tt = alldays(prisess).tt; % trial table

    minfr = zeros(length(neurs{tbin,q}),1);
    maxfr = zeros(length(neurs{tbin,q}),1);

    %%% Get preferred directions from tuning descriptions %%%
    pds = cell(length(neurs{tbin,q}),1);
    for i = 1:length(neurs{tbin,q})

        minfr(i) = min(neurs{tbin,q}{i}.tuning);
        maxfr(i) = max(neurs{tbin,q}{i}.tuning);

       [heights,pds{i}] = findpeaks(neurs{tbin,q}{i}.tuning);

    end

    %%% Only use single-peaked neurs %%%
    numpeaks = cellfun(@length,pds);
    well_tuned = find(numpeaks==1); % single peaked tuning curve
    prefd = cellfun(@mean,pds);

    %%% Extract movement directions
    move_dirs = alldays(prisess).tt(:,10);

    %%% Of single-peaked neurs, eliminate those with too low FR %%%
    maxtc = maxfr(well_tuned);
    mintc = minfr(well_tuned);
    badn = find(maxfr - minfr <= 0);
    vwell_tuned = well_tuned(~ismember(well_tuned,badn));

    % pds
    all_pds{q} = wrapped_cents(prefd(vwell_tuned));

end

%%
figure; hold on;

xs = 0:0.01:2*pi;
plot(cos(xs),sin(xs),'k');

for q = 1:length(session_nums)
    plot(cos(all_pds{q}),sin(all_pds{q}),'o','Color',plotcolors{q});
end


%%
for i = 1:length(neurs{1})
    
    mod1(i) = max(neurs{1}{i}.tuning) - min(neurs{1}{i}.tuning);
    mod2(i) = max(neurs{2}{i}.tuning) - min(neurs{2}{i}.tuning);
    
    tune_change(i,:) = (neurs{2}{i}.tuning - neurs{1}{i}.tuning)/(mean([mod1(i) mod2(i)]));
    
    [~,pds1{i}] = findpeaks(neurons{1}{i}.tuning);
    [~,pds2{i}] = findpeaks(neurons{2}{i}.tuning);
end

figure; hold on;
plot(wrapped_cents,nanmean(tune_change,1),'b');

double_peaks = find(cellfun(@length,pds2)>1);

figure; hold on;
for i = 1:length(double_peaks)
    
    plot(wrapped_cents(pds2{double_peaks(i)}),i*ones(1,length(pds2{double_peaks(i)})),'ro');
    plot(wrapped_cents(pds1{double_peaks(i)}),i*ones(1,length(pds1{double_peaks(i)})),'b.');
end

good_1 = find(mod1 > 1);
%plot(wrapped_cents,mean(tune_change(good_1,:),1),'r.');