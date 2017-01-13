brain_area = 'PMd';
loop_alignments = {6};
loop_ranges = {[500 1000]};

baseline_window = 200;

CO_index = 1;
prediction_day = 2;
reachdir_col = 17;

dt = diff(loop_ranges{1})/1000;

co_units = eval(sprintf('alldays(1).%s_units',brain_area));
tt_pred = alldays(prediction_day).tt;
coangs = alldays(CO_index).tt(:,reachdir_col);
sessangs = alldays(prediction_day).tt(:,reachdir_col);
%% Get Baselines
baseraster = cell(length(co_units),1);
for i = 1:length(co_units) % Loop through units
    clc; fprintf('Calculating Baselines: (%d/%d)\n',i,length(co_units));

    baseraster{i} = raster_plot(co_units{i},alldays(CO_index).tt,[-baseline_window,0]./1000,6,0,'none');
end
baseline_levels = cell2mat(cellfun(@(x) nanmean(nansum(x{1},2)./(size(x{1},2)/1000),1),baseraster,'UniformOutput',0));  
tbt_baselines = repmat(baseline_levels',size(alldays(prediction_day).tt,1),1);
%% Get Tuning Curves
fullcurves = cell(length(co_units),1);
varcurves = zeros(629,length(co_units));
TC = zeros(629,length(co_units));
 for i = 1:length(co_units) % Loop through units
    clc; fprintf('%d/%d\n',i,length(co_units));

    % set time bin edges
    t1 = loop_ranges{1}(1);
    t2 = loop_ranges{1}(2);

    [cents,tune,~,~,~,~,fullcm] = co_tuning(co_units,alldays(CO_index).tt(:,reachdir_col),...
                             alldays(CO_index).tt,t1,t2,loop_alignments{1},i);
%     fullc = cell(1,size(fullcm,2));
%     for z = 1:size(fullcm,2)
%         fullc{z} = fullcm{:,z};
%     end
    fullc = fullcm;
                         
    fullcurves{i} = cellfun(@(x) sum(x,2),fullc,'UniformOutput',0);
    cents_rep = [cents-2*pi cents cents+2*pi]; tune_rep = [tune' tune' tune'];
    TC(:,i) = interp1(cents_rep,tune_rep,0:0.01:2*pi);
    binvars = cellfun(@(x) var(x),fullcurves{i});
    varcurves(:,i) = interp1(cents_rep,[binvars' binvars' binvars'],0:0.01:2*pi);
 end

%% Get Counts
unit_rates = zeros(size(tt_pred,1),length(co_units));
for q = 1:length(co_units)

    clc; fprintf('Computing Rasters...\n');

    out = raster_get(co_units{q},tt_pred,loop_ranges{1}./1000,loop_alignments{1});

    unit_rates(:,q) = nansum(out,2)./(diff(loop_ranges{1})./1000); % in FIRING RATE
end
%%
Araw = unit_rates*dt; % Spike Counts
TC_hz = TC;
TC_counts = TC*dt;
thets = 0:0.01:2*pi;

paren = exp(-TC_hz*dt)
for i = 1:size(Araw,1)

    TCnonzero = TC_hz; TCnonzero(TCnonzero==0)=NaN;
    pdensA = exp(nansum(repmat(Araw(i,:),629,1).*log(TCnonzero)-TCnonzero*dt,2));

    bootinds = bootstrp(1000,@(x) x,1:length(co_units));
  
    pdens = zeros(1000,629);
    for b = 1:size(bootinds,1)
        bootn = bootinds(b,:);
        
        part1 = 
        pdens(b,:) = exp(nansum(repmat(Araw(i,bootn),629,1).*log(TCnonzero(:,bootn))-TCnonzero(:,bootn)*dt,2));
    end
    
end

%%
for i = 1:size(Araw,1)
    
    for q = 1:length(co_units)
    
        nonscale = poisspdf(Araw(i,q),TC_counts(:,q));
        probq(:,q) = 1./sum(nonscale)*nonscale; 
    end
end

%%
[Bhat,Bhat_base,KLd] = deal(zeros(length(co_units),8));
[Bhatcurves,Bhatcurves_base,KLdcurves] = deal(zeros(629,length(co_units)));
for q = 1:length(co_units)
    for i = 1:8
        
        samedir = fullcurves{q}{i};
        otherdirs = fullcurves{q}; 
        otherdirs(i) = [];
        rnge = 0:max(cell2mat(fullcurves{q}));
        
       
        [y,x] = hist(cell2mat(otherdirs),rnge);
        
        Qx = y./sum(y);
        
        [y2,x2] = hist(samedir,rnge);
        
        Px = y2./sum(y2);
        
        [y3,x3] = hist(nansum(baseraster{q}{1},2),rnge);
        
        Bx = y3./sum(y3);
        
        Bhat(q,i) = -log(sum(sqrt(Px.*Qx)));
        Bhat_base(q,i) = -log(sum(sqrt(Px.*Bx)));
        
        kldstep1 = Px.*log(Px./Qx);
        KLd(q,i) = nansum(kldstep1(~isinf(kldstep1)));
    end
    
    Bhatcurves(:,q) = interp1(cents_rep,repmat(Bhat(q,:),1,3),0:0.01:2*pi);
    Bhatcurves_base(:,q) = interp1(cents_rep,repmat(Bhat_base(q,:),1,3),0:0.01:2*pi);
    KLdcurves(:,q) = interp1(cents_rep,repmat(KLd(q,:),1,3),0:0.01:2*pi);
    
end

%%
PDi = zeros(length(co_units),1);
for i = 1:length(co_units)
    PDi(i) = find(TC_counts(:,i)==max(TC_counts(:,i)),1,'first');
end
%%
[Dbhat,Dbhat_base,DKLd,DPD] = deal(zeros(size(tt_pred,1),length(co_units)));
dif_act = zeros(size(tt_pred,1),length(co_units));
for i = 1:size(tt_pred,1)
    rd = tt_pred(i,reachdir_col);
    ri = find(abs(circ_dist(rd,thets))==min(abs(circ_dist(rd,thets))));
    
%     expect = TC_counts(ri,:);
    expect = tbt_baselines(i,:);
    actual = Araw(i,:);
    
    dif_act(i,:) = actual-expect;
    Dbhat(i,:) = Bhatcurves(ri,:);
    Dbhat_base(i,:) = Bhatcurves_base(ri,:);
    DKLd(i,:) = KLdcurves(ri,:);
    Dvar(i,:) = varcurves(ri,:);
    Dmean(i,:) = TC_counts(ri,:);
    Dfano(i,:) = fanodir(ri,:);
%     Dovervar(i,:) = overVAR(ri,:);

    for q = 1:length(co_units)
        DPD(i,q) = abs(circ_dist(thets(ri),thets(PDi(q)))); 
    end
    
    
end

