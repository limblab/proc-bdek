
for bin = 1:size(neurons,2)
    figure; hold on ; 
    for i = 1:22
        for j = 1:size(neurons{i,bin}.raster,2)
            plot(nanmean(neurons{i,bin}.raster(:,j)),nanvar(neurons{i,bin}.raster(:,j)),'b.')
        end
    end
    title(sprintf('Bin %d',bin));
end

%% Variance Estimate
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),7));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);
    
    errorbar(times,nanmean(popvar(:,sortinds(aa:bb)),2),1.96*nanstd(popvar(:,sortinds(aa:bb)),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);

%% End Prob
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),5));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);

    errorbar(times,nanmean(el(:,sortinds(aa:bb)),2),1.96*nanstd(el(:,sortinds(aa:bb)),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);


%% Target Prob
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),5));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);

    errorbar(times,nanmean(tl(:,sortinds(aa:bb)),2),1.96*nanstd(tl(:,sortinds(aa:bb)),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);


%% Endpoint Preference
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),4));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);
    
    endpt_pref = end_peak - targ_peak;

    errorbar(times,nanmean(endpt_pref(:,sortinds(aa:bb)),2),1.96*nanstd(endpt_pref(:,sortinds(aa:bb)),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);

%% Target Preference
times = (time_bins(2:end)+time_bins(1:end-1))./2;
%c = colormap;

devs = abs(adj_tt{2}(:,10)-adj_tt{2}(:,9));

[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),4));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);
    
    targ_pref = tl.*v_estimate;

    errorbar(times,nanmean(tl(:,sortinds(aa:bb)),2),1.96*nanstd(tl(:,sortinds(aa:bb)),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);

%% Endpoint Preference Over time
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortdevs(abs(sortdevs)>pi/2) = [];
sortinds(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),4));

figure; hold on; 
for i = 1:length(devedges)-1
    
    aa = devedges(i);
    bb = devedges(i+1);
    
    targ_pref = el.*v_estimate;

    errorbar(times,nanmean(targ_pref(:,aa:bb),2),1.96*nanstd(targ_pref(:,aa:bb),[],2)./sqrt(bb-aa),'Color',c(i*8,:));
    
    leg{i} = num2str(i);
end
legend(leg);

%% Var over Time
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortinds(abs(sortdevs)>pi/2) = [];
sortdevs(abs(sortdevs)>pi/2) = [];

%devedges = round(linspace(1,length(sortinds),4));
devedges = [1 165 591 756];

figure; hold on; 
for i = 1:length(devedges)-1
    
    if ismember(i,[1 3])
    aa = devedges(i);
    bb = devedges(i+1);
    
    errorbar(times,nanmean(popvar3(:,aa:bb),2),...
        1.96*nanstd(popvar3(:,aa:bb),[],2)./sqrt(bb-aa),'Color',c(i*5,:));
    
    leg{i} = num2str(i);
    end
end
legend(leg);
xlabel('Time from Target','FontSize',14);
ylabel('1/max(P(theta))','FontSize',14);


%%
binv = cell(length(conds),1);
numbins = 5;

figure; hold on; 
for i = 1:length(conds)
    
    [binv{i,1} binv{i,2}] = bin_array(v_estimate(i,:),1,numbins);
    errorbar(100/5:100/5:100,binv{i,1},1.96*binv{i,2}./sqrt(size(v_estimate,2)),'Color',c(i*10,:))
    
end
    

%% Endpoint Preference Over time
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

figure; hold on; 
    
targ_pref = tl2-el2;

errorbar(times,nanmean(targ_pref(:,highdevinds_early),2),...
    1.96*nanstd(targ_pref(:,highdevinds_early),[],2)./sqrt(length(highdevinds_early)),'Color',c(1,:));

errorbar(times,nanmean(targ_pref(:,highdevinds_late),2),...
    1.96*nanstd(targ_pref(:,highdevinds_late),[],2)./sqrt(length(highdevinds_late)),'Color',c(50,:));

legend('Early','Late');

%% Target preference over Time
times = (time_bins(2:end)+time_bins(1:end-1))./2;
c = colormap;

devs = adj_tt{2}(:,10)-adj_tt{2}(:,9);
[sortdevs,sortinds] = sortrows(devs);
sortinds(abs(sortdevs)>pi/2) = [];
sortdevs(abs(sortdevs)>pi/2) = [];

devedges = round(linspace(1,length(sortinds),4));
%devedges = [1 165 591 756];  
targ_pref = tl2 - el2;

figure; hold on; 
for i = 1:length(devedges)-1
    
    if ismember(i,[1 3])
    aa = devedges(i);
    bb = devedges(i+1);
    
    errorbar(times,nanmean(targ_pref(:,aa:bb),2),...
        1.96*nanstd(targ_pref(:,aa:bb),[],2)./sqrt(bb-aa),'Color',c(i*5,:));
    
    leg{i} = num2str(i);
    end
end
legend(leg);
xlabel('Time from Target','FontSize',14);
ylabel('1/max(P(theta))','FontSize',14);

%%
figure; hold on;
errorbar(times,nanmean(popconc,2),...
        1.96*nanstd(popconc,[],2)./sqrt(size(popconc,2)),'b');
%%
pc1 = popconc(:,1:300);
pcend = popconc(:,end-299:end);
figure; hold on;
errorbar(times,nanmean(pc1,2),...
    1.96*nanstd(pc1,[],2)./sqrt(size(pc1,2)),'b');

errorbar(times,nanmean(pcend,2),...
    1.96*nanstd(pcend,[],2)./sqrt(size(pcend,2)),'r');

%%
c = colormap;
figure; hold on; 
for i = 1:size(popconc,1)
    
    Y = popconc(i,:)';
    X = [ones(size(popconc,2),1) (1:size(popconc,2))'];
    
    [B,Bint] = regress(Y,X);
    
    slope_of_conc(i) = B(2);
    param_conf(:,i) = Bint(2,:)';
    

    
    %plot(mod(popmean(i,:),pi),'.','Color',c(i*8,:));
end

plot(slope_of_conc,'b')
plot(param_conf(1,:),'r'); plot(param_conf(2,:),'r')    

%%
figure; hold on;
for i = 1:length(cos_tunes)
    
    plot(wrapped_cents,cos_tunes{i}.tuning);
end


%%
figure; hold on;
for i = 1:length(cos_tunes)
pref = wrapped_cents(find(cos_tunes{i}.norm == max(cos_tunes{i}.norm),1,'first'));
plot([0 cos(pref)],[0 sin(pref)],'b-');

end