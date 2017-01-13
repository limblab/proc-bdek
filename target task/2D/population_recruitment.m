%%
session_nums = [2];
tbin = 2;
symmetric_plot = 0;

[binndfmincL,binndfmincH,binndpmincL,binndpmincH,binndepmincL,binndepmincH,...
 averaged_fmL,averaged_fmH,averaged_pmL,averaged_pmH,averaged_epmL,averaged_epmH,...
 boot_pmL_low,boot_pmL_high,boot_epmL_low,boot_epmL_high,...
 boot_pmH_low,boot_pmH_high,boot_epmH_low,boot_epmH_high]...
 = deal(cell(length(session_nums),1));

plotcolors = {'r','b'};
figure; 
for q = 1:length(session_nums)
   
    fprintf('%d/%d\n',q,length(session_nums));
    
    prisess = session_nums(q);
    
    cur_tt = alldays(prisess).tt;
    dub_tt = [cur_tt; cur_tt];
    lind = find(cur_tt(:,3)==max(cur_tt(:,3)));
    hind = find(cur_tt(:,3)==min(cur_tt(:,3)));
    %lctt = length(cur_tt);
    %lind = lind(round((2/3)*lctt):round((3/3)*lctt));
    
    ecounts = expect_counts{prisess-1}{tbin};
    acounts = trial_counts{prisess-1}{tbin};

    minfr = zeros(length(neurons{tbin}),1);
    maxfr = zeros(length(neurons{tbin}),1);

    pds = cell(length(neurons{tbin}),1);
    for i = 1:length(neurons{tbin})

        minfr(i) = min(neurons{tbin}{i}.tuning);
        maxfr(i) = max(neurons{tbin}{i}.tuning);

       [heights,pds{i}] = findpeaks(neurons{tbin}{i}.tuning);

    end

    numpeaks = cellfun(@length,pds);
    well_tuned = find(numpeaks==1); % single peaked tuning curve
    prefd = cellfun(@mean,pds);

%
    move_dirs = alldays(prisess).tt(:,10);

    maxtc = maxfr(well_tuned);
    mintc = minfr(well_tuned);

    badn = find(maxfr - minfr < 0);
    vwell_tuned = well_tuned(~ismember(well_tuned,badn));

    perc_max = zeros(length(vwell_tuned),length(move_dirs));
    from_move = zeros(length(vwell_tuned),length(move_dirs));
    eperc_max = zeros(length(vwell_tuned),length(move_dirs));
    for i = 1:length(move_dirs)

        from_move(:,i) = wrapped_cents(prefd(vwell_tuned)) - move_dirs(i);

        from_move(from_move > pi) = from_move(from_move > pi) - 2*pi;
        
        %from_move = abs(from_move);

        perc_max(:,i) = (acounts(i,vwell_tuned)' - minfr(vwell_tuned))./(maxfr(vwell_tuned)-minfr(vwell_tuned));
        eperc_max(:,i) = (ecounts(i,vwell_tuned)' - minfr(vwell_tuned))./(maxfr(vwell_tuned)-minfr(vwell_tuned));

    end
    if symmetric_plot == 1
        lind = find(dub_tt(:,3)==max(dub_tt(:,3)));
        hind = find(dub_tt(:,3)==min(dub_tt(:,3)));
    end
    
    from_move = [from_move -from_move];
    perc_max = [perc_max perc_max];
    eperc_max = [eperc_max eperc_max];

    allfm = reshape(from_move,numel(from_move),1);
    allpm = reshape(perc_max,numel(perc_max),1);
    allepm = reshape(eperc_max,numel(eperc_max),1);

    allfmL = reshape(from_move(:,lind),numel(from_move(:,lind)),1);
    allfmH = reshape(from_move(:,hind),numel(from_move(:,hind)),1);
    allpmL = reshape(perc_max(:,lind),numel(perc_max(:,lind)),1);
    allpmH = reshape(perc_max(:,hind),numel(perc_max(:,hind)),1);
    allepmL = reshape(eperc_max(:,lind),numel(eperc_max(:,lind)),1);
    allepmH = reshape(eperc_max(:,hind),numel(eperc_max(:,hind)),1);

    [fm_increasingL,orderfmL] = sortrows(allfmL);
    pm_increasingL = allpmL(orderfmL);
    epm_increasingL = allepmL(orderfmL);

    [fm_increasingH,orderfmH] = sortrows(allfmH);
    pm_increasingH = allpmH(orderfmH);
    epm_increasingH = allepmH(orderfmH);

    binndfmincL{q} = bin_array(fm_increasingL,5,1);
    binndfmincH{q} = bin_array(fm_increasingH,5,1);

    binndpmincL{q} = bin_array(pm_increasingL,5,1);
    binndpmincH{q} = bin_array(pm_increasingH,5,1);

    binndepmincL{q} = bin_array(epm_increasingL,5,1);
    binndepmincH{q} = bin_array(epm_increasingH,5,1);
    
    binwindow = 10000;
    mask = ones(1,binwindow)/binwindow;
    
    averaged_fmL{q} = conv(fm_increasingL,mask,'valid');
    averaged_fmH{q} = conv(fm_increasingH,mask,'valid');
    
    averaged_pmL{q} = conv(pm_increasingL,mask,'valid');
    boot_wrap = @(x) conv(pm_increasingL,hist(x,1:binwindow)/binwindow,'valid');
    [boot_pmL_low{q}, boot_pmL_high{q}] = boot_bounds(100,boot_wrap,1:binwindow,2.5,97.5);
    
    averaged_pmH{q} = conv(pm_increasingH,mask,'valid');
    boot_wrap = @(x) conv(pm_increasingH,hist(x,1:binwindow)/binwindow,'valid');
    [boot_pmH_low{q}, boot_pmH_high{q}] = boot_bounds(100,boot_wrap,1:binwindow,2.5,97.5);
    
    averaged_epmL{q} = conv(epm_increasingL,mask,'valid');
%     boot_wrap = @(x) conv(epm_increasingL,hist(x,1:binwindow)/binwindow,'valid');
%     [boot_epmL_low{q}, boot_epmL_high{q}] = boot_bounds(1000,boot_wrap,1:binwindow,2.5,97.5);
    
    averaged_epmH{q} = conv(epm_increasingH,mask,'valid');
%     boot_wrap = @(x) conv(epm_increasingH,hist(x,1:binwindow)/binwindow,'valid');
%     [boot_epmH_low{q}, boot_epmH_high{q}] = boot_bounds(1000,boot_wrap,1:binwindow,2.5,97.5);
    
%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,1); hold on;
    plot(averaged_fmL{q},averaged_pmL{q},'Color',plotcolors{q},'LineWidth',3);
%     plot(-averaged_fmL{q},averaged_pmL{q},'--','Color',plotcolors{q},'LineWidth',3);
    
    patch([averaged_fmL{q}' fliplr(averaged_fmL{q}')],[boot_pmL_low{q}' fliplr(boot_pmL_high{q}')],...
        plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0);
    
%     patch([-averaged_fmL{q}' fliplr(-averaged_fmL{q}')],[boot_pmL_low{q}' fliplr(boot_pmL_high{q}')],...
%         plotcolors{q},'FaceAlpha',.25,'EdgeAlpha',0);
    
    plot(averaged_fmL{q},averaged_epmL{q},'k--','LineWidth',2);
%     plot(-averaged_fmL{q},averaged_epmL{q},'k--','LineWidth',2);
%     interpL = interp1([fliplr(-averaged_fmL{q}'),averaged_fmL{q}'],...
%     [fliplr(averaged_epmL{q}') averaged_epmL{q}'],-min(averaged_fmL{q}):.01:min(averaged_fmL{q}),'spline');
%     plot(-min(averaged_fmL{q}):.01:min(averaged_fmL{q}),interpL,'k--','LineWidth',2);
    
    ylim([0 1.2]);
    title('Low Likelihood Unc');

    
    subplot(1,2,2); hold on;
    plot(averaged_fmH{q},averaged_pmH{q},'Color',plotcolors{q},'LineWidth',3);
%     plot(-averaged_fmH{q},averaged_pmH{q},'--','Color',plotcolors{q},'LineWidth',3);
    
    patch([averaged_fmH{q}' fliplr(averaged_fmH{q}')],[boot_pmH_low{q}' fliplr(boot_pmH_high{q}')],...
        plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0)
    
%     patch([-averaged_fmH{q}' fliplr(-averaged_fmH{q}')],[boot_pmH_low{q}' fliplr(boot_pmH_high{q}')],...
%         plotcolors{q},'FaceAlpha',.25,'EdgeAlpha',0)
    
    plot(averaged_fmH{q},averaged_epmH{q},'k--','LineWidth',2);
%     plot(-averaged_fmH{q},averaged_epmH{q},'k--','LineWidth',2);
%     interpH = interp1([fliplr(-averaged_fmH{q}'),averaged_fmH{q}'],...
%     [fliplr(averaged_epmH{q}') averaged_epmH{q}'],-min(averaged_fmH{q}):.01:min(averaged_fmH{q}),'spline');
%     plot(-min(averaged_fmH{q}):.01:min(averaged_fmH{q}),interpH,'k--','LineWidth',2);

    ylim([0 1.2]);
    title('High Likelihood Unc');
end


