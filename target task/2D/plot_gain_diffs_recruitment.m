figure; hold on; 
for q = 1:2
    subplot(1,2,1); hold on;
%     plot(averaged_fmL{q},averaged_pmL{q},'Color',plotcolors{q},'LineWidth',3);
%    
%     patch([averaged_fmL{q}' fliplr(averaged_fmL{q}')],[boot_pmL_low{q}' fliplr(boot_pmL_high{q}')],...
%         plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0);
%         
%     plot(averaged_fmL{q},averaged_epmL{q},'k--','LineWidth',2);

    %%% PLOT GAIN CHANGE %%%%
    plot(averaged_fmL{q},averaged_pmL{q}./averaged_epmL{q},'Color',plotcolors{q},'LineWidth',3);
    patch([averaged_fmL{q}' fliplr(averaged_fmL{q}')],...
        [boot_pmL_low{q}'./averaged_epmL{q}' fliplr(boot_pmL_high{q}'./averaged_epmL{q}')],...
        plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0);


    ylim([0 3]);
    title('Low Likelihood Unc');

    
    subplot(1,2,2); hold on;
%     plot(averaged_fmH{q},averaged_pmH{q},'Color',plotcolors{q},'LineWidth',3);
% 
%     patch([averaged_fmH{q}' fliplr(averaged_fmH{q}')],[boot_pmH_low{q}' fliplr(boot_pmH_high{q}')],...
%         plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0)
% 
%     plot(averaged_fmH{q},averaged_epmH{q},'k--','LineWidth',2);

    %%% PLOT GAIN CHANGE %%%%
    plot(averaged_fmH{q},averaged_pmH{q}./averaged_epmH{q},'Color',plotcolors{q},'LineWidth',3);
    patch([averaged_fmH{q}' fliplr(averaged_fmH{q}')],...
        [boot_pmH_low{q}'./averaged_epmH{q}' fliplr(boot_pmH_high{q}'./averaged_epmH{q}')],...
        plotcolors{q},'FaceAlpha',.5,'EdgeAlpha',0);


    ylim([0 3]);
    title('High Likelihood Unc');
end