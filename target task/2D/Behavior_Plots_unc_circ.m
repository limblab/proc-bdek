%% Plots
% Nature Plot

Linds = find(plot_tt(:,3) == max(plot_tt(:,3)));
Hinds = find(plot_tt(:,3) == min(plot_tt(:,3)));

dump = 0;

linds = Linds(Linds>dump);
hinds = Hinds(Hinds>dump);

Nature_low = plot_tt(linds,9:10);
Nature_high = plot_tt(hinds,9:10);

Nature_all = [Nature_low;Nature_high];

[BL,NatureLsort] = sortrows(Nature_low(:,1));
[BH,NatureHsort] = sortrows(Nature_high(:,1));

figure; plot(plot_tt(linds,9),plot_tt(linds,10),'b.');
hold on; title(sprintf('%s:  %s/%s/%s',monkey,MO,DA,YE),'FontSize',18);
xlabel('Centroid Location (rad)','FontSize',20); ylabel('End Cursor Location (rad)','FontSize',20);
plot(plot_tt(hinds,9),plot_tt(hinds,10),'r.');

[slope_low,SL] = polyfit(Nature_low(:,1),Nature_low(:,2),1);
[slope_high,SH] = polyfit(Nature_high(:,1),Nature_high(:,2),1);

[YL,DELTAL] = polyconf(slope_low,Nature_low(:,1),SL,'predopt','curve');
[YH,DELTAH] = polyconf(slope_high,Nature_high(:,1),SH,'predopt','curve');

%plot([min(Nature_all(:,1)) max(Nature_all(:,1))],...
    %[slope_low(1)*min(Nature_all(:,1))+ slope_low(2) slope_low(1)*max(Nature_all(:,1)) + slope_low(2)],'b','LineWidth',3);

plot(Nature_low(:,1),YL,'b','LineWidth',2);
plot(Nature_high(:,1),YH,'r','LineWidth',2);

legend(sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_low(1),min(plot_tt(:,3))),...
sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_high(1),max(plot_tt(:,3))));

plot(repmat(BL,1,2),[YL(NatureLsort)-DELTAL(NatureLsort) ...
                      YL(NatureLsort)+DELTAL(NatureLsort)],'b--','LineWidth',1);
plot(repmat(BH,1,2),[YH(NatureHsort)-DELTAH(NatureHsort) ...
                      YH(NatureHsort)+DELTAH(NatureHsort)],'r--','LineWidth',1);

%Plot Errorbars
%     plot([min(Nature_all(:,1)) max(Nature_all(:,1))],...
%         [slope_low(1)*min(Nature_all(:,1))+ slope_low(2) slope_low(1)*max(Nature_all(:,1)) + slope_low(2)],'b','LineWidth',3);
% 
%     plot([min(Nature_all(:,1)) max(Nature_all(:,1))],...
%         [slope_high(1)*min(Nature_all(:,1))+ slope_high(2) slope_high(1)*max(Nature_all(:,1)) + slope_high(2)],'r','LineWidth',3);
% 
%     legend(sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_low(1),min(plot_tt(:,3))),...
%            sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_high(1),max(plot_tt(:,3))));

%%
% % time to target
% [nlow,xout] = hist(kin.ttt_L,0.1:0.1:2);
% nhigh = hist(kin.ttt_H,0.1:0.1:2);
% 
% figure; hold on;
% bar(xout-.025,nlow,'b');
% bar(xout+.025,nhigh,'r');
% title('Time To Target','FontSize',18);
% xlabel('Time of Movement (s)');
% ylabel('Count');
% legend('Low variance','High variance');
% 
% % Path length
% [pllow,xplout] = hist(kin.path_L,target_rad:0.5:2*target_rad);
% plhigh = hist(kin.path_H,target_rad:0.5:2*target_rad);
% 
% figure; hold on;
% bar(xplout-.1,pllow,'b');
% bar(xplout+.1,plhigh,'r');
% title('Hand Path Distance','FontSize',18);
% xlabel('Hand Path Distance (cm)');
% ylabel('Count');
% legend('Low variance','High variance');
