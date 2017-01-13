%--------------------------------------------------------------------------
dump = [];
BDF = bdfM;
BDF.pos(:,2) = BDF.pos(:,2) + x_off;
BDF.pos(:,3) = BDF.pos(:,3) + y_off; 

tt(isnan(tt(:,5)),:) = [];
tt(tt(:,3)<0,:) = [];

comp_tt = tt(tt(:,8)==32 | tt(:,8)==34,:);

% Low indices in comp_tt
low_var_inds = find(comp_tt(:,3)== min(comp_tt(:,3)));
high_var_inds = find(comp_tt(:,3)== max(comp_tt(:,3)));

L_start_times = comp_tt(low_var_inds,6);
H_start_times = comp_tt(high_var_inds,6);

L_end_times = comp_tt(low_var_inds,7);
H_end_times = comp_tt(high_var_inds,7);

% Target Shifts
low_var_shifts = comp_tt(low_var_inds,2);
high_var_shifts = comp_tt(high_var_inds,2);

% Time to target
low_var_ttt = comp_tt(low_var_inds,7)-comp_tt(low_var_inds,6);
high_var_ttt = comp_tt(high_var_inds,7)-comp_tt(high_var_inds,6);

low_var_pathdist = zeros(length(low_var_inds),1);
high_var_pathdist = zeros(length(high_var_inds),1);

% Get slice positions
bursttimes = vertcat(bdfM.databursts{:,1});
slice_x = zeros(length(comp_tt),numslices);
for i = 1:length(comp_tt)
    burst = bdfM.databursts{comp_tt(i,1)==bursttimes,2};
    xypos = bytes2float(burst(end-(4*numslices*2-1):end));
    slice_x(i,:) = xypos(1:2:end);
    
    centroid(i) = mean(slice_x(i,:));
end
low_cents = centroid(low_var_inds);
high_cents = centroid(high_var_inds);
% Find endpoint cursor positions
L_end_pos = zeros(length(low_var_inds),1);
for i = 1:length(low_var_inds)
    bdf_start_ind = find(BDF.pos(:,1)>L_start_times(i),1,'first');
    bdf_end_ind = find(BDF.pos(:,1)<L_end_times(i),1,'last');
    
    L_end_pos(i) = BDF.pos(bdf_end_ind,2); % final position
    
    movex = BDF.pos(bdf_start_ind:bdf_end_ind,2);
    movey = BDF.pos(bdf_start_ind:bdf_end_ind,3);
    
    distx = sum(abs(diff(movex)));
    disty = sum(abs(diff(movey)));
    
    low_var_pathdist(i) = distx+disty;
    
end
H_end_pos = zeros(length(high_var_inds),1);
for i = 1:length(high_var_inds)
    bdf_start_ind = find(BDF.pos(:,1)>H_start_times(i),1,'first');
    bdf_end_ind = find(BDF.pos(:,1)<H_end_times(i),1,'last');
    
    H_end_pos(i) = BDF.pos(bdf_end_ind,2); %final position
    
    movex = BDF.pos(bdf_start_ind:bdf_end_ind,2);
    movey = BDF.pos(bdf_start_ind:bdf_end_ind,3);
    
    distx = sum(abs(diff(movex)));
    disty = sum(abs(diff(movey)));
    
    high_var_pathdist(i) = distx+disty;
end

if NATURE_PLOT == 1
    
    % Nature Plot
    badlow = find(abs(low_cents)>100);
    badhigh = find(abs(high_cents)>100);
    
    Nature_low = [low_cents' L_end_pos];
    Nature_high = [high_cents' H_end_pos];
    
    low_cents(badlow) = [];
    high_cents(badhigh) = [];
    L_end_pos(badlow) = [];
    H_end_pos(badhigh) = [];
    Nature_high(badhigh,:) = [];
    Nature_low(badlow,:) = [];

    
    low_cents(dump) = [];
    high_cents(dump) = [];
    L_end_pos(dump) = [];
    H_end_pos(dump) = [];
    Nature_high(dump,:) = [];
    Nature_low(dump,:) = [];
    
    Nature_all = [Nature_low;Nature_high];
    
    [BL,NatureLsort] = sortrows(Nature_low(:,1));
    [BH,NatureHsort] = sortrows(Nature_high(:,1));

    figure; plot(low_cents,Nature_low(:,2),'b.','MarkerSize',4);
    hold on; title(sprintf('%s:  %s/%s/%s',monkey,MO,DA,YE),'FontSize',18);
    xlabel('Centroid Location'); ylabel('End Cursor Location');
    plot(high_cents,Nature_high(:,2),'r.','MarkerSize',4);

    [slope_low,SL] = polyfit(Nature_low(:,1),Nature_low(:,2),1);
    [slope_high,SH] = polyfit(Nature_high(:,1),Nature_high(:,2),1);
    
    [YL,DELTAL] = polyconf(slope_low,Nature_low(:,1),SL,'predopt','curve');
    [YH,DELTAH] = polyconf(slope_high,Nature_high(:,1),SH,'predopt','curve');

    %plot([min(Nature_all(:,1)) max(Nature_all(:,1))],...
        %[slope_low(1)*min(Nature_all(:,1))+ slope_low(2) slope_low(1)*max(Nature_all(:,1)) + slope_low(2)],'b','LineWidth',3);
    
    plot(Nature_low(:,1),YL,'b','LineWidth',5);
    plot(Nature_high(:,1),YH,'r','LineWidth',5);
    
    legend(sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_low(1),min(comp_tt(:,3))),...
    sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_high(1),max(comp_tt(:,3))));
    
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
%     legend(sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_low(1),min(comp_tt(:,3))),...
%            sprintf('Fit Slope = %.3f\nCloud var = %.2f',slope_high(1),max(comp_tt(:,3))));
       
    axis equal
    
    % time to target
    [nlow,xout] = hist(low_var_ttt,0.1:0.1:2);
    nhigh = hist(high_var_ttt,0.1:0.1:2);
    
    figure; hold on;
    bar(xout-.025,nlow,'b');
    bar(xout+.025,nhigh,'r');
    title('Time To Target','FontSize',18);
    xlabel('Time of Movement (s)');
    ylabel('Count');
    legend('Low variance','High variance');
    
        
    % Path length
    [pllow,xplout] = hist(low_var_pathdist,target_rad:0.5:2*target_rad);
    plhigh = hist(high_var_pathdist,target_rad:0.5:2*target_rad);
    
    figure; hold on;
    bar(xplout-.1,pllow,'b');
    bar(xplout+.1,plhigh,'r');
    title('Hand Path Distance','FontSize',18);
    xlabel('Hand Path Distance (cm)');
    ylabel('Count');
    legend('Low variance','High variance');
    
end