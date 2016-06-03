index = 13;
time_length = 200; % In ms
discard_thresh = 1; % In cm
VEL_PLOTS = 0;

unit = trains{index};
len = length(Nature_low);

dist = @(x1,x2) sqrt((x1(1)-x2(1)).^2 + (x1(2)-x2(2)).^2);

firings = cell(length(Nature_low),3);
all_distances = zeros(length(Nature_low),2);
movement_inds = cell(length(Nature_low),3);

for i = 1:length(Nature_low)
    
    p1_ind = low_var(i);
    
%%% Find nearest point (high var)
    p1 = Nature_low(i,:);
    
    p2s = [];
    d1s = [];
    for j = 1:length(Nature_high)
        
        p2s = [p2s; Nature_high(j,:)];
        d1s = [d1s ;dist(p1,p2s(j,:))];
        
    end
    
    p2_i = find(d1s == min(d1s));
    d1 = min(d1s);
    p2_ind = high_var(p2_i);
    
%%% Find nearest point (low var)

    p3s = [];
    d2s = [];
    for j = 1:length(Nature_low)
        
        p3s = [p3s; Nature_low(j,:)];
        d2s = [d2s ;dist(p1,p3s(j,:))];
        
    end
    d2s(d2s==0) = 1000;
    
    p3_i = find(d2s == min(d2s));
    d2 = min(d2s); 
    p3_ind = low_var(p3_i); 
    
    startend_1 = S_M_F_low(i,:);
    startend_2 = S_M_F_high(p2_i,:);
    startend_3 = S_M_F_low(p3_i,:);
    
    % Movement indexes for [trial_start feedback_turns_on trial_end];
    y_2move_1 = BDF.pos(startend_1(2):startend_1(3),3);
    y_2move_2 = BDF.pos(startend_2(2):startend_2(3),3);
    y_2move_3 = BDF.pos(startend_3(2):startend_3(3),3);
    
    dy_1 = y_2move_1(2:end) - y_2move_1(1:end-1);
    dy_2 = y_2move_2(2:end) - y_2move_2(1:end-1);  
    dy_3 = y_2move_3(2:end) - y_2move_3(1:end-1); 
    
    if length(dy_1) < 500
        dy_1 = [dy_1; 1000*ones(500-length(dy_1),1)];
    end
    if length(dy_2) < 500
        dy_2 = [dy_2; 1000*ones(500-length(dy_2),1)];
    end
    if length(dy_3) < 500
        dy_3 = [dy_3; 1000*ones(500-length(dy_3),1)];
    end
    
    low_vel_1 = find(dy_1==min(dy_1(1:500)),1,'first') + startend_1(2);
    low_vel_2 = find(dy_2==min(dy_2(1:500)),1,'first') + startend_2(2);
    low_vel_3 = find(dy_3==min(dy_3(1:500)),1,'first') + startend_3(2);
    

%     figure; 
%     plot(BDF.vel(startend_1(2):startend_1(3),1),BDF.vel(startend_1(2):startend_1(3),3)); hold on;
%     plot(BDF.vel(low_vel_1,1),BDF.vel(low_vel_1,3),'ro'); 
%     
%     figure; 
%     plot(BDF.vel(startend_2(2):startend_2(3),1),BDF.vel(startend_2(2):startend_2(3),3)); hold on;
%     plot(BDF.vel(low_vel_2,1),BDF.vel(low_vel_2,3),'ro');
%     
%     figure;
%     plot(BDF.vel(startend_3(2):startend_3(3),1),BDF.vel(startend_3(2):startend_3(3),3)); hold on;
%     plot(BDF.vel(low_vel_3,1),BDF.vel(low_vel_3,3),'ro');
    distances = [d1 d2];
    all_distances(i,:) = distances; 
    
    movement_inds{i,1} = [low_vel_1 low_vel_1 + time_length];
    movement_inds{i,2} = [low_vel_2 low_vel_2 + time_length];
    movement_inds{i,3} = [low_vel_3 low_vel_3 + time_length];
    
    firing1 = sum(unit > BDF.pos(low_vel_1,1) & unit < BDF.pos(low_vel_1 + time_length, 1));
    firing2 = sum(unit > BDF.pos(low_vel_2,1) & unit < BDF.pos(low_vel_2 + time_length, 1));
    firing3 = sum(unit > BDF.pos(low_vel_3,1) & unit < BDF.pos(low_vel_3 + time_length, 1));
    
    firings{i,1} = firing1; firings{i,2} = firing2; firings{i,3} = firing3;
end

con = [vertcat(firings{:,1}) vertcat(firings{:,3}) vertcat(firings{:,2})];
% [Low_1 Low_2 High]
bads = find(all_distances(:,1) > discard_thresh | all_distances(:,2) > discard_thresh);

con(bads,:) = [];
all_distances(bads,:) = [];

difs21 = con(:,2) - con(:,1);
difs31 = con(:,3) - con(:,1);

[d21,xout] = hist(difs21, min([difs21;difs31]):max([difs21;difs31]));
d31 = hist(difs31, min([difs21;difs31]):max([difs21;difs31]));

figure; hold on; bar(xout,[d21' d31']);
title(sprintf('Unit %.0f: Spiking Differences',index),'FontSize',14);
xlabel('Change in Spikes');
ylabel('Count');
legend(sprintf('Same Condition (Low Var)\nmean: %.3f  var: %.3f',mean(difs21),var(difs21)),...
    sprintf('Opposite Condition (High Var)\nmean: %.3f  var: %.3f',mean(difs31),var(difs31)));

if VEL_PLOTS == 1

    figure; subplot(1,2,1); hold on;
    title('X velocity after low vel');
    for i = 1:length(Nature_low)
        plot(0:time_length, BDF.vel(movement_inds{i,1}(1):movement_inds{i,1}(2),2),'g');
        plot(0:time_length, BDF.vel(movement_inds{i,3}(1):movement_inds{i,3}(2),2),'b');
        plot(0:time_length, BDF.vel(movement_inds{i,2}(1):movement_inds{i,2}(2),2),'r');
        if i == 1
            legend('Low Var 1','Low Var 2','High Var');
        end
    end

    subplot(1,2,2); hold on;
    title('Y velocity after low vel');
    for i = 1:length(Nature_low)
        plot(0:time_length, BDF.vel(movement_inds{i,1}(1):movement_inds{i,1}(2),3),'g');
        plot(0:time_length, BDF.vel(movement_inds{i,3}(1):movement_inds{i,3}(2),3),'b');
        plot(0:time_length, BDF.vel(movement_inds{i,2}(1):movement_inds{i,2}(2),3),'r');
    end

    set(gcf,'NextPlot','add');
    axes;
    h = title(sprintf('Unit: %.0f',index),'FontSize',16);
    set(gca,'Visible','off');
    set(h,'Visible','on');

end