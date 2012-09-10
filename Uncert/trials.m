pof = 3; % Position of Feedback (Where in y the cursor appears)

x_off = -2; % Screen offsets
y_off = 32.5;

BDF = bdf;
BDF.pos(:,2) = BDF.pos(:,2) + x_off;
BDF.pos(:,3) = BDF.pos(:,3) + y_off; 

comp_tt = tt(tt(:,8)==32 | tt(:,8)==34,:);

low_var = find(comp_tt(:,3)== min(comp_tt(:,3)));
high_var = find(comp_tt(:,3)== max(comp_tt(:,3)));

low_trials = comp_tt(low_var,6:7);
high_trials = comp_tt(high_var,6:7);

figure; hold on;

PC_low = zeros(length(low_var),1);
for i = 1:length(low_var)
    % find indices
    start_v = find(BDF.vel(:,1) < low_trials(i,1),1,'last');
    last_v = find(BDF.vel(:,1) > low_trials(i,2),2,'first');
       
    vx = BDF.vel(start_v:last_v,2);
    vy = BDF.vel(start_v:last_v,3);
    v = sqrt(vx.^2 + vy.^2);
    %plot(v,'b');
    
    pert = comp_tt(low_var(i),2); 
    px = BDF.pos(start_v:last_v,2) + pert;
    py = BDF.pos(start_v:last_v,3);
    
    mid_p = find(py > pof , 1, 'first');
    init_shift = px(mid_p);
    final_shift = px(end);
    
    nec_correction = -init_shift;
    act_correction = final_shift - init_shift;
    perc_corr = act_correction/nec_correction;
    
    PC_low(i) = perc_corr;
    
    plot([px(mid_p) px(end)],[py(mid_p) py(end)],'b');
end

PC_high = zeros(length(high_var),1);
for i = 1:length(high_var)
    % find indices
    start_v = find(BDF.vel(:,1) < high_trials(i,1),1,'last');
    last_v = find(BDF.vel(:,1) > high_trials(i,2),2,'first');
  
    vx = BDF.vel(start_v:last_v,2);
    vy = BDF.vel(start_v:last_v,3);
    v = sqrt(vx.^2 + vy.^2);
    
    
    pert = comp_tt(high_var(i),2); 
    px = BDF.pos(start_v:last_v,2) + pert;
    py = BDF.pos(start_v:last_v,3);
    
    mid_p = find(py > pof , 1, 'first');
    init_shift = px(mid_p);
    final_shift = px(end);
    
    nec_correction = -init_shift;
    act_correction = final_shift - init_shift;
    perc_corr = act_correction/nec_correction;
    
    PC_high(i) = perc_corr;
    
    plot([px(mid_p) px(end)],[py(mid_p) py(end)],'r');
        
end

%figure; plot(ones(length(low_var),1),PC_low,'b.'); hold on;
%plot(2*ones(length(high_var),1),PC_high,'r.');
