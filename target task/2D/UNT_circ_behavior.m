function[BDF,comp_tt,tt_labels,Linds,Hinds,kin,slices] = UNT_circ_behavior(bdfM,tt,x_off,y_off)
%--------------------------------------------------------------------------
BDF = bdfM;
BDF.pos(:,2) = BDF.pos(:,2) + x_off;
BDF.pos(:,3) = BDF.pos(:,3) + y_off; 

%tt(isnan(tt(:,5)),:) = [];
badlikes = find(tt(:,3)<0);
for i = 1:length(badlikes)
    bl = badlikes(i);   
    tt(bl,3) = tt(bl-1,3);
end
    
%tt(tt(:,3)<0,:) = []; 

comp_tt =  tt(ismember(tt(:,8),[32 34]),:);
%comp_tt(comp_tt(:,3)<0.01,:) = [];

%% Get slice Angles
bursttimes = vertcat(bdfM.databursts{:,1});
centroid = zeros(size(comp_tt,1),1);
slice_k = zeros(size(comp_tt,1),1);
slices = nan(size(comp_tt,1),10);
nsls = zeros(size(comp_tt,1),1);

for i = 1:size(comp_tt,1)
    if ~isnan(comp_tt(i,1))
        burst = bdfM.databursts{comp_tt(i,1)==bursttimes,2};
    else
        burst = nan(size(bdfM.databursts{1,2}));
    end
    
    if length(burst) == 61 
        newburst = zeros(1,65); 
        newburst([1:17 22:65]) = burst; 
        burst = newburst;
    end
    nsls(i) = bytes2float(burst(22:25)); 
    
    if ~ismember(nsls(i),0:10); 
        if i > 1
            nsls(i) = nsls(i-1); 
        else
            nsls(i) = 1;
        end
    end
    slice_angs = comp_tt(i,2)+bytes2float(burst(26:65));
    %slice_angs = comp_tt(i,2)+bytes2float(burst(end-(2*numslices*2-1):end));
    disp_angs = slice_angs(1:nsls(i));
    
    if sum(abs(disp_angs) > 1000) > 0; disp_angs = nan; end;
    
    slices(i,1:nsls(i)) = disp_angs;
    
    
    if nsls(i) > 1
        if max(abs(disp_angs))>4*pi
            centroid(i) = NaN;
            slice_k(i) = NaN;
        elseif sum(isnan(disp_angs))==0
            centroid(i) = circ_mean(disp_angs);
            %centroid(i) = mean_angle(slice_angs,'rads');
            %slice_k(i) = kappa_calc(slice_angs,0.05);
            %[~,slice_k(i)] = circ_var(slice_angs);
            slice_k(i) = calc_kappa(disp_angs);
        else
            centroid(i) = NaN;
            slice_k(i) = NaN;
        end

    else

        if max(abs(disp_angs))>4*pi
            centroid(i) = nan;
            slice_k(i) = nan;
        else
            centroid(i) = disp_angs;
            slice_k(i) = nan;
        end
    end
end
% Eliminate junk
bad_trials = [];%isnan(centroid);
comp_tt(bad_trials,:) = [];
centroid(bad_trials) = [];
slice_k(bad_trials) = [];
slices(bad_trials,:) = [];

%% Time to target
% Low/High indices in comp_tt
Linds = find(comp_tt(:,3)== max(comp_tt(:,3)));
Hinds = find(comp_tt(:,3)== min(comp_tt(:,3)));

start_times = comp_tt(:,6);
end_times = comp_tt(:,7);

kin.ttt_L = comp_tt(Linds,7)-comp_tt(Linds,6);
kin.ttt_H = comp_tt(Hinds,7)-comp_tt(Hinds,6);

kin.path = zeros(size(comp_tt,1),1);

%% Find endpoint cursor positions/angles
end_ang = zeros(size(comp_tt,1),1);
for i = 1:size(comp_tt,1)
    bdf_start_ind = find(BDF.pos(:,1)>start_times(i),1,'first');
    bdf_end_ind = find(BDF.pos(:,1)<end_times(i),1,'last');
    
    end_pos = BDF.pos(bdf_end_ind,2:3); % final position
    end_ang(i) = atan2(end_pos(2),end_pos(1));
    
    movex = BDF.pos(bdf_start_ind:bdf_end_ind,2);
    movey = BDF.pos(bdf_start_ind:bdf_end_ind,3);
    
    distx = sum(abs(diff(movex)));
    disty = sum(abs(diff(movey)));
    
    kin.path(i) = distx+disty;

end

%% Adjust trial table
comp_tt = [comp_tt zeros(size(comp_tt,1),2)];
% Add centroid
comp_tt(:,9) = mod(centroid,2*pi);
% Add guess
comp_tt(:,10)= mod(end_ang,2*pi);
% Add slice kappa
comp_tt(:,11) = slice_k;

tt_labels = {'DB time','Target Loc','Var','CT ON','OT ON','GO','END','Result','Centroid','Estimate','Slice Kappa'};
