%% Spiking Behavior
sess_ind = [4,5];

range = [-1000 10000];



comp_tt = session{sess_ind(1), sess_ind(2)}.tt;
range_rad = range.*pi./180;

comp_tt = comp_tt(comp_tt(:,10)> range_rad(1) & comp_tt(:,10) < range_rad(2), :);
BDF = session{sess_ind(1),sess_ind(2)}.kin;
PMd_units = session{sess_ind(1),sess_ind(2)}.PMd_units;

save_plots = 0;
do_raster_plot = 0;

t1 = -0.75;
t2 = 0.75;

%goodPMD_1232013 = [1 2 5 6 9 10 11 15 16 17 22 24 25 26 27 29 30 32 36 37 43 44];
%goodPMD_2122013 = [3 4 6 7 8 11 12 13 15 16 17 19 21 22 24 25 26 27 29 ...
%           30 33 36 37 38 39 41 42 43 45 46 47 48 49 50 51 52 57 60 63 66];

% 1_23_2013 
%TARGET_PMD = [15 16 19 22 24 25 26 27 32 36 41 42 44];
%PLAN_PMD = [2 6 9 11 17 26 27 29 37];
%MOVE_PMD = [5 6 10 11 17 35 41];

allunits = 1:length(PMd_units);
use_set = 6;
for neu = use_set
    neuron_ind = find(use_set==neu);
    if save_plots==1
        tic;
        close;
    end

    train = PMd_units{neu};
    [OUTPUTS,Lrast,Hrast] = raster_UNT_circ(BDF,train,comp_tt,t1,t2);
        %OUTPUTS.x = [final_pos, target_pos, maxspeed, maxspeed_time, direction, UNC];
        %OUTPUTS.y = [mean_fr];    
    title(sprintf('PMd: %d (%.1f) (Circ Targ)',neu,train(1)),'FontSize',18);
    xlabel('Time from Target (ms)'); ylabel('Endpoint Location');

    if save_plots ==1
        h = gca;
        cd('C:\Users\limblab\Desktop\figures\4_04_2013\Go_aligned\Rasters\Direction_ranked\');
        saveas(h,sprintf('PMd_%s_circgo_rast_dir_%d (%d_%d)',...
            [MO DA],neu,floor(train(1)),round(10*(train(1)-floor(train(1))))),'png');
        cd(orig_place);
        dt = toc;
        clc;
        fprintf('Neuron %d/%d (Time remaining: %.1f min)\n',neu,length(use_set),...
            toc*(length(use_set)-neuron_ind)./60);
    end
end

%% Raster Plot
if do_raster_plot == 1
    t1 = -.2;
    t2 = .8;

    goodPMD = [1 2 5 6 9 10 11 15 16 17 22 24 25 26 27 29 30 32 36 37 43 44];
    handset = 1:length(PMd_units);%

    for neu = 37%1:length(PMd_units)%handset%goodPMD
    %     close
        train = PMd_units{neu};
        rasterplot(BDF,train,comp_tt,t1,t2);
        title(sprintf('PMd: %d',neu),'FontSize',18);
    %     h = gca;
    %     cd('C:\Users\limblab\Desktop\figures\2_14_2013\Go_align\');
    %     saveas(h,sprintf('PMd_0214_GA_%d',neu),'png');
    %     cd(orig_place);
    end
end