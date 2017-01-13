%% Setup/User Input
computer = 'work';
monkey = 'Mihili';
MO = '07';
DA = '10';
YE = '2013';
TASK = 'CO';
fileindx = '001-01';
sortedM1 = 1;
sortedPMd = 0;
save_plots = 0; % 1 to save plots (check "Spiking Behavior" for directory)

%- BEHAVIOR PARAMETERS ----------------------------------------------------
shift_m = 0; % Shift mean
shift_v = 20; % Shift variance
numslices = 10;

NATURE_PLOT = 1;
yoff = -3; % Global y offset of center target
target_rad = 10; % Target radius
target_size = 2; % Target size

% Screen offsets
x_off = -2; 
y_off = 32.5;

%% Load files
modaye = [MO DA YE];
orig_place = cd;

% Go to folder with BDFs %
if strcmp(computer,'work')
    cd(sprintf('C:\\Users\\limblab\\Desktop\\%s_bdf',monkey));
else   
    cd('../');
    cd([cd '/bdf'])
end

% Load M1/Behavior bdf
if sortedM1==1
    M1_file = sprintf('bdf_sorted_%s_M1_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(M1_file);
    bdfM = bdf; clear bdf;
else
    M1_file = sprintf('bdf_%s_M1_%s_%s_%s.mat',monkey,modaye,TASK,fileindx(1:3));
    load(M1_file);
    bdfM = bdf; clear bdf;
end

if sortedPMd ==1
    % Load PMd bdf
    PMd_file = sprintf('bdf_sorted_%s_PMD_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(PMd_file);
    bdfP = bdf; clear bdf;
end

% Go back to original folder and clean up
cd(orig_place);
clear M1_file PMd_file fileindx modaye TASK
%% Get Unit information
if sortedM1==1
    M1_units = spiketrains(bdfM,1);
end

if sortedPMd==1
    PMd_units = spiketrains(bdfP,1);
end
%% Get Behavior
tt = getTT_UNT1D(bdfM);
UNT1D_behavior;

%% Spiking Behavior

t1 = 0;
t2 = 0.8;

goodPMD_1232013 = [1 2 5 6 9 10 11 15 16 17 22 24 25 26 27 29 30 32 36 37 43 44];
goodPMD_2122013 = [3 4 6 7 8 11 12 13 15 16 17 19 21 22 24 25 26 27 29 ...
           30 33 36 37 38 39 41 42 43 45 46 47 48 49 50 51 52 57 60 63 66];

% 1_23_2013 
TARGET_PMD = [15 16 19 22 24 25 26 27 32 36 41 42 44];
PLAN_PMD = [2 6 9 11 17 26 27 29 37];
MOVE_PMD = [5 6 10 11 17 35 41];

allunits = 1:length(PMd_units);

use_set = 1:length(PMd_units);
for neu = use_set
    neuron_ind = find(use_set==neu);
    if save_plots==1
        tic;
        close;
    end

    train = PMd_units{neu};
    [OUTPUTS] = raster_UNT1D(BDF,train,comp_tt,t1,t2);
        %OUTPUTS.x = [final_pos, target_pos, maxspeed, maxspeed_time, direction, UNC];
        %OUTPUTS.y = [mean_fr];    
    title(sprintf('PMd: %d (%.1f) (Target)',neu,train(1)),'FontSize',18);

    if save_plots ==1
        h = gca;
        cd('C:\Users\limblab\Desktop\figures\2_17_2013\Target Aligned\');
        saveas(h,sprintf('PMd_%s_target_%d (%d_%d)',...
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
%% GLM fits



