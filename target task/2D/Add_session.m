function[session] = Add_session(monkey,mmddyy,task,sortedM1,sortedPMd,show_plots,SESSION)
%%
numslices = 10;

%% Setup/User Input
computer = 'work';

MO = mmddyy(1:2);
DA = mmddyy(3:4);
YE = mmddyy(5:end); %#ok<NASGU>

fileindx = '001-01';

%- Lab PARAMETERS ----------------------------------------------------
% Screen offsets
x_off = -2; 
y_off = 32.5;

%% Load files
orig_place = cd;

% Go to folder with BDFs %
if strcmp(computer,'work')
    cd('C:\Users\limblab\Desktop\MrT_bdf');
else   
    cd('../../');
    cd([cd '/bdf'])
end

% Load M1/Behavior bdf
if sortedM1==1
    M1_file = sprintf('bdf_sorted_%s_M1_%s_%s_%s.mat',monkey,mmddyy,task,fileindx);
    load(M1_file);
    bdfM = bdf; clear bdf;
else
    M1_file = sprintf('bdf_%s_M1_%s_%s_%s.mat',monkey,mmddyy,task,fileindx(1:3));
    load(M1_file);
    bdfM = bdf; clear bdf;
end

% Load PMd bdf
if sortedPMd == 1
    PMd_file = sprintf('bdf_sorted_%s_PMD_%s_%s_%s.mat',monkey,mmddyy,task,fileindx);
    load(PMd_file);
    bdfP = bdf; clear bdf;
else
    bdfP = [];
end

% Go back to original folder and clean up
cd(orig_place);
clear M1_file PMd_file fileindx mmddyy task

%% Get Unit information
if sortedM1==1
    M1_units = spiketrains(bdfM,1);
else
    M1_units = [];
end

if sortedPMd==1  
    PMd_units = spiketrains(bdfP,1);
else
    PMd_units = [];
end

%% Get Behavior
tt = getTT_UNT_circ(bdfM);
[BDF, comp_tt, tt_labels, Linds, Hinds, kin, slices] = UNT_circ_behavior(bdfM,tt,x_off,y_off,numslices);

[corrected_ts,lags] = photo_correct(comp_tt,bdfM,[4 5 6 7],'mean');
comp_tt(:,[4 5 6 7]) = corrected_ts;

if show_plots == 1
    Behavior_Plots_unc_circ;
end

bdf_add_waveform;

%% Set up new variables for cross-day comparison
if isempty(SESSION)
    session = cell(12,31);
else
    session = SESSION;
end

session{str2double(MO),str2double(DA)}.bdfM = bdfM;
session{str2double(MO),str2double(DA)}.bdfP = bdfP;
session{str2double(MO),str2double(DA)}.tt = comp_tt;
session{str2double(MO),str2double(DA)}.PMd_units = PMd_units;
session{str2double(MO),str2double(DA)}.M1_units = M1_units;
session{str2double(MO),str2double(DA)}.kin = BDF;
session{str2double(MO),str2double(DA)}.vars.Linds = Linds;
session{str2double(MO),str2double(DA)}.vars.Hinds = Hinds;
session{str2double(MO),str2double(DA)}.vars.kin = kin;
session{str2double(MO),str2double(DA)}.vars.slices = slices;
session{str2double(MO),str2double(DA)}.vars.lags = lags;

end



