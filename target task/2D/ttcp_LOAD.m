%-% Setup/User Input
computer = 'work';
monkey = 'Mihili';
MO = '04';  %%7-19-2013  %%5-21-2014 %9-28-2015 %10-07-2013
DA = '09';
YE = '2014';
TASK = 'UNT2D'; %UNT2D
fileindx = '001-01';

sortedM1 = 1;
Behavior_plots = 1;
behavior_only = 0; % 1 if only loading behavior

%% Behavior Parameters
shift_m = 5.5; % Shift mean
shift_v = 20; % Shift variance

target_rad = 10; % Target radius
target_size = 2; % Target size

% Screen offsets
x_off = -3-0.16; 
y_off = 33+0.18;

%-% Load files
fprintf('Loading files...\n');
modaye = [MO DA YE];
orig_place = cd;

% Go to folder with BDFs %
if strcmp(computer,'work')
    cd(sprintf('C:\\Users\\limblab\\Desktop\\%s_bdf',monkey));
elseif strcmp(computer,'home')
    cd(sprintf('/Users/brian/Desktop/Northwestern/Data/Mihili/'));
else
    cd('../../');
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

% Load PMd bdf
if behavior_only ~= 1
    PMd_file = sprintf('bdf_sorted_%s_PMd_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(PMd_file);
    bdfP = bdf; clear bdf;
end

% Go back to original folder and clean up
cd(orig_place);
clear M1_file PMd_file fileindx %M1_units PMd_units
%-% Get Unit information
id_func = @(x) x(1);
if behavior_only ~= 1
    if sortedM1==1
        M1_units = spiketrains(bdfM,1);
        M1_ids = cellfun(id_func,M1_units);
    else
        M1_units = [];
        M1_ids = [];
    end
    PMd_units = spiketrains(bdfP,1);
    PMd_ids = cellfun(id_func,PMd_units);
else
    M1_units = [];
    M1_ids = [];
    PMd_units = [];
    PMd_ids = [];
end
%% Get Behavior
fprintf('Getting behavior information...\n');

[tt,tts,priors] = getTT_UNT_circ(bdfM);
ttslengths = cellfun(@(x) size(x,1), tts);
bad_block = find(ttslengths < 2); 
if ~isempty(bad_block); tts(bad_block) = []; priors(bad_block) = []; end

if length(tts) > 5
    clear tt tts priors
%switch datab_format
   
    %case 'old'
    [tt] = getTT_UNT_circ_OLD(bdfM);
    figure; plot(tt(:,2),'.'); title('Use to find prior blocks','FontSize',20);
    drawnow; pause;
    numps = input('How many prior blocks: ');
    
    [tts,priors] = deal(cell(numps,1));
    for pb = 1:numps
        fprintf('\n--- Prior Block %d ---\n',pb);
        
        if pb==1
            firstindofp = 1;
        else
            firstindofp = lastindofp + 1;
        end
        
        %firstindofp = input('First trial index: ');
        if pb==numps
            lastindofp = size(tt,1);
        else
            lastindofp = input('Last trial index: ');
        end
        prival = input('Prior kappa (0 for CO or leave empty for auto-find): ');
        if prival==0; prival=10000; end
        tts{pb} = tt(firstindofp:lastindofp,:);
        if isempty(prival)
            prival = circ_kappa(tts{pb}(:,2));
        end
        
        priors{pb}.val = prival;
        priors{pb}.inds = firstindofp:lastindofp;
    end
    close;
    
%     case 'new'
%     [tt,tts,priors] = getTT_UNT_circ(bdfM);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BDF ,comp_tt, tt_labels, Linds, Hinds, kin, slices] = UNT_circ_behavior(bdfM,tt,x_off,y_off);

[c_tt,full_c_tt,slice_table,full_slice_table] = deal(cell(length(tts),1));

for i = 1:length(tts)
    [~,c_tt{i},~,~,~,~,slice_table{i}] = UNT_circ_behavior(bdfM,tts{i},x_off,y_off);
end

fprintf('Correcting for photo-delay and catch trials...\n');
%% Photo Correction
columns2fix = [5];
if strcmp(YE,'2013') && str2double(MO)<8
    average_lag = 81;
else
    average_lag = 96;
end
orig_comp_tt = comp_tt; 
[corrected_ts,lags] = photo_correct(comp_tt,bdfM,columns2fix,'mean');

figure; hist(1000*lags(:),100); xlabel('Lag (ms)','FontSize',16); 
% if sum(~isnan(lags(:))) > 0.9*length(lags(:)) && 1000*nanmean(lags(:)) > 80 && 1000*nanmean(lags(:)) < 110;
%     keeplag = 'y';
% else
%     keeplag = 'n';
%     monkey = ['REDO' monkey];
% end

ylabel('Count','FontSize',16); 
title('Photo Lags','FontSize',18);
keeplag = input('Are lags good? y/n [y]: ','s');
if isempty(keeplag); keeplag = 'y'; end
switch keeplag
    case 'n'
        manlag = input('Input average lag in ms (leave empty for known average): ');
        if isempty(manlag); manlag = 10000; fprintf('Using average lag of %d ms\n',average_lag); end
        if manlag == 10000
            comp_tt(:,columns2fix) = comp_tt(:,columns2fix) + average_lag/1000;
            if exist('c_tt','var')
                o_c_tt = cell(length(c_tt),1);
                c_ts = cell(length(c_tt),1);
                for i = 1:length(c_tt)
                    o_c_tt{i} = c_tt{i};
                    c_ts{i} = o_c_tt{i};
                    c_ts{i}(:,columns2fix) = c_ts{i}(:,columns2fix) + average_lag/1000;
                end
            end
        else
            comp_tt(:,columns2fix) = comp_tt(:,columns2fix) + manlag/1000;
            if exist('c_tt','var')
                o_c_tt = cell(length(c_tt),1);
                c_ts = cell(length(c_tt),1);
                for i = 1:length(c_tt)
                    o_c_tt{i} = c_tt{i};
                    c_ts{i} = o_c_tt{i};
                    c_ts{i}(:,columns2fix) = c_ts{i}(:,columns2fix) + manlag/1000;
                end
            end
        end
    case 'y'
        comp_tt(:,columns2fix) = corrected_ts; 
        if exist('c_tt','var')
            o_c_tt = cell(length(c_tt),1);
            c_ts = cell(length(c_tt),1);
            for i = 1:length(c_tt)
                o_c_tt{i} = c_tt{i};
                [c_ts{i}] = photo_correct(c_tt{i},bdfM,columns2fix,'mean');
                c_tt{i}(:,columns2fix) = c_ts{i};
            end
        end
end
%% Drop Catch trials
[nocatch_tt, catch_tt] = dropCatchTrials(comp_tt,0.5);
comp_tt = nocatch_tt; 
good_inds = cell(length(tts),1);

if exist('c_tt','var')

    good_tt = cell(length(c_tt),1);
    bad_tt = cell(length(c_tt),1);
    for i = 1:length(c_tt)

       [good_tt{i},bad_tt{i},good_inds{i}] = dropCatchTrials(c_tt{i},0.5);
       full_c_tt{i} = c_tt{i};
       
       c_tt{i} = good_tt{i};
       slice_table{i} = slice_table{i}(good_inds{i},:);
       
       full_c_tt{i}(~good_inds{i},3) = full_c_tt{i}(~good_inds{i},3)+0.99;
       
       full_slice_table{i} = slice_table{i};

   end
end
%% Set up new variables for cross-day comparison
if ~exist('session_nokin','var'); session_nokin = cell(12,31); end
if ~exist('bdfP','var'); bdfP = 0; end
if ~exist('PMd_units','var'); PMd_units = 0; end

session_nokin{str2double(MO),str2double(DA)}.tt = comp_tt;
session_nokin{str2double(MO),str2double(DA)}.labels = tt_labels;
session_nokin{str2double(MO),str2double(DA)}.PMd_units = PMd_units;
session_nokin{str2double(MO),str2double(DA)}.M1_units = M1_units;
session_full = session_nokin;

if exist('c_tt','var')

    clear session_nokin;
    clear session_full;
    session_nokin = cell(length(c_tt),1);
    session_full = cell(length(c_tt),1);
    for i = 1:length(c_tt)
        if i==1
            session_nokin{i}.PMd_units = PMd_units;
            session_nokin{i}.M1_units = M1_units;
        else
            session_nokin{i}.PMd_units = [];
            session_nokin{i}.M1_units = [];
            
        end
        session_full{i} = session_nokin{i};
        
        session_nokin{i}.tt = c_tt{i};
        session_nokin{i}.slices = slice_table{i};
        
        session_full{i}.tt = full_c_tt{i};
        session_full{i}.slices = full_slice_table{i};
    end
end
alldays_nokin = vertcat(session_nokin{:}); % Reduce 'session' by concatenating
alldays_full = vertcat(session_full{:});
for blcks = 1:length(alldays_nokin); 
    alldays_nokin(blcks).slices(:,sum(isnan(alldays_nokin(blcks).slices))==size(alldays_nokin(blcks).slices,1)) = []; 
end
%-% Set up no-kin variables for cross-day comparison
if ~exist('session','var'); session = cell(12,31); end
if ~exist('bdfP','var'); bdfP = 0; end
if ~exist('PMd_units','var'); PMd_units = 0; end

session{str2double(MO),str2double(DA)}.bdfM = bdfM;
session{str2double(MO),str2double(DA)}.bdfP = bdfP;
session{str2double(MO),str2double(DA)}.tt = comp_tt;
session{str2double(MO),str2double(DA)}.PMd_units = PMd_units;
session{str2double(MO),str2double(DA)}.M1_units = M1_units;
session{str2double(MO),str2double(DA)}.kin = BDF;
session{str2double(MO),str2double(DA)}.vars.Linds = Linds;
session{str2double(MO),str2double(DA)}.vars.Hinds = Hinds;
session{str2double(MO),str2double(DA)}.slices = slices;
session{str2double(MO),str2double(DA)}.labels = tt_labels;


kinBDF.pos = BDF.pos;
kinBDF.vel = BDF.vel;
kinBDF.acc = BDF.acc;
if exist('c_tt','var')

    clear session;
    session = cell(length(c_tt),1);
    
    session{1}.tt = c_tt{1};
    session{1}.slices = slice_table{1};
    %session{1}.bdfM = bdfM;
    %session{1}.bdfP = bdfP;
    session{1}.PMd_units = PMd_units;
    session{1}.M1_units = M1_units;
    session{1}.kin = kinBDF;
    session{1}.vars.Linds = Linds;
    session{1}.vars.Hinds = Hinds;
    session{1}.labels = tt_labels;
    
    for i = 2:length(c_tt)
        session{i}.tt = c_tt{i};
        session{i}.slices = slice_table{i};
        
        %session{i}.bdfM = [];
        %session{i}.bdfP = [];
        session{i}.PMd_units = [];
        session{i}.M1_units = [];
        session{i}.kin = [];
        session{i}.vars.Linds = [];
        session{i}.vars.Hinds = [];
        session{i}.labels = [];
    end
end

alldays = vertcat(session{:}); % Reduce 'session' by concatenating
for blcks = 1:length(alldays); 
    alldays(blcks).slices(:,sum(isnan(alldays(blcks).slices))==size(alldays(blcks).slices,1)) = []; 
end
close;
fprintf('Done\n');