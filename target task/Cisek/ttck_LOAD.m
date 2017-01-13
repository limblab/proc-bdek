%-% Setup/User Input
computer = 'work';
monkey = 'Mihili';
MO = '07';  %%7-19-2013
DA = '08';
YE = '2014';
TASK = 'UCK';
fileindx = '001-01';

sortedM1 = 1;
Behavior_plots = 1;
behavior_only = 0; % 1 if only loading behavior

%- BEHAVIOR PARAMETERS ----------------------------------------------------
% Screen offsets
x_off = -3; 
y_off = 33;

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
    M1_file = sprintf('bdf_%s_M1SORTED_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
    load(M1_file);
    bdfM = bdf; clear bdf;
else
    M1_file = sprintf('bdf_%s_M1_%s_%s_%s.mat',monkey,modaye,TASK,fileindx(1:3));
    load(M1_file);
    bdfM = bdf; clear bdf;
end

% Load PMd bdf
if behavior_only ~= 1
    PMd_file = sprintf('bdf_%s_PMDSORTED_%s_%s_%s.mat',monkey,modaye,TASK,fileindx);
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
%%
%-% Get trial table
fprintf('Getting behavior information...\n');

[tt,~,labels] = getTT_CK(bdfM);
%-% Drop meaningless ratio values
weird_vals = find(abs(tt(:,14))>10000 | abs(tt(:,15))>10000);
tt(weird_vals,:) = [];

%-% Drop trials with missing timestamps
bad_times = find(prod(diff(tt(:,[4 5 6 7 9 10]),[],2),2)<0);
tt(bad_times,:) = [];

%-% Calculate endpoints
endpos = zeros(size(tt,1),1);
for i = 1:size(tt,1)
    endts = tt(i,10);
    indend = find(bdfM.pos(:,1)>endts,1,'first');
    endpos(i) = atan2(bdfM.pos(indend,3)+y_off,bdfM.pos(indend,2)+x_off);
end
tt = [tt endpos];

fprintf('Correcting for photo-delay and catch trials...\n');
%-% Photo Correction
columns2fix = [4:7 10];
if strcmp(YE,'2013') && str2double(MO)<8
    average_lag = 81;
else
    average_lag = 96;
end
orig_tt = tt; 
[corrected_ts,lags] = photo_correct(tt,bdfM,columns2fix,'mean');

figure; hist(1000*lags(:),100); xlabel('Lag (ms)','FontSize',16); 
ylabel('Count','FontSize',16); 
title('Photo Lags','FontSize',18);
keeplag = input('Are lags good? y/n [y]: ','s');
if isempty(keeplag); keeplag = 'y'; end
switch keeplag
    case 'n'
        manlag = input('Input average lag in ms (leave empty for known average): ');
        if isempty(manlag)
            fprintf('Using average lag of %d ms\n',average_lag);
            tt(:,columns2fix) = tt(:,columns2fix) + average_lag/1000;
        else
            tt(:,columns2fix) = tt(:,columns2fix) + manlag/1000;
        end
        
    case 'y'
        tt(:,columns2fix) = corrected_ts; 
end

%-% Identify center-out trials and fill in block_rats variable
trial_rat = zeros(size(tt,1),2);
for i = 1:size(tt,1)
    if round(10000*(10*tt(i,14) - round(10*tt(i,14))))==1 
        trial_rat(i,:) = [10000 10000];
        tt(i,14:15) = 10000;
    else
        trial_rat(i,:) = tt(i,14:15);
    end
end

%-% Add column for correct target
cor1 = find(abs(tt(:,2)==tt(:,13)));
cor2 = find(abs(tt(:,3)==tt(:,13)));

tt(cor1,18) = tt(cor1,14) + 1;
tt(cor2,18) = tt(cor2,14) + 2;

missed = find(tt(:,18)==0);

%-% Add column for chosen target
distances = abs(circ_dist(repmat(tt(:,17),1,2),tt(:,[2 3])));
choice = (distances(:,2)<distances(:,1)) + 1;

for i = 1:size(tt,1)
    tt(i,19) = choice(i) + trial_rat(i,choice(i));
end
%-% Cut into blocks based on target probability
tt(missed,:) = [];
switches = [0; find(diff(tt(:,14))~=0); size(tt,1)];
alldays = struct;
alldaysC = struct;
for i = 1:(length(switches)-1)
    blockinds = (switches(i)+1):switches(i+1);
    length(blockinds)
    alldays(i).tt = tt(blockinds,:);
    
    alldays(i).ratios = trial_rat(blockinds(i),:);
    
    catches = find(alldays(i).tt(:,8)<0);
    
    alldaysC(i).tt = alldays(i).tt(catches,:);
    alldays(i).tt(catches,:) = [];
    
    alldaysC(i).ratios = alldays(i).ratios;
    
end

alldays(1).PMd_units = PMd_units;
alldays(1).M1_units = M1_units;
alldays(1).labels = [labels; '17: Endpoint'; '18: Correct Target'; '19: Chosen Target'];
alldays(1).kin.pos = bdfM.pos; 
alldays(1).kin.pos(:,2) = alldays(1).kin.pos(:,2)-3;
alldays(1).kin.pos(:,3) = alldays(1).kin.pos(:,3)+33;
alldays(1).kin.vel = bdfM.pos;

alldaysC(1).PMd_units = PMd_units;
alldaysC(1).M1_units = M1_units;
alldaysC(1).labels = [labels; '17: Endpoint'; '18: Correct Target'; '19: Chosen Target'];

close;
fprintf('Done\n');