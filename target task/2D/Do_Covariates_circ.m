function [X,spikes,Xnames,trial_inds]=Do_Covariates_circ(pva,Spiketrains,comp_tt,BINSIZE,Uncert,BASIS)

%% Extract Temporal basis function sets
acBasis = BASIS{1};
cBasis = BASIS{2};
ckinBasis = BASIS{3};
aceventBasis = BASIS{4};
ceventBasis = BASIS{5};

%% Extract Input (X) and output (Y)
% Kinematics
speedvec = sqrt(pva.vel(:,3).^2 + pva.vel(:,2).^2);
thetavec = atan2(pva.vel(:,3), pva.vel(:,2));
        %accvec = sqrt(pva.acc(:,2).^2 + pva.acc(:,3).^2);

% Downsample to account for bin size
speedvec = downsample(speedvec,BINSIZE);
thetavec = downsample(thetavec,BINSIZE);
        %accvec = downsample(accvec,BINSIZE);
cosvec = cos(thetavec);
sinvec = sin(thetavec);
dirvec = [cosvec sinvec];
vx = speedvec.*dirvec(:,1);
vy = speedvec.*dirvec(:,2);

% Assemble spikes
spiketimes = 1000*Spiketrains;
spikes = histc(spiketimes, 1000*pva.vel(1,1):BINSIZE:1000*pva.vel(end,1));

%% Initialize
%Uncertainty Vectors
LOvec = zeros(length(pva.vel),1);
HIvec = zeros(length(pva.vel),1);
%Events
Tgo = zeros(ceil(length(pva.vel)/BINSIZE),1);
Fon_L = zeros(ceil(length(pva.vel)/BINSIZE),1);
Fon_H = zeros(ceil(length(pva.vel)/BINSIZE),1);
Tend = zeros(ceil(length(pva.vel)/BINSIZE),1);

trial_inds.wait = zeros(ceil(length(pva.vel)/BINSIZE),1);
trial_inds.move = zeros(ceil(length(pva.vel)/BINSIZE),1);
trial_inds.all = zeros(ceil(length(pva.vel)/BINSIZE),1);
%% Fill in Events Vectors
% First extract trial information
theta_shift = comp_tt(:,2);
target_on = 1000*comp_tt(:,5);
go_cue = 1000*comp_tt(:,6);
trialend = 1000*comp_tt(:,7);

for tr=1:length(comp_tt) 
    % mark trial indices
    trial_inds.wait(round((target_on(tr)-999)./BINSIZE):round((go_cue(tr)-999)./BINSIZE)) = tr;
    trial_inds.move(round((go_cue(tr)-999)./BINSIZE):round((trialend(tr)-999)./BINSIZE)) = tr;
    trial_inds.all(round((target_on(tr)-999)./BINSIZE):round((trialend(tr)-999)./BINSIZE)) = tr;
    
    % Add constant to either LO or HI uncertainty vector
    if(comp_tt(tr,3) == max(comp_tt(:,3))) % LOW Uncertainty (Kappa high)
        LOvec(round(target_on(tr)-999):round(trialend(tr)-999)) = 1; 
        Fon_L(round((target_on(tr)-999)./BINSIZE)) = 1; % Delta function (target on)
    else % HIGH Uncertainty (Kappa low)
        HIvec(round(target_on(tr)-999):round(trialend(tr)-999)) = 1;
        Fon_H(round((target_on(tr)-999)./BINSIZE)) = 1;  % Delta function (target on)
    end
    % Mark events with delta function
    Tgo(round((go_cue(tr)-999)./BINSIZE)) = 1; % Go cue
    Tend(round((trialend(tr)-999)./BINSIZE)) = 1; % End of trial
end

% Assemble vectors and downsample as necessary
Fon = Fon_L + Fon_H;
LOvec = downsample(LOvec,BINSIZE);
HIvec = downsample(HIvec,BINSIZE);
ANYvec = LOvec + HIvec;

%% Set temporal basis functions
nbasis.temp.ac=size(BASIS{1},2);        % Motor Neuronal activity
nbasis.temp.c=size(BASIS{2},2);         % Visual feedback
nbasis.temp.ckin = size(BASIS{3},2);    % Proprio/somato feedback
nbasis.temp.acevent = size(BASIS{4},2); % Anticipatory event activity
nbasis.temp.cevent = size(BASIS{5},2);  % Reactionary event activity

%% Assign predictors to desired basis functions
if Uncert==1
    % Motor activity
    COV.AC    = {vx,vy,speedvec,cosvec,sinvec}; 
    Xnames.ac = {'Vx','Vy','Speed','Sin','Cos'};
    
    % Visual feedback
    COV.C    = {ANYvec};                        
    Xnames.c = {'ANYvec'};
    
    % Proprio/somato feedback
    COV.Ckin    = {speedvec};                 
    Xnames.ckin = {'Speed_b'};
    
    % Anticipatory event activity
    COV.ACevent    = {};                      
    Xnames.acevent = {};
    
    % Reactionary event activity
    COV.Cevent    = {Fon,Tgo};                
    Xnames.cevent = {'Feedback','Go Cue'};

elseif Uncert==2
    % Motor activity
    COV.AC    = {vx,vy,speedvec,cosvec,sinvec}; 
    Xnames.ac = {'Vx','Vy','Speed','Sin','Cos'};
    
    % Visual feedback
    COV.C    = {LOvec,HIvec};                        
    Xnames.c = {'Low Unc','High Unc'};
    
    % Proprio/somato feedback
    COV.Ckin    = {speedvec};                 
    Xnames.ckin = {'Speed_b'};
    
    % Anticipatory event activity
    COV.ACevent    = {};                      
    Xnames.acevent = {};
    
    % Reactionary event activity
    COV.Cevent    = {Fon,Tgo};                
    Xnames.cevent = {'Feedback','Go Cue'};
end

%% Set sizes of inputs for different Basis functions
nbasis.sp.ac=size(Xnames.ac,2);
nbasis.sp.c=size(Xnames.c,2);
nbasis.sp.ckin = size(Xnames.ckin,2);
nbasis.sp.acevent = size(Xnames.acevent,2);
nbasis.sp.cevent = size(Xnames.cevent,2);

%Adjust forms
nbs = cell2mat(struct2cell(nbasis.sp));
nbt = cell2mat(struct2cell(nbasis.temp));
nbt(nbs==0)=0;
nbasis.temp = cell2struct(num2cell(nbt),{'ac','c','ckin','acevent','cevent'},1);

% Define lengths for use in indexing
lac = nbasis.sp.ac*nbasis.temp.ac;
lc = nbasis.sp.c*nbasis.temp.c;
lckin = nbasis.sp.ckin*nbasis.temp.ckin;
lacevent = nbasis.sp.acevent*nbasis.temp.acevent;
lcevent = nbasis.sp.cevent*nbasis.temp.cevent;

% Preallocate X
X = zeros(length(speedvec),(lac+lc+lckin+lacevent+lcevent));

% AC (Motor activity)
for j=1:nbasis.sp.ac  
    for i=1:nbasis.temp.ac
         temp=filter(acBasis(:,i), 1, COV.AC{j}(end:-1:1,1));
         X(:,nbasis.temp.ac*(j-1)+i) = temp(end:-1:1);
    end
end
% C (Visual feedback)
for j=1:nbasis.sp.c  
    for i=1:nbasis.temp.c
        % Uncertainty must lead neural activity so filter forward
        X(:,lac + nbasis.temp.c*(j-1)+i) = filter(cBasis(:,i), 1, COV.C{j});
    end
end
% Ckin (Proprio/somato feedback)
for j = 1:nbasis.sp.ckin
    for i = 1:nbasis.temp.ckin
        X(:,lac + lc + ...
            nbasis.temp.ckin*(j-1)+i) = filter(ckinBasis(:,i),1,COV.Ckin{j});
    end
end
% ACevent (Anticipatory event activity)
for j = 1:nbasis.sp.acevent
    for i = 1:nbasis.temp.acevent
        temp=filter(aceventBasis(:,i), 1, COV.ACevent{j}(end:-1:1,1));
        X(:,lac+lc+lckin + nbasis.temp.acevent*(j-1)+i) = temp(end:-1:1);
    end
end
% Cevent (Reactionary event activity)
for j = 1:nbasis.sp.cevent
    for i = 1:nbasis.temp.cevent
        X(:,lac+lc+lckin+lacevent+ nbasis.temp.cevent*(j-1)+i) = filter(ceventBasis(:,i),1,COV.Cevent{j});
    end
end

