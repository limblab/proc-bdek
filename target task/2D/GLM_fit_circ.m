function [b,dev,stats,yhat,psth,LL,AIC]=GLM_fit_circ(Trials,X,spikes,T,CROSSNUM,BINSIZE,Range,FitOnRange)

opts.dt=0.001;

%% Output (Y)
y = spikes; % Discrete
 
%% Get shuffled trial order information
trial_inds = T.trial_inds;
trialnums = T.trialnums;
shuff_trials = T.shuff_trials;
CVranges=T.CVranges;
hinds = T.highinds;
linds = T.lowinds;

%% Initialize
yhat = cell(CROSSNUM,1);
psth = cell(CROSSNUM,1);

b = cell(CROSSNUM,1);
dev = cell(CROSSNUM,1);
stats = cell(CROSSNUM,1);
%% Cycle through cross-validation
for cv=1:CROSSNUM

    % Set test set for current cross validation
    testSet = shuff_trials(CVranges(cv):CVranges(cv+1));
    
    % For no cross validation, setup train and test sets
    if CROSSNUM==1
        tranSet=testSet';
    else
        tranSet=setdiff(trialnums,testSet)';
    end
    
    % Find indices for training/test sets
    feedbackindst = find(ismember(trial_inds,testSet)); % Training
    feedbackindsr = find(ismember(trial_inds,tranSet)); % Test
    
    % Find indices for HIGH or LOW uncertainty trials
    feedbackH = find(ismember(trial_inds,hinds));
    feedbackL = find(ismember(trial_inds,linds));
    
    %% Establish the Region of Interest from input -----------------------%
    STARTROI = round(Range(1)/BINSIZE);
    ENDROI = round(Range(2)/BINSIZE);
    sizeROI = ENDROI-STARTROI+1;

    regionOI = zeros(length(Trials.idx),sizeROI);
    for i = 1:length(Trials.idx)
        regionOI(i,:) = round(Trials.idx(i,2)/BINSIZE+STARTROI):round(Trials.idx(i,2)/BINSIZE+ENDROI);
    end
    ALLregionOI = regionOI'; ALLregionOI = ALLregionOI(:);
    
    % Region of Interest (HIGH UNCERTAINTY) ------------------------------%
    regionOIH = zeros(length(hinds),sizeROI);
    tbt_yH = zeros(length(hinds),sizeROI);
    for i = 1:length(hinds)
        regionOIH(i,:) = round(Trials.idx(hinds(i),2)/BINSIZE+STARTROI):round(Trials.idx(hinds(i),2)/BINSIZE+ENDROI);
        tbt_yH(i,:) = y(regionOIH(i,:));
    end
    HregionOI = regionOIH'; HregionOI = HregionOI(:);
    psth_yH = mean(tbt_yH,1); % Average across trials (Angles)
    
    % Region of Interest (LOW UNCERTAINTY) ------------------------------%
    regionOIL = zeros(length(linds),sizeROI);
    tbt_yL = zeros(length(linds),sizeROI);
    for i = 1:length(linds)
        regionOIL(i,:) = round(Trials.idx(linds(i),2)/BINSIZE+STARTROI):round(Trials.idx(linds(i),2)/BINSIZE+ENDROI);
        tbt_yL(i,:) = y(regionOIL(i,:));
    end
    LregionOI = regionOIL'; LregionOI = LregionOI(:);
    psth_yL = mean(tbt_yL,1); % Average across trials (Angles)
    
    % Find regions of interest within trials for Train/Test/High Unc/Low Unc 
    ROITinds = feedbackindst(ismember(feedbackindst,ALLregionOI));
    ROIRinds = feedbackindsr(ismember(feedbackindsr,ALLregionOI));
    ROI_H = feedbackH(ismember(feedbackH,HregionOI));
    ROI_L = feedbackL(ismember(feedbackL,LregionOI));
    
    % Test/Train sets
    Xt = X(ismember(trial_inds,testSet),:);
    Xr = X(ismember(trial_inds,tranSet),:);
    
    % Set Inputs for ROI
    Xtf = X(ROITinds,:);
    XH = X(ROI_H,:);
    XL = X(ROI_L,:);
    
    % Set Output
    yr = y(ismember(trial_inds,tranSet),:);
    ytf = y(ROITinds,:);
    yH = y(ROI_H,:);
    yL = y(ROI_L,:);
    
    % IF Fitting only on Region of Interest, set INPUTS/OUTPUTS
    if FitOnRange == 1
        
        Xt = X(ROITinds,:);
        Xr = X(ROIRinds,:);
        Xtf = X(ROITinds,:);
        
        yt = y(ROITinds,:);
        yr = y(ROIRinds,:);

    end

%% GLM 
[b{cv},dev{cv},stats{cv}] = glmfit(Xr,yr,'poisson');

%% Find Predictions from GLM output
yhat{cv}.train = exp([ones(size(Xr,1),1) Xr]*b)+eps;
yhat{cv}.test = exp([ones(size(Xt,1),1) Xt]*b)+eps;
yhat{cv}.all = exp([ones(size(Xtf,1),1) Xtf]*b)+eps;
yhat{cv}.HI = exp([ones(size(XH,1),1) XH]*b)+eps;
yhat{cv}.LO = exp([ones(size(XL,1),1) XL]*b)+eps;

%% Create PSTHs
psth_lamHtbt = reshape(yhat.HI,sizeROI,length(yhat.HI)/sizeROI);
psth{cv}.yHtbt = mean(psth_lamHtbt,2)';

psth_lamLtbt = reshape(yhat.LO,sizeROI,length(yhat.LO)/sizeROI);
psth{cv}.yLtbt = mean(psth_lamLtbt,2)';

psth{cv}.yH = mean(psth.yHtbt);
psth{cv}.yL = mean(psth.yLtbt);

lam_NULL = mean(ytf);

%% Akaike   

LL_Htbt = tbt_yH.*log(psth_lamHtbt') - psth_lamHtbt' - log(factorial(tbt_yH));
LL_Ltbt = tbt_yL.*log(psth_lamLtbt') - psth_lamLtbt' - log(factorial(tbt_yL));

LL_nullH = tbt_yH.*log(lam_NULL) - lam_NULL - log(factorial(tbt_yH));
LL_nullL = tbt_yL.*log(lam_NULL) - lam_NULL - log(factorial(tbt_yL));

LLH = sum(LL_Htbt,1);
LLL = sum(LL_Ltbt,1);

LL.H.all = LL_Htbt;
LL.H.sum = LLH;

LL.L.all = LL_Ltbt;
LL.L.sum = LLL;

LL.FULL.all = vertcat(LL_Htbt,LL_Ltbt);
LL.FULL.sum = sum(LLH+LLL);

N = length(ytf);

k = size(X,2);

AIC_th = 2*k - 2.*LL.FULL.sum;
AIC = AIC_th + (2*k*(k+1))/(N-k-1);
    
end


