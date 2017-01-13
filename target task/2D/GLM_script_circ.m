%% Options
BINSIZE = 10;
Brain_Area='PMd';
neu=13;

%% Set up BASIS 
[iht, ihbas, ihbasis] = makeRaisedCosBasis(5, BINSIZE/1000, [2/1000 25/1000]*7, 1*1.5*1e-2, 1);
Basis = ihbasis;

BASIS{1} = Basis; %ac
BASIS{2} = Basis; %c
BASIS{3} = Basis; %ckin
BASIS{4} = Basis; %acevent
BASIS{5} = Basis; %cevent

%% choose neuron
if  strcmp(Brain_Area,'PMd')
    Spiketrains=PMd_units{neu};
elseif strcmp(Brain_Area,'M1')
    Spiketrains=M1_units{neu};
else
    fprintf('choose between "PMd" or "M1"')
end

%% Build Covariates
[X1,spikes1,Xnames1,trial_inds]=Do_Covariates_circ(BDF,Spiketrains,comp_tt,BINSIZE,1,BASIS);
[X2,spikes2,Xnames2,trial_inds2]=Do_Covariates_circ(BDF,Spiketrains,comp_tt,BINSIZE,2,BASIS);

%% GLM (GLM_fit)
DoPlot=0;
CROSSNUM = 1;
Range = [0 500]; %Range from Feedback
fit_on_range = 0; %0 for evaluation only. 1 for Fit and Eval
 
shuffleTrials_circ(Trials,comp_tt,trial_inds,CROSSNUM)

[b1,dev1,stats1,lambda_pred1,LL1,AIC1]=GLM_fit_circ(Trials,X1,spikes1,Xnames1,BASIS,T,CROSSNUM,BINSIZE,Range,fit_on_range);
[b2,dev2,stats2,lambda_pred2,LL2,AIC2]=GLM_fit_circ(Trials,X2,spikes2,Xnames2,BASIS,T,CROSSNUM,BINSIZE,Range,fit_on_range);

% Get CHI^2 
X2_p = 1 - chi2cdf(dev1 - dev2,size(X2,2)-size(X1,2));

% AIC
AIC_diff = exp((AIC1-AIC2)./2);

%fprintf('Model with uncertainty is %.2f times more likely (AIC)\nX2 p-value: %.3f\n',...
%    AIC_diff,X2_p);

lambda_pred1{1} = lambda_pred1{1}*10;
lambda_pred2{1} = lambda_pred2{1}*10;

%% plot predicted PSTHs
Onset=2;

PSTHLIM=[0,500];
nbins = 20;
%nbins = diff(PSTHLIM)/;

cv = 0;
figure;
subplot(1,2,1); hold on; cla; title(sprintf('%s:%d - No Uncertainty',Brain_Area,neu));
Do_PSTH_predicted(pva,Trials,T,trials2,Spiketrains,neu, BINSIZE, PSTHLIM, Onset, lambda_pred1, cv, 0,nbins);
subplot(1,2,2); hold on; cla; title(sprintf('%s:%d - Uncertainty',Brain_Area,neu));
Do_PSTH_predicted(pva,Trials,T,trials2,Spiketrains,neu, BINSIZE, PSTHLIM, Onset, lambda_pred2, cv, 0,nbins); 


