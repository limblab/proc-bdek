function [b,dev,stats,lambda_pred,LL,AIC]=GLM_fit_circ(Trials,X,spikes,Xnames,BASIS,T,CROSSNUM,BINSIZE,Range,FitOnRange)
    
%% 
acBasis = BASIS{1};
cBasis = BASIS{2};
ckinBasis = BASIS{3};
aceventBasis = BASIS{4};
ceventBasis = BASIS{5};

%%
clear aux R2Dt R2Dr
[x,y]=meshgrid(-99:1:100); y=-y;
ang=atan2(y,x);
rad2=sqrt(x.^2+y.^2);
min_rad=3;
max_rad=60;

BINZ=(rad2>min_rad).*(rad2<max_rad);
SB{1}=cos(ang).*BINZ;
SB{2}=sin(ang).*BINZ;


tempbasis.ac=acBasis;
tempbasis.c=cBasis;
tempbasis.ckin = ckinBasis;
tempbasis.acevent=aceventBasis;
tempbasis.cevent=ceventBasis;


nbasis.temp.ac=size(acBasis,2);
nbasis.temp.c=size(cBasis,2);
nbasis.temp.ckin = size(ckinBasis,2);
nbasis.temp.acevent = size(aceventBasis,2);
nbasis.temp.cevent = size(ceventBasis,2);

nbasis.sp.ac=length(Xnames.ac);
nbasis.sp.c=length(Xnames.c);
nbasis.sp.ckin = length(Xnames.ckin);
nbasis.sp.acevent = length(Xnames.acevent);
nbasis.sp.cevent = length(Xnames.cevent);

nbs = cell2mat(struct2cell(nbasis.sp));
nbt = cell2mat(struct2cell(nbasis.temp));
nbt(nbs==0)=0;
nbasis.temp = cell2struct(num2cell(nbt),{'ac','c','ckin','acevent','cevent'},1);

opts.dt=0.001;
% takethis= floor(Trials.idx(end,end)/BINSIZE);
% takethis= size(X,1);
% neu=14
% y=spikevec_PMd(:,neu);
y=spikes;
yc=(train2cont(y,50/BINSIZE))./BINSIZE + eps;%./(BINSIZE.^2);

% X=X(1:takethis,14);

% trialnums = 1:length(Trials.idx);
% shuff_trials = randperm(length(Trials.idx));
% CVranges=ceil(linspace(1,length(Trials.idx),CROSSNUM+1));
trial_inds = T.trial_inds;
trialnums = T.trialnums;
shuff_trials = T.shuff_trials;
CVranges=T.CVranges;
hinds = T.highinds;
linds = T.lowinds;

% CVranges=ceil(linspace(1,takethis,11));
% Shuf_trials=randperm(n_trials);
%fprintf('GLM for ...\n')
% CVranges
%fprintf('10-fold Cross-Validation...')

testpred = cell(CROSSNUM,1);
for cv=1:CROSSNUM%(length(CVranges)-1)
%     fprintf(num2str(cv));
%     if cv == CROSSNUM
%         fprintf('\n');
%     end

    testSet = shuff_trials(CVranges(cv):CVranges(cv+1));
    if CROSSNUM==1
        tranSet=testSet';
    else
        tranSet=setdiff(trialnums,testSet)';
    end
    
    feedbackindst = find(ismember(trial_inds,testSet));
    feedbackindsr = find(ismember(trial_inds,tranSet));
    
    feedbackH = find(ismember(trial_inds,hinds));
    feedbackL = find(ismember(trial_inds,linds));
    
    STARTROI = round(Range(1)/BINSIZE);
    ENDROI = round(Range(2)/BINSIZE);
    sizeROI = ENDROI-STARTROI+1;
    
    regionOI = zeros(length(Trials.idx),sizeROI);
    for i = 1:length(Trials.idx)
        regionOI(i,:) = round(Trials.idx(i,2)/BINSIZE+STARTROI):round(Trials.idx(i,2)/BINSIZE+ENDROI);
    end
    ALLregionOI = regionOI'; ALLregionOI = ALLregionOI(:);
    
    regionOIH = zeros(length(hinds),sizeROI);
    for i = 1:length(hinds)
        regionOIH(i,:) = round(Trials.idx(hinds(i),2)/BINSIZE+STARTROI):round(Trials.idx(hinds(i),2)/BINSIZE+ENDROI);
        tbt_yH(i,:) = y(regionOIH(i,:));
    end
    HregionOI = regionOIH'; HregionOI = HregionOI(:);
    psth_yH = mean(tbt_yH,1);
    
    regionOIL = zeros(length(linds),sizeROI);
    for i = 1:length(linds)
        regionOIL(i,:) = round(Trials.idx(linds(i),2)/BINSIZE+STARTROI):round(Trials.idx(linds(i),2)/BINSIZE+ENDROI);
        tbt_yL(i,:) = y(regionOIL(i,:));
    end
    LregionOI = regionOIL'; LregionOI = LregionOI(:);
    psth_yL = mean(tbt_yL,1);
    
    ROITinds = feedbackindst(ismember(feedbackindst,ALLregionOI));
    ROIRinds = feedbackindsr(ismember(feedbackindsr,ALLregionOI));
    ROI_H = feedbackH(ismember(feedbackH,HregionOI));
    ROI_L = feedbackL(ismember(feedbackL,LregionOI));
    
    Xt = X(ismember(trial_inds,testSet),:);
    Xr = X(ismember(trial_inds,tranSet),:);
    Xtr = X(ismember(trial_inds,unique([testSet tranSet'])),:);
    Xtf = X(ROITinds,:);
    XH = X(ROI_H,:);
    XL = X(ROI_L,:);
    
    yt = y(ismember(trial_inds,testSet),:);
    yr = y(ismember(trial_inds,tranSet),:);
    ytf = y(ROITinds,:);
    yH = y(ROI_H,:);
    yL = y(ROI_L,:);
    
    ytc = yc(ismember(trial_inds,testSet),:);
    yrc = yc(ismember(trial_inds,tranSet),:);
    ytfc = yc(ROITinds,:);
    ycH = yc(ROI_H,:);
    ycL = yc(ROI_L,:);
    
    if FitOnRange == 1
        
        Xt = X(ROITinds,:);
        Xr = X(ROIRinds,:);
        Xtr = X([ROITinds; ROIRinds],:);
        Xtf = X(ROITinds,:);
        
        yt = y(ROITinds,:);
        yr = y(ROIRinds,:);
    
        ytc = yc(ROITinds,:);
        yrc = yc(ROIRinds,:);
        
    end
    

%     Xt = X(testSet,:);
%     Xr = X(tranSet,:);
%     yt = y(testSet);
%     yr = y(tranSet);
    
%     yr=exp([ones(size(Xr,1),1) Xr]*bb_dummy_sal);
%     yt=exp([ones(size(Xt,1),1) Xt]*bb_dummy_sal);
    
    Z = Xr;
    Zt  = Xt;
    
    Z_save=Z;
%     B=ones(nbasis.temp.ac+nbasis.temp.c+1,1);

%% GLM ALL PARAMS

[b,dev,stats] = glmfit(Z,yr,'poisson');
bb_dummy_sal = b;

lambdar=exp([ones(size(Xr,1),1) Xr]*bb_dummy_sal)+eps;
lambdat=exp([ones(size(Xt,1),1) Xt]*bb_dummy_sal)+eps;
lambdatf = exp([ones(size(Xtf,1),1) Xtf]*bb_dummy_sal)+eps;
lambdaH = exp([ones(size(XH,1),1) XH]*bb_dummy_sal)+eps;
lambdaL = exp([ones(size(XL,1),1) XL]*bb_dummy_sal)+eps;

psth_lamHtbt = reshape(lambdaH,sizeROI,length(lambdaH)/sizeROI);
psth_lamH = mean(psth_lamHtbt,2)';

psth_lamLtbt = reshape(lambdaL,sizeROI,length(lambdaL)/sizeROI);
psth_lamL = mean(psth_lamLtbt,2)';

psth_lamHNULL = mean(psth_yH);
psth_lamLNULL = mean(psth_yL);

lam_NULL = mean(ytf);


lambda_pred{1} =(exp([ones(size(X,1),1) X]*bb_dummy_sal)+eps)./(BINSIZE*opts.dt);

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

% numk = struct2cell(Xnames(:));
% knames = horzcat(numk{:,:});
% k = length(knames);

k = size(X,2);

AIC_th = 2*k - 2.*LL.FULL.sum;
AIC = AIC_th + (2*k*(k+1))/(N-k-1);

    
end


