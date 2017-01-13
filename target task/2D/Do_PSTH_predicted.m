function varargout=Do_PSTH_predicted(pva,Trials,T,trials2,Spiketrains,neu, BINSIZE, PSTHLIM, Onset, lambda_pred, cv, eq_mag,nbins)
% % load('TrialsIdx','Trials','TrialTimesidx')
% load('MrT_data_9_24_2012')
% load('trials2')
%%
% aux=0.001:0.001:1;
pva.vel_ext=downsample(pva.vel,BINSIZE);
% PSTHLIM=[-2000,3000];
PSTHLIM(1)=-PSTHLIM(1);
% neu=14;
% figure,
spiketimes = 1000*Spiketrains;

spikes = histc(spiketimes, [1000*pva.vel(1,1):BINSIZE:1000*pva.vel(end,1)]);
% spikevec_PMd(:,neu) = spikes;
Spikes{neu}=spikes;

trial_inds = T.trial_inds;
trialnums = T.trialnums;
shuff_trials = T.shuff_trials;
CVranges=T.CVranges;

if cv == 0 
    testSet = 1:size(Trials.idx,1);
    uselam = mean(cell2mat(lambda_pred'),2);
else
    testSet = shuff_trials(CVranges(cv):CVranges(cv+1));
    uselam = lambda_pred{cv};
end


% forPSTH_spikes=zeros(size(Trials.idx,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
% forPSTH_spikes_pred=forPSTH_spikes;
% forPSTH_theta=zeros(size(Trials.idx,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
% forPSTH_speed=zeros(size(Trials.idx,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
% IdxU1=zeros(1,size(Trials.idx,1));
forPSTH_spikes=zeros(size(testSet,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
forPSTH_spikes_pred=zeros(size(testSet,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
forPSTH_theta=zeros(size(Trials.idx,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
forPSTH_speed=zeros(size(Trials.idx,1),(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1);
IdxU1=zeros(1,size(testSet,1));

% N_trials=size(Trials.idx,1);
% for trial=1:size(Trials.idx,1)-1
ind = 0;
for trial = testSet
    ind = ind+1;
    Times.Go=round(Trials.idx(trial,1)/BINSIZE);
    Times.Cloud=round(Trials.idx(trial,2)/BINSIZE);
    Times.End=round(Trials.idx(trial,3)/BINSIZE);
    if Onset==1
        Time.ONSET=Times.Go;
    elseif Onset==2
        Time.ONSET=Times.Cloud;
    elseif Onset==3
        Time.ONSET=Times.End;
    end
    
%     Times.Cloud=Times.Go;
    
    if trials2(trial,2)<2;
        IdxU1(ind)=1;
    else
        IdxU1(ind)=0;
    end
    forPSTH_spikes_pred(ind,:)=uselam(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE);
    forPSTH_spikes(ind,:)=Spikes{neu}(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE);
    forPSTH_theta(ind,:)=atan2(pva.vel_ext(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE,3),pva.vel_ext(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE,2));
    forPSTH_speed(ind,:)=sqrt(pva.vel_ext(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE,3).^2+pva.vel_ext(Time.ONSET-PSTHLIM(1)/BINSIZE:Time.ONSET+PSTHLIM(2)/BINSIZE,2).^2);
end


%%

edges=round(linspace(1,(PSTHLIM(1)+PSTHLIM(2))/BINSIZE+1,nbins+1));
for i=1:(length(edges)-1)
    N.u1(i)=length(find(IdxU1));
    aux=forPSTH_spikes(find(IdxU1),edges(i):edges(i+1))*1000/BINSIZE;
    Spike_mean.u1(i)=mean(aux(:));
    Spike_std.u1(i)=std(aux(:))/sqrt(length(aux(:)));
    
    N.u2(i)=length(find(~IdxU1));
    aux=forPSTH_spikes(find(~IdxU1),edges(i):edges(i+1))*1000/BINSIZE;
    Spike_mean.u2(i)=mean(aux(:));
    Spike_std.u2(i)=std(aux(:))/sqrt(length(aux(:)));
end

%%
%figure
%subplot(1,2,1)
%hold on
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u1,'k','LineWidth',2)
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u1+Spike_std.u1,'k','LineWidth',1)
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u1-Spike_std.u1,'k','LineWidth',1)
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u2,'r','LineWidth',2)
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u2+Spike_std.u2,'r','LineWidth',1)
plot(BINSIZE*(edges(1:end-1)-PSTHLIM(1)/BINSIZE+1+10/BINSIZE+1-1),Spike_mean.u2-Spike_std.u2,'r','LineWidth',1)

U1x = BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1);
%U1y = smooth(mean(forPSTH_spikes_pred(find(IdxU1==1),:))/BINSIZE,1);
U1y = mean(forPSTH_spikes_pred(find(IdxU1==1),:))/BINSIZE;

U2x = BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1);
%U2y = smooth(mean(forPSTH_spikes_pred(find(IdxU1==0),:))/BINSIZE,1);
U2y = mean(forPSTH_spikes_pred(find(IdxU1==0),:))/BINSIZE;

scaleu1u2 = max([max(Spike_mean.u1) max(Spike_mean.u2)])./ max([max(U1y) max(U2y)]);

if eq_mag == 1
    U1y = U1y * scaleu1u2;
    U2y = U2y * scaleu1u2;
end

plot(U1x,U1y,'-k','LineWidth',3,'Color',[55,192,252]/255)
plot(U2x,U2y,'-r','LineWidth',3,'Color',[252,189,0]/255)

% subplot(2,2,2)
% plot(BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1),((circ_mean(forPSTH_theta(find(IdxU1),:))))*180/pi,'-k','LineWidth',1)
% hold on
% plot(BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1),((circ_mean(forPSTH_theta(find(~IdxU1),:))))*180/pi,'-r','LineWidth',1)
% ylabel('angle')
% 
% subplot(2,2,4)
% plot(BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1),smooth(mean(forPSTH_speed(find(IdxU1),:)*180/pi,1)),'-k','LineWidth',1)
% hold on
% plot(BINSIZE*((1:(PSTHLIM(2)+PSTHLIM(1))/BINSIZE+1)-PSTHLIM(1)/BINSIZE+1),smooth(mean(forPSTH_speed(find(~IdxU1),:)*180/pi,1)),'-r','LineWidth',1)
% ylabel('speed')
