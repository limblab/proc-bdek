%%
figure; hold on;
for i = low_devind2'
    
    bdfstart = find(alldays(2).bdfM.pos(:,1) > tt2(i,6),1,'first');
    bdfend = find(alldays(2).bdfM.pos(:,1) < tt2(i,7),1,'last');
    
    plot(alldays(2).bdfM.pos(bdfstart:bdfend,2),alldays(2).bdfM.pos(bdfstart:bdfend,3),'b');
end

for i = high_devind2'
    
    bdfstart = find(alldays(2).bdfM.pos(:,1) > tt2(i,6),1,'first');
    bdfend = find(alldays(2).bdfM.pos(:,1) < tt2(i,7),1,'last');
    
    plot(alldays(2).bdfM.pos(bdfstart:bdfend,2),alldays(2).bdfM.pos(bdfstart:bdfend,3),'r');
end

%%
figure; 
hold on;
errorbar(-300:200:700,nanmean(ep2(:,feedind2),2),1.96*nanstd(ep2(:,feedind2),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ep2(:,priorind2),2),1.96*nanstd(ep2(:,priorind2),[],2)/10,'b.-');

figure; hold on;
errorbar(-300:200:700,nanmean(ep2(:,high_devind2),2),1.96*nanstd(ep2(:,high_devind2),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ep2(:,low_devind2),2),1.96*nanstd(ep2(:,low_devind2),[],2)/10,'b.-');

figure; hold on;
errorbar(-300:200:700,nanmean(ep3(:,high_devind3),2),1.96*nanstd(ep3(:,high_devind3),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ep3(:,low_devind3),2),1.96*nanstd(ep3(:,low_devind3),[],2)/10,'b.-');

%%
figure; 
hold on;
errorbar(-300:200:700,nanmean(el2(:,high_devind2),2),1.96*nanstd(el2(:,high_devind2),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(tl2(:,high_devind2),2),1.96*nanstd(tl2(:,high_devind2),[],2)/10,'b.-');

figure; 
hold on;
errorbar(-300:200:700,nanmean(ve3(:,high_devind3),2),1.96*nanstd(ve3(:,high_devind3),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ve3(:,low_devind3),2),1.96*nanstd(ve3(:,low_devind3),[],2)/10,'b.-');

%%
figure; 
hold on;
errorbar(-300:200:700,nanmean(ep2(:,priorleft),2),1.96*nanstd(ep2(:,priorleft),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ep2(:,priorright),2),1.96*nanstd(ep2(:,priorright),[],2)/10,'k.-');

%%
figure; hold on;
errorbar(-300:200:700,nanmean(ve2(:,prior_driven),2),1*nanstd(ve2(:,prior_driven),[],2)/10,'r.-');
errorbar(-300:200:700,nanmean(ve2(:,targ_driven),2),1*nanstd(ve2(:,targ_driven),[],2)/10,'b.-');

%% Day2 Target likelihood
move2prior = find(abs(tt2(:,10)-pi/2) < 2.5/180*pi);
lindsp2 = move2prior(tt2(move2prior,3)==max(tt2(move2prior,3)));
hindsp2 = move2prior(tt2(move2prior,3)==min(tt2(move2prior,3)));

[ig, tfromp] = sortrows(abs(tt2(move2prior,9)-pi/2));

prior_driven = move2prior(tfromp(end/2+1:end));
targ_driven = move2prior(tfromp(1:end/2));

figure; hold on;
errorbar(-300:200:700,nanmean(tl2(:,prior_driven),2),1*nanstd(tl2(:,prior_driven),[],2)/sqrt(length(prior_driven)),'r.-');
errorbar(-300:200:700,nanmean(tl2(:,targ_driven),2),1*nanstd(tl2(:,targ_driven),[],2)/sqrt(length(targ_driven)),'b.-');

figure; hold on;
errorbar(-300:200:700,nanmean(el2(:,prior_driven),2),1*nanstd(el2(:,prior_driven),[],2)/sqrt(length(prior_driven)),'r.-');
errorbar(-300:200:700,nanmean(el2(:,targ_driven),2),1*nanstd(el2(:,targ_driven),[],2)/sqrt(length(targ_driven)),'b.-');

%% Day3 Target Likelihood
move2prior3 = find(abs(tt3(:,10)-pi/2) < 5/180*pi);
lindsp3 = move2prior3(tt3(move2prior3,3)==max(tt3(move2prior3,3)));
hindsp3 = move2prior3(tt3(move2prior3,3)==min(tt3(move2prior3,3)));

[ig3, tfromp3] = sortrows(abs(tt3(move2prior3,9)-pi/2));

prior_driven3 = move2prior3(tfromp3(end/2+1:end));
targ_driven3 = move2prior3(tfromp3(1:end/2));

figure; hold on;
errorbar(-300:200:700,nanmean(tl3(:,prior_driven3),2),1*nanstd(tl3(:,prior_driven3),[],2)/sqrt(length(prior_driven3)),'r.-');
errorbar(-300:200:700,nanmean(tl3(:,targ_driven3),2),1*nanstd(tl3(:,targ_driven3),[],2)/sqrt(length(targ_driven3)),'b.-');

figure; hold on;
errorbar(-300:200:700,nanmean(el3(:,prior_driven3),2),1*nanstd(el3(:,prior_driven3),[],2)/sqrt(length(prior_driven3)),'r.-');
errorbar(-300:200:700,nanmean(el3(:,targ_driven3),2),1*nanstd(el3(:,targ_driven3),[],2)/sqrt(length(targ_driven3)),'b.-');


%%
figure; hold on;
errorbar(-300:200:700,nanmean(targ_or_end(:,high_devind2),2),...
    1.96*nanstd(targ_or_end(:,high_devind2),[],2)/sqrt(length(high_devind2)),'b.-')
errorbar(-300:200:700,nanmean(targ_or_end(:,low_devind2),2),...
    1.96*nanstd(targ_or_end(:,low_devind2),[],2)/sqrt(length(low_devind2)),'r.-')

figure; hold on;
errorbar(-300:200:700,nanmean(targ_or_end(:,prior_driven),2),...
    1.96*nanstd(targ_or_end(:,prior_driven),[],2)/sqrt(length(prior_driven)),'b.-')
errorbar(-300:200:700,nanmean(targ_or_end(:,targ_driven),2),...
    1.96*nanstd(targ_or_end(:,targ_driven),[],2)/sqrt(length(targ_driven)),'r.-')


