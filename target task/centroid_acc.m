numslices = 10;

comp_tt = tt(tt(:,end)==32 | tt(:,end)==34,:);

bursttimes = vertcat(BDF.databursts{:,1});
slice_x = zeros(length(comp_tt),numslices);

centroid = zeros(length(comp_tt),1);
estimate = zeros(length(comp_tt),1);
targetp = zeros(length(comp_tt),1);
for i = 1:length(comp_tt)
    burst = BDF.databursts{comp_tt(i,1)==bursttimes,2};
    xypos = bytes2float(burst(end-(4*numslices*2-1):end));
    slice_x(i,:) = xypos(1:2:end);
    
    endtime = comp_tt(i,7);
    endind = find(BDF.pos(:,1)>endtime,1,'first');
    
    centroid(i) = mean(slice_x(i,:));
    estimate(i) = BDF.pos(endind,2);
    targetp(i) = comp_tt(i,2);
    
end

err.act = estimate-targetp;
err.the = centroid-targetp;

badtrials = find(abs(err.the)>100);

err.act(badtrials) = [];
err.the(badtrials) = [];
comp_tt(badtrials,:) = [];
targetp(badtrials) = [];
centroid(badtrials) = [];
estimate(badtrials) = [];

ind1 = find(comp_tt(:,3)==1);
ind2 = find(comp_tt(:,3)==2);
ind3 = find(comp_tt(:,3)==3);
ind4 = find(comp_tt(:,3)==4);

mean_err.act = [mean(err.act(ind1));mean(err.act(ind2));mean(err.act(ind3));mean(err.act(ind4))];
std_err.act = [std(err.act(ind1));std(err.act(ind2));std(err.act(ind3));std(err.act(ind4))];

mean_err.the = [mean(err.the(ind1));mean(err.the(ind2));mean(err.the(ind3));mean(err.the(ind4))];
std_err.the = [std(err.the(ind1));std(err.the(ind2));std(err.the(ind3));std(err.the(ind4))];

figure; hold on; 
subplot(1,2,1); plot(mean_err.the,mean_err.act,'b.'); 
xlabel('Theoretical Mean Error','FontSize',14);
ylabel('Actual Mean Error','FontSize',14);

subplot(1,2,2); plot(std_err.the,std_err.act,'r.-');
hold on; plot([std_err.the(1) std_err.the(end)],[std_err.the(1) std_err.the(end)],'k--');
hold on; plot([1 4]./sqrt(10),[1 4]./sqrt(10),'g--');
xlabel('Theoretical SEM','FontSize',14);
ylabel('Actual SEM','FontSize',14);

%% Nature
figure; hold on;
subplot(1,4,1); hold on;
plot(centroid(ind1),estimate(ind1),'b.');
[slope.f1,SL.f1] = polyfit(centroid(ind1),estimate(ind1),1);
[YL.f1,DELTAL.f1] = polyconf(slope.f1,centroid(ind1),SL.f1,'predopt','curve');
plot(centroid(ind1),YL.f1,'b','LineWidth',2);
title('SD 1','FontSize',14);
legend(sprintf('Fit Slope = %.3f',slope.f1(1)));

subplot(1,4,2); hold on;
plot(centroid(ind2),estimate(ind2),'r.');
[slope.f2,SL.f2] = polyfit(centroid(ind2),estimate(ind2),1);
[YL.f2,DELTAL.f2] = polyconf(slope.f2,centroid(ind2),SL.f2,'predopt','curve');
plot(centroid(ind2),YL.f2,'r','LineWidth',2);
title('SD 2','FontSize',14);
legend(sprintf('Fit Slope = %.3f',slope.f2(1)));

subplot(1,4,3); hold on;
plot(centroid(ind3),estimate(ind3),'g.');
[slope.f3,SL.f3] = polyfit(centroid(ind3),estimate(ind3),1);
[YL.f3,DELTAL.f3] = polyconf(slope.f3,centroid(ind3),SL.f3,'predopt','curve');
plot(centroid(ind3),YL.f3,'g','LineWidth',2);
title('SD 3','FontSize',14);
legend(sprintf('Fit Slope = %.3f',slope.f3(1)));

subplot(1,4,4); hold on;
plot(centroid(ind4),estimate(ind4),'k.');
[slope.f4,SL.f4] = polyfit(centroid(ind4),estimate(ind4),1);
[YL.f4,DELTAL.f4] = polyconf(slope.f4,centroid(ind4),SL.f4,'predopt','curve');
plot(centroid(ind4),YL.f4,'k','LineWidth',2);
title('SD 4','FontSize',14);
legend(sprintf('Fit Slope = %.3f',slope.f4(1)));




