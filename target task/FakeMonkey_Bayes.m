%% Task Parameters
numslices = 10;
targetsize = 2;

mup = 2; % Prior mean
sigp = 3; % Prior std

% Feedback standard deviations
s1 = 1;%2.2; %1
s2 = 4;%7; %4

ntrials = 1000; %Trials per session
reps = 500; %Total number of sessions       
sigp_hats = 0:0.1:4; %Range of prior estimates (strategies)
strategies = [3]; %Strategies for which you want scatterplot with slopes

%% Mr. F (the fake monkey)
centroid_estimation_error = 0; % percent var
motor_error = 0; % std (constant)

%% Initialize arrays
correct1 = zeros(reps,length(sigp_hats));
correct2 = zeros(reps,length(sigp_hats));
totcorrect = zeros(reps,length(sigp_hats));

sse1 = zeros(reps,length(sigp_hats));
sse2 = zeros(reps,length(sigp_hats));
sset = zeros(reps,length(sigp_hats));

sumsig1 = zeros(reps,length(sigp_hats));
sumsig2 = zeros(reps,length(sigp_hats));
sumsigt = zeros(reps,length(sigp_hats));

theor = zeros(length(sigp_hats),2);

sigp_hats = 1/1000*round(1000*sigp_hats); %Range of prior estimates (strategies)
strategies = 1/1000*round(1000*strategies); %Strategies for which you want scatterplot with slopes
[throwaway,plotstrategies] = ismember(strategies,sigp_hats);
%% Do Simulation
for rep = 1:reps %Cycle through sessions

    % Obtain true target positions and separate into high/low
    targets = mup + randn(ntrials,1).*sigp;
    targetsL = targets(1:end/2);
    targetsH = targets(end/2+1:end);

    % Represent the targets with slices (drawn from feedback stds)
    slicesf1 = repmat(targetsL,1,numslices) + s1*randn(ntrials/2,numslices);
    slicesf2 = repmat(targetsH,1,numslices) + s2*randn(ntrials/2,numslices);

    % Calculate the centroid positions and stds of the slices
    mu_hats1 = mean(slicesf1,2);
    sig_hats1 = std(slicesf1,0,2)./sqrt(numslices);
    
    mu_hats1 = mu_hats1 + centroid_estimation_error.*sig_hats1.*randn(length(mu_hats1),1);

    mu_hats2 = mean(slicesf2,2);
    sig_hats2 = std(slicesf2,0,2)./sqrt(numslices);

    mu_hats2 = mu_hats2 + centroid_estimation_error.*sig_hats2.*randn(length(mu_hats2),1);
    
    if rep==1; figure; hold on; end
    
    for i = 1:length(sigp_hats) % Cycle through strategies
        sigp_hat = sigp_hats(i);

        % Bayesian integration to calculate final target position estimate
        pos1 = (mu_hats1.*sigp_hat.^2 + mup.*sig_hats1.^2)./(sigp_hat.^2+sig_hats1.^2);
         pos1 = pos1 + motor_error*randn(length(pos1),1);
        pos2 = (mu_hats2.*sigp_hat.^2 + mup.*sig_hats2.^2)./(sigp_hat.^2+sig_hats2.^2);
         pos2 = pos2 + motor_error*randn(length(pos2),1);
        if rep == 1 
            
            % Calculate Theoretical Slopes
            theor(i,1) = sigp_hat^2/(sigp_hat^2+(s1/sqrt(numslices))^2);
            theor(i,2) = sigp_hat^2/(sigp_hat^2+(s2/sqrt(numslices))^2);
            
            % Plot 2004 Nature plot for strategies of interest (see TASK
            % PARAMETERS)
            if ismember(i,plotstrategies)
                subplot(1,length(plotstrategies),find(plotstrategies == i)); hold on;
                plot(mu_hats1,pos1,'b.');
                plot(mu_hats2,pos2,'r.');

                [slope1 S1] = polyfit(mu_hats1,pos1,1); 
                [slope2 S2] = polyfit(mu_hats2,pos2,1);

                plot([-25 25],polyval(slope1,[-25 25]),'b');
                legend(sprintf('Slope: %.3f',slope1(1)),sprintf('Slope: %.3f',slope2(1)))
                plot([-25 25],polyval(slope2,[-25 25]),'r');
                title(sprintf('Prior Est: %.1f',sigp_hat));

                axis([min(targets) max(targets) min(targets) max(targets)]);
            end
        end
        
        % Calculate standard deviations of posteriors
        sig1 = sqrt((sigp_hat.^2*sig_hats1.^2)./(sigp_hat.^2+sig_hats1.^2));
        sig2 = sqrt((sigp_hat.^2*sig_hats2.^2)./(sigp_hat.^2+sig_hats2.^2));
        
        % Find average posterior std
        sumsig1(rep,i) = mean(sig1);
        sumsig2(rep,i) = mean(sig2);
        sumsigt(rep,i) = mean([sig1;sig2]);
        
        % Calculate errors of MAP estimate
        err1 = abs(pos1 - targetsL);
        err2 = abs(pos2 - targetsH);
        
        % SSE for separate cases and also combined
        sse1(rep,i) = sum(err1.^2);
        sse2(rep,i) = sum(err2.^2);
        sset(rep,i) = sse1(rep,i)+sse2(rep,i);
        
        % Calculate success rate (% targets hit)
        correct1(rep,i) = sum(err1<=targetsize/2)./(ntrials/2);
        correct2(rep,i) = sum(err2<=targetsize/2)./(ntrials/2);
        totcorrect(rep,i) = (correct1(rep,i)*(ntrials/2)+correct2(rep,i)*(ntrials/2))./ntrials;
    end
end

%% Calculate mean and stds of success rates
% Feedback 1
CP1 = mean(correct1,1);
cp1sd = std(correct1,0,1);
% Feedback 2
CP2 = mean(correct2,1);
cp2sd = std(correct2,0,1);
% Total
TC = mean(totcorrect,1);
tcsd = std(totcorrect,0,1);

%% Calculate mean and stds of SSE
% Feedback 1
SSE1 = mean(sse1,1);
SSE1std = std(sse1,0,1);
% Feedback 2
SSE2 = mean(sse2,1);
SSE2std = std(sse2,0,1);
% Total
SSET = mean(sset,1);
SSETstd = std(SSET,0,1);

%% Calculate mean and stds of Posterior Standard deviations
% Feedback 1
SIG1 = mean(sumsig1,1);
SIG1std = std(sumsig1,0,1);
% Feedback 2
SIG2 = mean(sumsig2,1);
SIG2std = std(sumsig2,0,1);
% Total
SIGT = mean(sumsigt,1);
SIGTstd = std(sumsigt,0,1);

%% Create Performance Figures
% Performance (Percent correct) across strategies
figure; 
subplot(3,2,1); hold on;
errorbar(sigp_hats,CP1,cp1sd,'bo'); title('Low Var Trials','FontSize',16);
plot([sigp sigp],[0 1],'g--')
subplot(3,2,3); hold on;
errorbar(sigp_hats,CP2,cp2sd,'ro'); title('High Var Trials','FontSize',16);
plot([sigp sigp],[0 1],'g--')
subplot(3,2,5); hold on;
errorbar(sigp_hats,TC,tcsd,'ko'); title('All Trials','FontSize',16);
plot([sigp sigp],[0 1],'g--')

CP1max = max(CP1);
CP2max = max(CP2);
TCmax = max(TC);

d1 = CP1max- CP1;
d2 = CP2max - CP2;
dt = TCmax - TC;

weights = 0.1:0.1:10;
for i = 1:length(weights)
    cost1(i,:) = weights(i)*d1 + SIG1;
        best1(i) = sigp_hats(cost1(i,:)==min(cost1(i,:)));
    cost2(i,:) = weights(i)*d2 + SIG2;
        best2(i) = sigp_hats(cost2(i,:)==min(cost2(i,:)));
    costt(i,:) = weights(i)*dt + SIGT;
        bestt(i) = sigp_hats(costt(i,:)==min(costt(i,:)));
end

xlabel('Estimated Prior SD','FontSize',14); 
ylabel('Performance (% correct)','FontSize',14);

% Performance (SSE) across strategies
subplot(3,2,2); hold on;
errorbar(sigp_hats,SSE1,SSE1std,'bo');
plot([sigp sigp],[0 500],'g--')
subplot(3,2,4); hold on;
errorbar(sigp_hats,SSE2,SSE2std,'ro'); 
plot([sigp sigp],[0 500],'g--')
subplot(3,2,6); hold on;
errorbar(sigp_hats,SSET,SSETstd,'ko');
plot([sigp sigp],[0 1000],'g--')

xlabel('Estimated Prior SD','FontSize',14); 
ylabel('Performance (SSE)','FontSize',14);

%% Create Robustness Figures
% Return the best Prior estimate (SHOULD BE THE TRUE PRIOR)
indbest = find(TC==max(TC));
fprintf('Optimal Prior: %.3f\n',sigp_hats(indbest));

% Perform ttest to compare performance of strategies wrt optimal prior
bestset = totcorrect(:,indbest);
comp = zeros(size(totcorrect,2));
pval = comp;
for i = 1:size(totcorrect,2)
    [comp(i) pval(i)] = ttest2(totcorrect(:,i),bestset);
end

% Plot P-values from comparison
figure; plot(sigp_hats,pval,'b.-');
xlabel('Prior Estimates','FontSize',14); 
ylabel('P value','FontSize',14);
title('Significance of difference from optimal');
%axis([sigp_hats(1) sigp_hats(end) 0 .06]);

% Plot Slopes across strategies
figure; hold on;
plot(sigp_hats,theor(:,1),'b');
plot(sigp_hats,theor(:,2),'r');
plot([sigp sigp],[0 1],'g--')
legend('Feedback 1','Feedback 2');
xlabel('Prior Estimates','FontSize',14);
ylabel('Theoretical Slopes','FontSize',14);
title('Prior Estimate Effect on Slope','FontSize',15);

dslope1 = theor(indbest,1) - theor(:,1);
dslope2 = theor(indbest,2) - theor(:,2);

drew1 = CP1(indbest) -CP1;
drew2 = CP2(indbest) -CP2;


