function plot_GLM(X,y,MODEL_1,MODEL_2,trials2,feed,min_inds,trialinds,aligntype,X2_p)

speedvec = X(:,3);

% Select alignment (feed_tbt for feedback or min_inds for minimum speed)
if aligntype == 1
    min_ind = feed;
else
    min_ind = min_inds;
end

for unit = 1:size(y,2)
    chun = MODEL_1{unit}.id;
    
    if X2_p{unit} < 0.05

    X_1 = [ones(length(X),1) X(:,MODEL_1{unit}.stats(2:end,1))];
    X_2 = [ones(length(X),1) X(:,MODEL_2{unit}.stats(2:end,1))];

    lambda_1 = exp(X_1*MODEL_1{unit}.stats(:,2));
    lambda_2 = exp(X_2*MODEL_2{unit}.stats(:,2));

    %Align by min_speed

    lam1_mat = nan(length(trials2),801);
    lam2_mat = nan(length(trials2),801);
    spike_mat = nan(length(trials2),801);
    speed_mat = nan(length(trials2),801);

    for i = 1:length(trials2)

        %fprintf('Trial: %d\n',i);

        current_bounds1 = lambda_1(trialinds==i);
        current_bounds2 = lambda_2(trialinds==i);
        current_boundspike = y(trialinds==i,unit);
        current_boundspeed = speedvec(trialinds==i);

        if (length(current_bounds1) - min_ind(i)) > 600
            lam1_mat(i,:) = current_bounds1(min_ind(i)-200:min_ind(i)+600)';

        else
            lam1_mat(i,1:(length(current_bounds1(min_ind(i)-200:end)))) = ...
                current_bounds1(min_ind(i)-200:end)';
        end


        if (length(current_bounds2) - min_ind(i)) > 600
            lam2_mat(i,:) = current_bounds2(min_ind(i)-200:min_ind(i)+600)';
        else 
            lam2_mat(i,1:(length(current_bounds2(min_ind(i)-200:end)))) = ...
                current_bounds2(min_ind(i)-200:end)';
        end

        if (length(current_boundspike) - min_ind(i)) > 600
            spike_mat(i,:) = current_boundspike(min_ind(i)-200:min_ind(i)+600)';
        else 
            spike_mat(i,1:(length(current_boundspike(min_ind(i)-200:end)))) = ...
                current_boundspike(min_ind(i)-200:end)';
        end

        if (length(current_boundspeed) - min_ind(i)) > 600
            speed_mat(i,:) = current_boundspeed(min_ind(i)-200:min_ind(i)+600)';
        else 
            speed_mat(i,1:(length(current_boundspeed(min_ind(i)-200:end)))) = ...
                current_boundspeed(min_ind(i)-200:end)';
        end

    end

    High_lam1_mat = lam1_mat(trials2(:,2) == max(trials2(:,2)),:);
    Low_lam1_mat = lam1_mat(trials2(:,2) == min(trials2(:,2)),:);

    High_lam2_mat = lam2_mat(trials2(:,2) == max(trials2(:,2)),:);
    Low_lam2_mat = lam2_mat(trials2(:,2) == min(trials2(:,2)),:);

    High_spike_mat = spike_mat(trials2(:,2) == max(trials2(:,2)),:);
    Low_spike_mat = spike_mat(trials2(:,2) == min(trials2(:,2)),:);

    High_speed_mat = speed_mat(trials2(:,2) == max(trials2(:,2)),:);
    Low_speed_mat = speed_mat(trials2(:,2) == min(trials2(:,2)),:);

    binstep = 25; 
    bincent = binstep:binstep:800 - binstep;
    bin_HSM = zeros(size(High_spike_mat,1),length(bincent));
    bin_LSM = zeros(size(Low_spike_mat,1),length(bincent));
    for i = 1:length(bincent)
        bin_HSM(:,i) = nanmean(High_spike_mat(:,bincent(i)-binstep+1:bincent(i)+binstep),2);
        bin_LSM(:,i) = nanmean(Low_spike_mat(:,bincent(i)-binstep+1:bincent(i)+binstep),2);
    end

    av_H1M = nanmean(High_lam1_mat,1);
    av_L1M = nanmean(Low_lam1_mat,1);

    av_H2M = nanmean(High_lam2_mat,1);
    av_L2M = nanmean(Low_lam2_mat,1);

    av_HSM = nanmean(bin_HSM,1);
    av_LSM = nanmean(bin_LSM,1);

    av_Hspeed = nanmean(High_speed_mat,1);
    av_Lspeed = nanmean(Low_speed_mat,1);

    figure; subplot(2,1,1);hold on;
    plot(-200:600,1000*av_H1M,'r--','LineWidth',1.5); 
    plot(-200:600,1000*av_L1M,'r','LineWidth',1.5);

    plot(-200:600,1000*av_H2M,'b--','LineWidth',1.5);
    plot(-200:600,1000*av_L2M,'b','LineWidth',1.5);

    plot(bincent - 200,1000*av_HSM,'k--');
    plot(bincent - 200,1000*av_LSM,'k');

    box off;
    %axis([0 500 15 45]);

    if aligntype == 1
        xlabel('Time from Feedback (ms)','FontSize',16);
    else
        xlabel('Time from Min Speed (ms)','FontSize',16);
    end
    ylabel('Firing Rate (spikes/sec)','FontSize',16);

    legend('Model 1','','Model 2','','Raw Data','');
    title(sprintf('ID %d (Channel: %d  Unit: %d)',...
        unit,floor(chun),round((chun-floor(chun))*10)),'FontSize',16);

    subplot(2,1,2); hold on;
    hold on; plot(-200:600, av_Hspeed, 'k--','LineWidth',2); 
    box off;
    axis([0 500 min(av_Hspeed)-10 max(av_Hspeed)+10]);
    plot(-200:600, av_Lspeed,'k','LineWidth',2);
    title('Speed Traces','FontSize',16);

    if aligntype == 1
        xlabel('Time from Feedback (ms)','FontSize',14);
    else
        xlabel('Time from Min Speed (ms)','FontSize',14);
    end

    ylabel('Hand Speed (cm/s)','FontSize',14);
    %legend('High Uncertainty','Low Uncertainty');

    end
end
% 