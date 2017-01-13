function[Recruit] = PPC_looper(alldays,PDI,TA,good_neurs,firing_perc,like_ind) 

%percnegfunc = @(F,mint,maxt) 2*F./(maxt-mint)-(maxt+mint)./(maxt-mint);
%gainfunc = @(TA,act,gn,dt) (act(gn)'/(dt/1000) - min(TA(gn,:),[],2))./(max(TA(gn,:),[],2)-min(TA(gn,:),[],2));

for i = 1:length(firing_perc{1})
    clc; fprintf('time bin: %d/%d\n',i,length(firing_perc{1}));
    nTA = (TA{end}(good_neurs{i},:) - repmat(min(TA{end}(good_neurs{i},:),[],2),1,628))./...
        (repmat(max(TA{end}(good_neurs{i},:),[],2) - min(TA{end}(good_neurs{i},:),[],2),1,628));
% 
%     nTA = (TA{end}(good_neurs,:) - repmat(min(TA{end}(good_neurs,:),[],2),1,628))./...
%         (repmat(max(TA{end}(good_neurs,:),[],2) - min(TA{end}(good_neurs,:),[],2),1,628));
    
%     nTA = (TA(good_neurs,:) - repmat(min(TA(good_neurs,:),[],2),1,628))./...
%         (repmat(max(TA(good_neurs,:),[],2) - min(TA(good_neurs,:),[],2),1,628));
    
    base = nansum(nTA);
    %% Do center out simulations to get estimated gains for all neurons
%     co_list = 1:628;
%     co_num = length(co_list);
% 
%     ec = good_neurs(:,co_list)';
%     %new_moduls = repmat(MODS{i}(1,:),co_num,1);
% %     e_gain = (ec(:,good_neurs) - ...
% %                 repmat(min(TA{i}(good_neurs,:),[],2)',co_num,1))./new_moduls(:,good_neurs);
% 
%     e_gain = percnegfunc(ec(:,good_neurs), repmat(min(good_neurs(good_neurs,:),[],2)',co_num,1),...
%         repmat(max(good_neurs(good_neurs,:),[],2)',co_num,1));
%                 
%             
%     ALLAF = zeros(co_num,628);
%     for loc = 1:co_num
% 
%         gains = e_gain(loc,:);
%         efunc = sum(repmat(gains',1,628).*nTA);
% 
%         move_ind = loc;
% 
%         atzero = [efunc(move_ind:end) efunc(1:(move_ind-1))];
%         ALLAF(loc,:) = [atzero(end/2:end) atzero(1:(end/2)-1)];
%     end
%     base_recruit = mean(ALLAF);
    %% Get UN recruitments
    [AF] = deal(zeros(size(alldays(PDI).tt,1),628));
    for trial = 1:size(alldays(PDI).tt,1)

        %gains = firing_diffs{1}{1}(trial,:);
        gains = firing_perc{1}{i}(trial,:); %gains(gains > 1) = 1; gains(gains < 0) = 0;
        
        %gains = gainfunc(TA,raw_spikes{1}{i}(trial,:),good_neurs,binsize_ms);
        
        %egains = efiring_perc{1}{1}(trial,:);

        func = nansum(repmat(gains',1,628).*nTA)./base; 
        %efunc = sum(repmat(egains',1,628).*nTA)./base;

        %add_func = func./efunc;
        add_func = func;%./base;

        %moveloc = alldays(PDI).tt(trial,10);
        moveloc = mean(alldays(PDI).tt(:,10));
        move_ind = round(100*moveloc);


        if move_ind == 0 
            atzero = add_func; 
            %eatzero = e
        else
            atzero = [add_func(move_ind:end) add_func(1:(move_ind-1))];
        end
        
        %EF(trial,:) = [
        AF(trial,:) = [atzero(end/2:end) atzero(1:(end/2)-1)];
    end

    % Clean up

    [bL,bU,av_rec,REC] = deal(cell(length(like_ind{1}),1));
    for lik = 1:length(like_ind{1})

        [bL{lik}, bU{lik}] = boot_bounds(1000,@nanmean,AF(like_ind{1}{lik},:),2.5,97.5);

        REC{lik} = AF(like_ind{1}{lik},:);
        av_rec{lik} = nanmean(REC{lik})';
    end
        
    Recruit.low(:,i) = bL;
    Recruit.high(:,i) = bU; 
    Recruit.av(:,i) = av_rec;
    %Recruit.baseline(:,i) = base_recruit';
    Recruit.all{i} = REC;
end

 