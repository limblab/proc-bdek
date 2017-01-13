%Activity_script2;

figure; hold on;
AG = cell(1,length(TA));
for bin = 1:length(TA)

    nTA = (TA{bin}(good_neurs{bin},:) - repmat(min(TA{bin}(good_neurs{bin},:),[],2),1,628))./...
    (repmat(max(TA{bin}(good_neurs{bin},:),[],2) - min(TA{bin}(good_neurs{bin},:),[],2),1,628));

    base = sum(nTA);
%% Do center out simulations to get estimated gains for all neurons
co_list = 1:628;
co_num = length(co_list);

ec = TA{bin}(:,co_list)';
new_moduls = repmat(MODS{bin}(1,:),co_num,1);
e_gain = (ec(:,good_neurs{bin}) - ...
            repmat(min(TA{bin}(good_neurs{bin},:),[],2)',co_num,1))./new_moduls(:,good_neurs{bin});

ALLAF = zeros(co_num,628);
for loc = 1:co_num

    gains = e_gain(loc,:);
    efunc = sum(repmat(gains',1,628).*nTA);
    
    moveloc = wrapped_cents(loc);
    move_ind = loc;
    
    atzero = [efunc(move_ind:end) efunc(1:(move_ind-1))];
    ALLAF(loc,:) = [atzero(end/2:end) atzero(1:(end/2)-1)];
end
base_recruit = mean(ALLAF);
%% Get UN recruitments
[AF, COAF, UNAF] = deal(zeros(size(alldays(2).tt,1),628));
[AFs, COAFs, UNAFs] = deal(zeros(size(alldays(2).tt,1),314));
baseline_func = sum(nTA);
slice_atzero = zeros(size(alldays(2).tt,1),5);
for trial = 1:size(alldays(2).tt,1)
    
    %gains = firing_diffs{1}{1}(trial,:);
    gains = firing_perc{1}{bin}(trial,:); %gains(gains > 1) = 1; gains(gains < 0) = 0;
    egains = efiring_perc{1}{bin}(trial,:);
    
    func = sum(repmat(gains',1,628).*nTA)./base; 
    efunc = sum(repmat(egains',1,628).*nTA)./base;

    un_func = func;
    co_func = efunc;
    add_func = func./efunc;

    moveloc = alldays(2).tt(trial,10);
    move_ind = round(100*moveloc);
    
%     delfun = zeros(1,628); delfun(move_ind) = 1000;
%     delfun = circ_vmpdf(linspace(0,2*pi,628),moveloc,100)';
%     tranfun = fft(delfun)./fft(efunc);
%     
%     add_func = ifft(fft(func).*tranfun);
    if move_ind == 0 
        atzero = add_func; 
        co_atzero = co_func;
        un_atzero = un_func;
    else
        atzero = [add_func(move_ind:end) add_func(1:(move_ind-1))];
        co_atzero = [co_func(move_ind:end) co_func(1:(move_ind-1))];  
        un_atzero = [un_func(move_ind:end) un_func(1:(move_ind-1))];
    end
    
    AF(trial,:) = [atzero(end/2:end) atzero(1:(end/2)-1)];
    AFs(trial,:) = nanmean([atzero(315:628); atzero(314:-1:1)], 1);
    
    COAF(trial,:) = [co_atzero(end/2:end) co_atzero(1:(end/2)-1)];
    COAFs(trial,:) = nanmean([co_atzero(315:628); co_atzero(314:-1:1)],1);
    
    UNAF(trial,:) = [un_atzero(end/2:end) un_atzero(1:(end/2)-1)];
    UNAFs(trial,:) = nanmean([un_atzero(315:628); un_atzero(314:-1:1)],1);
    
    slice_atzero(trial,:) = alldays(2).slices(trial,:) - moveloc;
end

%% Plotting
G = [AFs fliplr(AFs)];

[bL,bU,av_gain,co_recruit, un_recruit, SAZ,slice_notinco,histed, GAIN] = ...
    deal(cell(length(set_of_inds),1));

sthet = linspace(-180,180,628);
sthetD = [linspace(-180,180,628), linspace(180,-180,628)];

cop = {'b','g','r'};
subplot(1,length(TA),bin); hold on;

X = G;
for lik = 1:length(set_of_inds)
    
    SAZ{lik} = slice_atzero(like_ind{1}{lik},:); SAZ{lik} = SAZ{lik}(:)*180/pi;
    slice_notinco{lik} = SAZ{lik}(abs(SAZ{lik}) > 15);
    histed{lik} = histc(SAZ{lik},-180:10:180)./length(like_ind{1}{lik});
    %histed{lik} = histc(slice_notinco{lik},-180:10:180)./length(like_ind{1}{lik});
    
    [bL{lik}, bU{lik}] = boot_bounds(1000,@nanmean,X(like_ind{1}{lik},:),2.5,97.5);
%     [bL{lik}, bU{lik}] = boot_bounds(1000,@nanmean,AF(like_ind{1}{lik},:),2.5,97.5);
    
%     offsetter = mean(mean(X(like_ind{1}{lik},:)));
%     offsetter = mean(X(like_ind{1}{lik},314));
%     offsetter = min(mean(X(like_ind{1}{lik},:)));
    offsetter = 0;
    
%     X(like_ind{1}{lik},:) = X(like_ind{1}{lik},:);
    co_recruit{lik} = mean(COAF(like_ind{1}{lik},:));
    un_recruit{lik} = mean(UNAF(like_ind{1}{lik},:));
    av_gain{lik} = (mean(X(like_ind{1}{lik},:)) - offsetter);
    bL{lik} = (bL{lik} - offsetter);
    bU{lik} = (bU{lik} - offsetter);
    
    GAIN{lik} = AF(like_ind{1}{lik},:);
    
%     
    av_gain{lik} = av_gain{lik}.*base_recruit;
    bL{lik} = bL{lik}.*base_recruit';
    bU{lik} = bU{lik}.*base_recruit';
    
    plot(sthet,av_gain{lik},cop{lik});
    patch(sthetD,[bL{lik}' fliplr(bU{lik}')],cop{lik},'FaceAlpha',0.25,'EdgeAlpha',0.25);
    
    AG{bin}{lik} = [av_gain{lik}; bL{lik}' ;bU{lik}'];
    
end
yy = ylim;
plot([0 0],[yy(1) yy(2)],'k--');
xlim([-180 180])
end
%%
figure; hold on; 
[PPC_surf, h_surf] = deal(cell(length(set_of_inds),1));
for lik = 1:length(set_of_inds)
    
    alltimes = cellfun(@(x) x{lik}(1,:), AG,'UniformOutput',0);
    PPC_surf{lik} = vertcat(alltimes{:});

    h_surf{lik} = surf(1:length(TA),wrapped_cents,PPC_surf{lik}'); 
    set(h_surf{lik},'FaceColor',cop{lik},'FaceAlpha',0.5,'EdgeAlpha',0);
end
