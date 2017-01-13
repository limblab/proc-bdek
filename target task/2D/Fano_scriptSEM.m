%% Create directional bins
dbins = (0:pi/8:2*pi)+pi/16;

fanofunc = @(x) var(x)./mean(x);
crossfunc = @(x,y) cos(x).*sin(y) - sin(x).*cos(y);                     
inbinfunc = @(x,b1,b2) crossfunc(x,b1) < 0 & crossfunc(x,b2) > 0;

FRthresh = 0;
minTrials = 10;

LHvals = flipud(unique(alldays(prediction_day_indices).tt(:,3)));
movedirs = alldays(prediction_day_indices).tt(:,10);
[FANOS,FF] = deal(cell(length(LHvals),1));
[BF] = deal(zeros(length(LHvals),3));
colp = {'b','g','r'};
for TIMES = 1:length(spike_count_raw{1})
    timebinindex = TIMES;
    COUNTS = spike_count_raw{1}{timebinindex};
    
    av_rates = mean(COUNTS)./diff(tune_ranges(1:2)/1000);
    keepneur = find(av_rates > FRthresh);
    
    COUNTS = COUNTS(:,keepneur);
    
    [bininds,LBIN_counts] = deal(cell(length(LHvals),length(dbins)-1));
    for likei = 1:length(LHvals)
        likecondinds = find(alldays(prediction_day_indices).tt(:,3)==LHvals(likei));
        for i = 1:(length(dbins)-1) 
            fullinds = find(inbinfunc(movedirs,dbins(i),dbins(i+1)));
            bininds{likei,i} = fullinds(ismember(fullinds,likecondinds)); 
            LBIN_counts{likei,i} = COUNTS(bininds{likei,i},:);
        end
    end

    % Get bins with at least 10 trials
    tlengths = cellfun(@(x) size(x,1),LBIN_counts);
    basenum = 10;
    %basenum = min(tlengths(tlengths>10));
    normFF = cell(length(LHvals),1);
    for likei = 1:length(LHvals)

        bininds_L = bininds(likei,:);
        tlengths_L = cellfun(@(x) size(x,1),bininds_L);
        goodlengths = find(tlengths_L > minTrials);

        %% Separate neurons
        [fano_neur,fano_sem,fano_semsqr] = deal(zeros(size(COUNTS,2),length(goodlengths)));
        for j = 1:size(COUNTS,2)
            
            for i = 1:length(goodlengths) 
                fano_neur(j,i) = fanofunc(COUNTS(bininds_L{goodlengths(i)},j));
                [fl,fh] = boot_bounds(1000,fanofunc,COUNTS(bininds_L{goodlengths(i)},j),15.9,84.1);
                fano_sem(j,i) = (fh - fl)/2;
                fano_semsqr(j,i) = fano_sem(j,i).^2;
                clc; fprintf('Time Bin %d/%d - Likelihood %d/%d - neuron %d/%d\n',...
                    TIMES,length(spike_count_raw{1}),likei,length(LHvals),j,size(COUNTS,2));
            end
            
            
            normFF{likei}(j) = nansum(fano_neur(j,:)./fano_semsqr(j,:))./sum(1./fano_semsqr(j,:));
            
        end
        [BF(likei,1),BF(likei,2)] = boot_bounds(1000,@mean,normFF{likei}(~isnan(normFF{likei})),2.5,97.5);
        BF(likei,3) = nanmean(normFF{likei});
        
        FF{likei}(:,TIMES) = [BF(likei,[1 2 3])'];
    end
end
%%

figure; hold on; 
for likei = 1:length(LHvals)
    
    if length(tune_ranges)==2
        XS = [1 2];
        LB = repmat(FF{likei}(1,:),1,2);
        UB = repmat(FF{likei}(2,:),1,2);
        M = repmat(FF{likei}(3,:),1,2);
    else
        XS = xsforplot;
        LB = FF{likei}(1,:);
        UB = FF{likei}(2,:);
        M = FF{likei}(3,:);
    end

    plot(XS,M,'LineWidth',3,'Color',colp{likei});
    patch([XS fliplr(XS)],[LB fliplr(UB)],...
          colp{likei},'EdgeAlpha',0,'FaceAlpha',0.25);
end
    



