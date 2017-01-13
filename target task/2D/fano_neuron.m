%% Create directional bins
dbins = (0:pi/8:2*pi)+pi/16;

fanofunc = @(x) var(x)./mean(x);
crossfunc = @(x,y) cos(x).*sin(y) - sin(x).*cos(y);                     
inbinfunc = @(x,b1,b2) crossfunc(x,b1) < 0 & crossfunc(x,b2) > 0;

FRthresh = 0;
countthresh = 25;
minTrials = 0;    
basenum = 10;

LHvals = flipud(unique(alldays(prediction_day_indices).tt(:,3)));
movedirs = alldays(prediction_day_indices).tt(:,10);
[FANOS,FF] = deal(cell(length(LHvals),1));
[BF] = deal(zeros(length(LHvals),3));
colp = {'b','g','r'};

frintbin = cellfun(@(x) mean(x)./(diff(tune_ranges(1:2)/1000)),spike_count_raw{1},'UniformOutput',0);
allfrs = vertcat(frintbin{:});
maxfrs = max(allfrs,[],1);
for TIMES = 1:length(spike_count_raw{1})
    timebinindex = TIMES;
    COUNTS = spike_count_raw{1}{timebinindex};
    
%     av_rates = mean(COUNTS)./diff(tune_ranges(1:2)/1000);
%     keepneur = find(av_rates > FRthresh);
    keepneur = find(maxfrs > FRthresh);
    
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
    %basenum = min(tlengths(tlengths>10));
    for likei = 1:length(LHvals)
        bininds_L = bininds(likei,:);
        tlengths_L = cellfun(@(x) size(x,1),bininds_L);
        goodlengths = find(tlengths_L > minTrials);
        %% Set up "conditions"
        newcond = cell(length(goodlengths),1);
        for i = 1:length(goodlengths) % directional bins
            binindex = goodlengths(i);
            newnumconds = 1;%floor(tlengths_L(binindex)/basenum);
            shuffd = bininds{likei,binindex}(randperm(length(bininds{likei,binindex})));
            for j = 1:newnumconds
                if j<newnumconds
                    newcond{i}{j} = shuffd(((j-1)*basenum + 1):(j*basenum));
                else
                    newcond{i}{j} = shuffd(((j-1)*basenum + 1):end);
                end
            end
        end
        allconds = horzcat(newcond{:});

        %% Separate neurons
        fano_cond = zeros(length(allconds),size(COUNTS,2));
        for i = 1:length(allconds) 
            for j = 1:size(COUNTS,2)      
                if sum(COUNTS(allconds{i},j)) > countthresh
                    fano_cond(i,j) = fanofunc(COUNTS(allconds{i},j));
                else
                    fano_cond(i,j) = NaN;
                end
            end
        end
        FANOS{likei} = reshape(fano_cond,numel(fano_cond),1);
        [BF(likei,1),BF(likei,2)] = boot_bounds(1000,@nanmean,FANOS{likei},2.5,97.5);
        BF(likei,3) = nanmean(FANOS{likei});
        
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
    



