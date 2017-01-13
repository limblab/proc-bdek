fanofunc = @(x) var(x)./mean(x);
crossfunc = @(x,y) cos(x).*sin(y) - sin(x).*cos(y);                     
inbinfunc = @(x,b1,b2) crossfunc(x,b1) < 0 & crossfunc(x,b2) > 0;
cobins = (0:pi/4:2*pi) + pi/8;

FRthresh = 0;
NEL = 10;

TT = alldays(prediction_day_indices).tt;
priord = circ_mean(TT(:,2));

LHvals = flipud(unique(TT(:,3)));
movedirs = TT(:,10);
[FANOS,FF] = deal(cell(length(LHvals),1));
[BF] = deal(zeros(length(LHvals),3));

if sum(LHvals > 1000)
    colp = {'k','b','g','r'};
else
    colp = {'b','g','r'};
end

for TIMES = 1:length(spike_count_raw{1})
    timebinindex = TIMES;
    COUNTS = spike_count_raw{1}{timebinindex};
    
    [LBIN_counts,FANOS,ALLFANOS] = deal(cell(length(LHvals),1));
    for likei = 1:length(LHvals)
        
        % Find indices for current likelihood condition
        likecondinds = find(TT(:,3)==LHvals(likei));
        COUNTS_like = COUNTS(likecondinds,:);
        TT_like = TT(likecondinds,:);
        
        if LHvals(likei) > 1000 % If we're in center out
            for cobin = 1:(length(cobins)-1) % loop through CO targets
                
                % Find CO trials within current target
                incobin = find(inbinfunc(TT_like(:,10),cobins(cobin),cobins(cobin+1)));
                TT_cobin = TT_like(incobin,:);
                COUNTS_cobin = COUNTS_like(incobin,:);
                
                % Sort the reaches by distance to center of the target
                [comoves_ordered,coorderind] = ...
                    sortrows(circ_dist(TT_cobin(:,10),cobins(cobin)+pi/8));
                
                % Separate into groups of desired trial length
                coorder_pad = [coorderind; nan(NEL-rem(length(coorderind),NEL),1)];
                cogroups = reshape(coorder_pad,NEL,[]);
                cogroups(:,isnan(sum(cogroups))) = [];
                for ci = 1:size(cogroups,2) % Loop through n-trial groups
                    cogroupmoves = TT_cobin(cogroups(:,ci),10);
                    cogroupspan = abs(circ_dist(cogroupmoves(1),cogroupmoves(end)));
                    
                    if cogroupspan < pi/8
                        % Assign spike counts
                        LBIN_counts{likei}{cobin,ci} = COUNTS_cobin(cogroups(:,ci),:);
                    else
                        LBIN_counts{likei}{cobin,ci} = NaN;
                    end
                end
                
            end
        else

            [moves_ordered,orderind] = sortrows(circ_dist(TT_like(:,10),priord));
            order_pad = [orderind; nan(NEL-rem(length(orderind),NEL),1)];
            groups = reshape(order_pad,NEL,[]);
            groups(:,isnan(sum(groups))) = [];

            for i = 1:size(groups,2)
                groupmoves = TT_like(groups(:,i),10); 
                groupspan = abs(circ_dist(groupmoves(1),groupmoves(end)));

                if groupspan < pi/8
                    LBIN_counts{likei}{i} = COUNTS_like(groups(:,i),:);
                else
                    LBIN_counts{likei}{i} = NaN;
                end
            end
        end

        FANOS{likei} = cellfun(fanofunc,LBIN_counts{likei},'UniformOutput',0);
        ALLFANOS{likei} = horzcat(FANOS{likei}{:});
    
        [BF(likei,1),BF(likei,2)] = boot_bounds(1000,@nanmean,ALLFANOS{likei},2.5,97.5);
        BF(likei,3) = nanmean(ALLFANOS{likei});

        FF{likei}(:,TIMES) = BF(likei,[1 2 3])';
    
    end
 
end
%%

figure; hold on; 
for likei = 1:length(LHvals)
 
    XS = xsforplot;
    LB = FF{likei}(1,:);
    UB = FF{likei}(2,:);
    M = FF{likei}(3,:);

    if length(tune_ranges)==2     
        XS = [1 2];
        LB = repmat(LB,1,2);
        UB = repmat(UB,1,2);
        M = repmat(M,1,2); 
    end

    plot(XS,M,'LineWidth',3,'Color',colp{likei});
    patch([XS fliplr(XS)],[LB fliplr(UB)],...
          colp{likei},'EdgeAlpha',0,'FaceAlpha',0.25);
end
    
%%
if sum(LHvals > 1000)
    figure; hold on;
    for likei = 2:length(LHvals)

        XS = xsforplot;
        LB = FF{likei}(1,:) - FF{1}(3,:);
        UB = FF{likei}(2,:) - FF{1}(3,:);
        M = FF{likei}(3,:) - FF{1}(3,:);    
        
        if length(tune_ranges)==2
            XS = [1 2];
            LB = repmat(LB,1,2);
            UB = repmat(UB,1,2);
            M = repmat(M,1,2);
        end

        plot(XS,M,'LineWidth',3,'Color',colp{likei});
        patch([XS fliplr(XS)],[LB fliplr(UB)],...
              colp{likei},'EdgeAlpha',0,'FaceAlpha',0.25);
    
    end
end
