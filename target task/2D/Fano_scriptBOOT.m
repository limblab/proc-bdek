%% Create directional bins
dbins = (0:pi/16:2*pi)+pi/32;
DBINS = cell(length(dbins)-1,1);
for i =1:length(dbins)-1
    DBINS{i} = dbins(i:(i+1));
end

fanofunc = @(x) nanvar(x,1)./nanmean(x,1);
crossfunc = @(x,y) cos(x).*sin(y) - sin(x).*cos(y);                     
inbinfunc = @(x,b1,b2) crossfunc(x,b1) < 0 & crossfunc(x,b2) > 0;

FRthresh = 0;
minTrials = 10;

LHvals = flipud(unique(alldays(prediction_day_indices).tt(:,3)));
movedirs = alldays(prediction_day_indices).tt(:,10);
[FAN,FF] = deal(cell(length(LHvals),1));
[BF] = deal(zeros(length(LHvals),3));
colp = {'b','g','r'};
for TIMES = 1:length(spike_count_raw{1})
    
    COUNTS = spike_count_raw{1}{TIMES};
    
    av_rates = mean(COUNTS)./diff(tune_ranges(1:2)/1000);
    keepneur = find(av_rates > FRthresh);
 
    for likei = 1:length(LHvals)
        likecondinds = find(alldays(prediction_day_indices).tt(:,3)==LHvals(likei));
        FULLINDS = cellfun(@(x) find(inbinfunc(movedirs,x(1),x(2))),DBINS,'UniformOutput',0);

        [~,booti] = bootstrp(1000,[],likecondinds);
        fano_run = zeros(size(booti,2),size(COUNTS,2));
        [fano_pop,fano_pop_varweight,fano_pop_ntweight,fano_pop_nelweight,...
            fano_pop_neurnelweight] = deal(zeros(size(booti,2),1));
        for rerun = 1:size(booti,2)
            lci = likecondinds(booti(:,rerun));
            LBIN_counts = cellfun(@(x) COUNTS(lci(ismember(lci,x)),:),FULLINDS,'UniformOutput',0);
            
            oklength = cellfun(@(x) size(x,1),LBIN_counts)>minTrials;
            fanos_bin = cellfun(fanofunc,LBIN_counts(oklength),'UniformOutput',0);
            fano_run(rerun,:) = nanmean(vertcat(fanos_bin{:})); 
            neurel = cellfun(@sum,LBIN_counts(oklength),'UniformOutput',0);
            
            fb_all = horzcat(fanos_bin{:});
            neurel_all = horzcat(neurel{:});
            
            fano_pop_mean = cellfun(@(x) nanmean(fanofunc(x),2), LBIN_counts(oklength),'UniformOutput',0);
            fpm = vertcat(fano_pop_mean{:})';
            fpv = cellfun(@(x) nanvar(fanofunc(x))/(sum(~isnan(fanofunc(x)))),LBIN_counts(oklength))';
            nt = cellfun(@(x) size(x,1),LBIN_counts(oklength))';
            nel = cellfun(@(x) size(x,1).*sum(~isnan(fanofunc(x))),LBIN_counts(oklength))';

            fano_pop_varweight(rerun) = nansum(fpm./fpv)/nansum(1./fpv);
            fano_pop_ntweight(rerun) = nansum(nt.*fpm)/nansum(nt);
            fano_pop_nelweight(rerun) = nansum(nel.*fpm)/nansum(nel);
            fano_pop_neurnelweight(rerun) = nansum((fb_all.*neurel_all))./(sum(neurel_all));

%             fanos_cnt = cellfun(@sum,LBIN_counts(oklength),'UniformOutput',0);
%             vectr_var = cellfun(@(x) (std(x)/sqrt(length(x))).^2,LBIN_counts(oklength),...
%                                                         'UniformOutput',0);
%             fanos_norm = cellfun(@(x) kstest((x-nanmean(x))/nanstd(x)),LBIN_counts(oklength),'UniformOutput',0);
% 
%             fanos_all = horzcat(fanos_bin{:});
%             cnts_all = horzcat(fanos_cnt{:});
%             var_all = horzcat(vectr_var{:}); var_all(var_all==0)=nan;
%             norm_all = horzcat(fanos_norm); norm_all = ~norm_all;
%             
%             FAN{likei}(rerun,TIMES) = nansum(fanos_all./var_all)./nansum(1./var_all);
%             FAN{likei}(rerun,TIMES) = nansum(fanos_all.*cnts_all)./nansum(cnts_all);
%             FAN{likei}(rerun,TIMES) = nanmean(fanos_all);
% 
%             FAN{likei}(rerun,TIMES) = nanmean(cellfun(@nanmean,fanos_bin));
              if rem(rerun,10)==0
                  clc; fprintf('Bin %d/%d - Likelihood %d/%d\n%d/%d\n',...
                    TIMES,length(spike_count_raw{1}),likei,length(LHvals),rerun,size(booti,2));
              end
        end
        FAN{likei}(:,TIMES) = fano_pop_varweight;
        %FAN{likei}(:,TIMES) = fano_pop_ntweight;
        %FAN{likei}(:,TIMES) = fano_pop_nelweight;
        %FAN{likei}(:,TIMES) = fano_pop_neurnelweight;
 
        Finorder = sortrows(FAN{likei}(:,TIMES));
        lfa = Finorder(round(0.025*length(Finorder)));
        ufa = Finorder(round(0.975*length(Finorder)));
        %[lfa,ufa] = boot_bounds(1000,@nanmean,FAN{likei}(:,TIMES),2.5,97.5);
        
        FF{likei}(:,TIMES) = [lfa; ufa; nanmean(FAN{likei}(:,TIMES))];
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
    



