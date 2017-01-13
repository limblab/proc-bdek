plot_unit = 31.1;
area = 'PMd';

t1 = -100;
t2 = 1500;

areaunits = eval(sprintf('alldays(1).%s_units',area));
unit_ids = zeros(length(areaunits),1);
for i = 1:length(areaunits)
    unit_ids(i) = areaunits{i}(1);
end
unitidx = find(unit_ids==plot_unit);

combtt = [alldays(1).tt ; alldays(2).tt];

[cents,~,~,~,~,~,dirrast,dirrast_inds] = co_tuning(areaunits,combtt(:,10),combtt,t1,t2,'target',unitidx);
trialnums = cellfun(@(x) size(x,1),dirrast);
trialnums(trialnums<10)=0; max_trialnum = max(sum(trialnums,2));
maxes = max(max(trialnums));

%%
tune_times = {[50 250], [300 700], [0 200]};
tune_alignments = {'target','target',12};
trial_threshold = 10;
[dirrastSUBTIME,tune_time_means,tune_time_lows,tune_time_highs] = deal(cell(length(tune_times),1));
for i = 1:(length(tune_times))
    [~,~,~,~,~,~,dirrastSUBTIME{i}] = co_tuning(areaunits,combtt(:,10),...
        combtt,tune_times{i}(1),tune_times{i}(2),tune_alignments{i},unitidx);
    subthreshold = double(cell2mat(cellfun(@(x) size(x,1)>trial_threshold,dirrastSUBTIME{i},'UniformOutput',0)));
    subthreshold(subthreshold==0)=NaN;
    
    tune_time_means{i} = cell2mat(cellfun(@(x) 1000*sum(sum(x))./numel(x),dirrastSUBTIME{i},'UniformOutput',0)).*subthreshold;
    
    trial_rates = cellfun(@(x) 1000*sum(x,2)./size(x,2),dirrastSUBTIME{i},'UniformOutput',0);

    for j = 1:size(trial_rates,1)
        for k = 1:size(trial_rates,2)
            if length(trial_rates{j,k})>trial_threshold
                [tune_time_lows{i}(j,k), tune_time_highs{i}(j,k)] = boot_bounds(1000,@mean,trial_rates{j,k},2.5,97.5);
            else
                tune_time_lows{i}(j,k) = nan;
                tune_time_highs{i}(j,k) = nan;
            end
        end
    end
end

%%
xarray = [t1:(t2-1)];
colplot = {'r','b','k'};
figure; hold on; 
plot_locations = [6 3 2 1 4 7 8 9];
for i = 1:8
    
    subplot(3,3,plot_locations(i)); hold on; 
    
    numtrials = cellfun(@(x) size(x,1),dirrast(i,:));
    indoffset = [nan maxes maxes];
    
    for j = 1:size(dirrast,2)
        
        currast = dirrast{i,j};
        currast(currast==0)=nan;
        currast_inds = dirrast_inds{i,j};
        
        xsgocue = round(1000*(combtt(currast_inds,6)-combtt(currast_inds,5)))';
        xsendtrial = round(1000*(combtt(currast_inds,7)-combtt(currast_inds,5)))';
        xsRTcue = round(1000*(combtt(currast_inds,12)-combtt(currast_inds,5)))';
        xsRTgo = round(1000*(combtt(currast_inds,12)-combtt(currast_inds,6)))';
        
        [~,gocuesort] = sortrows(xsgocue');
%         [~,gocuesort] = sortrows(xsRTgo');
        gocuesort = flipud(gocuesort);
        for k = 1:size(currast,1)
            
            xs = find(currast(k,:)==1);
            ys = find(gocuesort==k);
   
            if size(currast,1)>10
                plottick(xarray(xs),ys*ones(size(xs)),size(currast,1),maxes+2,j,colplot{j});
                plottick(xsgocue(k),ys,size(currast,1),maxes+2,j,'c');
                [YRT(1), YRT(2)] = plottick(xsRTcue(k),ys,size(currast,1),maxes+2,j,'m');
                plot([xsRTcue(k) xsRTcue(k)+200],[1 1]*mean(YRT),'Color',[.5 .5 .5]);
                plottick(xsendtrial(k),ys,size(currast,1),maxes+2,j,'y');
                
            end

        end
        plot([t1,t2],nansum(indoffset(1:j))*[1 1],'k-');
        xlim([t1 t2]);
        ylim([0 3*maxes]);
    end
end
%%
figure; hold on; 
for j = 1:length(tune_times)
    subplot(1,length(tune_times),j); hold on; 
    ext = zeros(size(tune_time_means,1),3);
    for i = 1:size(tune_time_means,1)
%         [ext(i,1)] = polygon_tuning(cents,tune_time_means{j}(:,i),colplot{i},'-',plot_unit);
%         [ext(i,2)] = polygon_tuning(cents,tune_time_lows{j}(:,i),colplot{i},'--',plot_unit);
%         [ext(i,3)] = polygon_tuning(cents,tune_time_highs{j}(:,i),colplot{i},'--',plot_unit);
%         
        [ext(i,3)] = polygon_tuning_patch(cents,...
            [tune_time_means{j}(:,i) tune_time_lows{j}(:,i) tune_time_highs{j}(:,i)],...
            colplot{i},'-',1,plot_unit); 
    end
    xlim([-max(ext(:)) max(ext(:))]);
    ylim([-max(ext(:)) max(ext(:))]);
    axis square;
end
% axis square;
% axis off;

