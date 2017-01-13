plot_unit = 31.1;
area = 'PMd';

t1 = -200;
t2 = 2000;

areaunits = eval(sprintf('alldays(1).%s_units',area));
unit_ids = zeros(length(areaunits),1);
for i = 1:length(areaunits)
    unit_ids(i) = areaunits{i}(1);
end
unitidx = find(unit_ids==plot_unit);

combtt = [alldays(1).tt ; alldays(2).tt];

[cents,~,~,~,~,~,dirrast,dirrast_inds] = co_tuning(areaunits,combtt(:,10),combtt,t1,t2,'target',unitidx);

tune_times = [50 250 600 800];
[~,~,~,~,~,~,dirrastSUBTIME] = co_tuning(areaunits,combtt(:,10),combtt,t1,t2,'target',unitidx);

trialnums = cellfun(@(x) size(x,1),dirrast);
trialnums(trialnums<10)=0; max_trialnum = max(sum(trialnums,2));
maxes = max(max(trialnums));
%%
xarray = [t1:(t2-1)];
colplot = {'r','b','k'};
figure; hold on; 
plot_locations = [6 3 2 1 4 7 8 9];
% tune_times = [0 200 800];
tune_periods = find(ismember([t1:t2],tune_times));
polytune = cell(size(dirrast,2),length(tune_periods)-1);
for i = 1:8
    
    subplot(3,3,plot_locations(i)); hold on; 
    
    numtrials = cellfun(@(x) size(x,1),dirrast(i,:));
    %indoffset = [nan numtrials(1:end-1)];
    indoffset = [nan maxes maxes];
    
    for j = 1:size(dirrast,2)
        
        currast = dirrast{i,j};
        currast(currast==0)=nan;
        currast_inds = dirrast_inds{i,j};
        
        for pp = 1:(length(tune_periods)-1)
            if size(currast,1)>10
                polytune{j,pp}(i) = 1000*nansum(nansum(currast(:,(tune_periods(pp)):(tune_periods(pp+1)-1))))./...
                    ((tune_periods(pp+1)-tune_periods(pp))*size(currast,1));
            else
                polytune{j,pp}(i) = nan;
            end
        end
        
        xsgocue = round(1000*(combtt(currast_inds,6)-combtt(currast_inds,5)))';
        xsendtrial = round(1000*(combtt(currast_inds,7)-combtt(currast_inds,5)))';
        
        [~,gocuesort] = sortrows(xsgocue');
        gocuesort = flipud(gocuesort);
        for k = 1:size(currast,1)
            
            xs = find(currast(k,:)==1);
            ys = find(gocuesort==k);
   
            if size(currast,1)>10
                %plot(xarray(xs),(ys+nansum(indoffset(1:j)))*ones(size(xs)),'square','Color',colplot{j},'MarkerSize',.25);
                plottick(xarray(xs),ys*ones(size(xs)),size(currast,1),maxes+2,j,colplot{j});
                plottick(xsgocue(k),ys,size(currast,1),maxes+2,j,'c');
                %plottick(xsendtrial(k),k,size(currast,1),maxes+2,j,'m');
            end

        end
        plot([t1,t2],nansum(indoffset(1:j))*[1 1],'k-');
        xlim([t1 t2]);
        ylim([0 3*maxes]);%sum(numtrials)+1]);
    end
end
%%
% figure; hold on;
% % subplot(3,3,5); hold on; 
% %figure; hold on; 
% for j = 1:3
%     subplot(1,3,j); hold on; 
%     for i = 1:size(polytune,1)
%         polygon_tuning(cents,polytune{i,1}',colplot{i},'-',plot_unit)
% %         polygon_tuning(cents,polytune{i,3}',colplot{i},'--',plot_unit)
%         axis square;
%     end
% axis square;
% axis off;

