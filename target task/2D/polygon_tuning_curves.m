path_name = 'C:\Users\limblab\Desktop\Polygons\MrT\M1\Run2\';
plot_unit = 'all';

area = 'M1';
t1 = -200;
t2 = 1500;
plot_periods = [1 3];

help_fun = @(x) 1000*nansum(nansum(x))./(size(x,2)*size(x,1));
if isnumeric(plot_unit)

    areaunits = eval(sprintf('alldays(1).%s_units',area));
    unit_ids = zeros(length(areaunits),1);
    for i = 1:length(areaunits)
        unit_ids(i) = areaunits{i}(1);
    end
    unitidx = find(unit_ids==plot_unit);

    combtt = [alldays(1).tt ; alldays(2).tt];

    [cents,~,~,~,~,~,dirrast] = co_tuning(areaunits,combtt(:,10),combtt,t1,t2,'target',unitidx);

    trialnums = cellfun(@(x) size(x,1),dirrast);
    trialnums(trialnums<10)=0; max_trialnum = max(sum(trialnums,2));
    maxes = max(max(trialnums));
    %%
    xarray = [t1:(t2-1)];
    colplot = {'r','b','k'};
    % tune_times = [0 200 800];
    tune_times = [50 250 600 800];
    tune_periods = find(ismember([t1:t2],tune_times));
    polytune = cell(size(dirrast,2),length(tune_periods)-1);
    for i = 1:8

        numtrials = cellfun(@(x) size(x,1),dirrast(i,:));
        indoffset = [nan maxes maxes];

        for j = 1:size(dirrast,2)

            currast = dirrast{i,j};
            currast(currast==0)=nan;

            for pp = 1:(length(tune_periods)-1)
                if size(currast,1)>10
                    polytune{j,pp}(i) = 1000*nansum(nansum(currast(:,(tune_periods(pp)+1):(tune_periods(pp+1)-1))))./...
                        ((tune_periods(pp+1)-tune_periods(pp))*size(currast,1));
                else
                    polytune{j,pp}(i) = nan;
                end
            end
        end
    end
    %%
    f = figure; hold on; title(sprintf('%.1f',plot_unit),'FontSize',16);
    %figure; hold on; 
    for i = 1:length(plot_periods)
        subplot(1,length(plot_periods),i); hold on;
        for j = 1:size(polytune,1)
            polygon_tuning(cents,polytune{j,plot_periods(i)}',colplot{j},'-',plot_unit)
            axis square;
        end
    end
    axis square;

else
    
    areaunits = eval(sprintf('alldays(1).%s_units',area));
    for units = 1:length(areaunits)
    
       
        unitidx = units;
        plot_unit = areaunits{units}(1);

        combtt = [alldays(1).tt ; alldays(2).tt];

        [cents,~,~,~,~,~,dirrast] = co_tuning(areaunits,combtt(:,10),combtt,t1,t2,'target',unitidx);

        trialnums = cellfun(@(x) size(x,1),dirrast);
        trialnums(trialnums<10)=0; max_trialnum = max(sum(trialnums,2));
        maxes = max(max(trialnums));
        %%
        xarray = [t1:(t2-1)];
        colplot = {'r','b','k'};
        % tune_times = [0 200 800];
        tune_times = [50 250 600 800];
        tune_periods = find(ismember([t1:t2],tune_times));
        [polytune,polytune_L,polytune_H] = deal(cell(size(dirrast,2),length(tune_periods)-1));
        for i = 1:8

            numtrials = cellfun(@(x) size(x,1),dirrast(i,:));
            indoffset = [nan maxes maxes];

            for j = 1:size(dirrast,2)

                currast = dirrast{i,j};
                currast(currast==0)=nan;

                for pp = 1:(length(tune_periods)-1)
                    if size(currast,1)>10
                        polytune{j,pp}(i) = help_fun(currast(:,(tune_periods(pp)+1):(tune_periods(pp+1))));
                        [polytune_L{j,pp}(i),polytune_H{j,pp}(i)] = ...
                            boot_bounds(1000,help_fun,currast(:,(tune_periods(pp)+1):(tune_periods(pp+1))),2.5,97.5);
                        
                  
%                         polytune{j,pp}(i) = 1000*nansum(nansum(currast(:,(tune_periods(pp)+1):(tune_periods(pp+1)-1))))./...
%                             ((tune_periods(pp+1)-tune_periods(pp))*size(currast,1));
                    else
                        polytune{j,pp}(i) = nan;
                        polytune_L{j,pp}(i) = nan;
                        polytune_H{j,pp}(i) = nan;
                    end
                end
            end
        end
        %%
        f = figure; hold on; title(sprintf('%.1f',plot_unit),'FontSize',16);
        %figure; hold on; 
        for i = 1:length(plot_periods)
            subplot(1,length(plot_periods),i); hold on;
            for j = 1:size(polytune,1)
                polygon_tuning(cents,polytune{j,plot_periods(i)}',colplot{j},'-',plot_unit)
                polygon_tuning(cents,polytune_L{j,plot_periods(i)}',colplot{j},'--',plot_unit)
                polygon_tuning(cents,polytune_H{j,plot_periods(i)}',colplot{j},'--',plot_unit)
                axis square;
                axis equal;
            end
        end
        axis square;
        axis equal;
        saveas(f,sprintf('%sUnit_%d_%d.png',path_name,floor(plot_unit),round(10*(plot_unit-floor(plot_unit)))));
        saveas(f,sprintf('%sUnit_%d_%d.eps',path_name,floor(plot_unit),round(10*(plot_unit-floor(plot_unit)))));
        close;
    end
end





