%%
% Set up list of all days for inclusion
% FileName = {'Mihili','07112013',  2        ;...
%             'Mihili','07152013', [2,3]     ;...
%             'Mihili','07192013', [2,3,4]   ;...
%             'Mihili','08062013', [2,3]     ;...
%             'Mihili','08122013', [2,3]     ;...
%             'Mihili','08152013',  2        ;...
%             'MrT'   ,'05042013',  2        ;...
%             'MrT'   ,'05052013',  2        ;...
%             'MrT'   ,'05062013',  2        };
        
% FileName = {'Mihili','09282015', 2};
%   
FileName = {'Mihili','07112013',  2        ;...
            'Mihili','07152013',  2        ;...
            'Mihili','07192013',  2        ;...
            'Mihili','08062013',  2        ;...
            'Mihili','08122013',  2        ;...
            'Mihili','08152013',  2        ;...
            'MrT'   ,'05042013',  2        ;...
            'MrT'   ,'05052013',  2        ;...
            'MrT'   ,'05062013',  2        };
%   
% FileName = {'Mihili','07192013',  2};
% 
% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10122015',  2        ;...
%             'Mihili','10272015',  2        };
        
% CONTROLS
% FileName = {'Mihili','10082015',  2        ;...
%             'Mihili','10272015',  2        ;...
%             'Mihili','11022015',  2        };

% FileName = {'Mihili','05062014', 2};
        
BRAIN_AREA = 'PMd';
[UDIFFS,UDIFFS_N, BPDS] = deal(cell(size(FileName,1),1)); %initialize
for daynum = 1:size(FileName,1)
    
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2},0); 
    alldays(1).tt(isnan(alldays(1).tt(:,3)),3) = alldays(1).tt(find(isnan(alldays(1).tt(:,3)))-1,3);
    if isfield(alldays,'bdfM')
        BDF = alldays(1).bdfM;
    elseif isfield(alldays,'kin')
        BDF = alldays(1).kin;
    end
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    alldays(2).tt(alldays(2).tt(:,3)>100,:) = [];
    
    if exist(sprintf('%s_PD90_%s_%s_full.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'file')
        load(sprintf('%s_PD90_%s_%s_full.mat',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}));
    else
        brain_area = BRAIN_AREA;
        CO_index = 1;
        reachdir_col = 10;
        loop_alignmentsPD = {'target','go'};
        loop_rangesPD = {[600 800],[50 250]};
        [~,pref,alltuned] = tuning_types_PD_func(alldays,brain_area,CO_index,reachdir_col,loop_alignmentsPD,loop_rangesPD);
        save(sprintf('%s_PD90_%s_%s_full',BRAIN_AREA,FileName{daynum,1},FileName{daynum,2}),'pref','alltuned');
    end
    prs = nan(size(alltuned)); 
    bpds = nan(size(alltuned,1),1);
    for i = 1:size(prs,1);
        if ~isnan(pref(i))
            prs(i,pref(i)) = alltuned(i,pref(i)); 
            bpds(i) = alltuned(i,pref(i));
        end
    end
    best_PDS = bpds;%prs(:,2);
    
    pretarget_plot_PDOD; close;
    
    UDIFFS{daynum} = uncdiff;
    UDIFFS_N{daynum} = uncdiff_n;
    
    clearvars -except FileName UDIFFS UDIFFS_N BPDS BRAIN_AREA;
    
end
%%
figure; hold on; 
for j = 1:length(UDIFFS{1})
    subplot(1,length(UDIFFS{1}),j); hold on; 
    for i = 1:length(UDIFFS)
        plot(UDIFFS{i}{j}(3,:),'b');
        patch([1:size(UDIFFS{i}{j},2) size(UDIFFS{i}{j},2):-1:1],...
            [UDIFFS{i}{j}(1,:) fliplr(UDIFFS{i}{j}(2,:))],'k','FaceAlpha',0.25,'EdgeAlpha',0);
    end
end

%%
figure; hold on; 
subplot(1,2,1); hold on; 
orders = [1 3 2];
for i = 1:length(UDIFFS_N)
    for j = 1:length(orders)
        if UDIFFS_N{i}{j}(1,2) > 0 || UDIFFS_N{i}{j}(2,2) < 0
            sigcol = 'k';
        else
            sigcol = [.5 .5 .5];
        end
        
        plot([1,1]*orders(j),UDIFFS_N{i}{j}(1:2,2),'Color',sigcol,'LineWidth',2);
    end
    ordered = cellfun(@(x) x(3,2),UDIFFS_N{i});
    plot(1:3,ordered([1 3 2]));
end

subplot(1,2,2); hold on; 
for i = 1:length(UDIFFS_N)
    for j = 1:length(orders)
        if UDIFFS_N{i}{j}(1,3) > 0 || UDIFFS_N{i}{j}(2,3) < 0
            sigcol = 'k';
        else
            sigcol = [.5 .5 .5];
        end
        plot([1,1]*orders(j),UDIFFS_N{i}{j}(1:2,3),'Color',sigcol,'LineWidth',2);
    end
    ordered = cellfun(@(x) x(3,3),UDIFFS_N{i});
    plot(1:3,ordered([1 3 2]));
end



%%
% %%
% figure; hold on; 
% 
% for i = 1:length(UDIFFS)
%     
%     plot3(UDIFFS{i}{1}(3,:),UDIFFS{i}{2}(3,:),UDIFFS{i}{3}(3,:),'.-');
%     xlabel('SD'); ylabel('OD'); zlabel('ORTH');
% end
% 
% 

