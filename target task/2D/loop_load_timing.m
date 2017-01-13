%%
% Set up list of all days for inclusion
FileName = {'Mihili','07112013',  2        ;...
            'Mihili','07152013', [2,3]     ;...
            'Mihili','07192013', [2,3,4]   ;...
            'Mihili','08062013', [2,3]     ;...
            'Mihili','08122013', [2,3]     ;...
            'Mihili','08152013',  2        ;...
            'MrT'   ,'05042013',  2        ;...
            'MrT'   ,'05052013',  2        ;...
            'MrT'   ,'05062013',  2        };
        
[TL, TU, TM, TC] = deal(cell(size(FileName,1),1)); %initialize
Deneves = zeros(size(FileName,1),2);
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2}); 
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    Activity_script2; % Obtain spike counts, etc. 
    hinds = find(alldays(2).tt(:,3)==min(alldays(2).tt(:,3)));
    
    [FDE,FDL] = deal(zeros(size(firing_diffs{1}{1},2),2));

    [FDE(:,1), FDE(:,2)] = boot_bounds(1000,@nanmean,firing_diffs{1}{1}(hinds,:),2.5,97.5);
    [FDL(:,1), FDL(:,2)] = boot_bounds(1000,@nanmean,firing_diffs{1}{2}(hinds,:),2.5,97.5);
    
    Deneves(daynum,1) = sum(FDE(:,2)<0 & FDL(:,1)>0);
    Deneves(daynum,2) = sum(~isnan(FDE(:,1)));
    
   
    % Fill variables
%     TL{daynum} = l_alltimes;
%     TU{daynum} = u_alltimes;
%     TM{daynum} = av_alltimes;
    
%     TL{daynum} = l_alldeltas;
%     TU{daynum} = u_alldeltas;
%     TM{daynum} = av_alldeltas;
%     for i = 1:length(like_ind{1})
%         TC{daynum,i} = T_change{i};
%     end
%     
% %     if length(like_ind{1})==3
% %         TL{daynum}(2:3) = TL{daynum}([3 2]);
% %         TU{daynum}(2:3) = TU{daynum}([3 2]);
% %         TM{daynum}(2:3) = TM{daynum}([3 2]);
% %     end

    fprintf('%d/%d\n',daynum,size(FileName,1)); 
end

%%
% Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
% MrT_inds = find(strcmp(FileName(:,1),'MrT'));
% 
% figure; hold on; 
% 
% cfp = {'b','g','r'};
% for daynum = 1:size(FileName,1)
%     if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
%     for i = 1:size(TM{daynum},1)
%         plot(TM{daynum}{i},daynum + x_off,'.','Color',cfp{i},'MarkerSize',10);
%         plot([TL{daynum}{i} TU{daynum}{i}],daynum*[1 1]+x_off,'Color',cfp{i},'LineWidth',1);
%     end
% end
% ylim([-2 daynum+3+3]);
% 
% 
% cfp = {'b','g','r'};
% for i = 1:2
%     %Mihili
%     allvals = vertcat(TC{Mihili_inds,i});
%     avvals = mean(allvals); 
%     [l_vals, u_vals] = boot_bounds(1000,@mean,allvals,2.5,97.5);
%  
%     plot(avvals*[1 1], 1.5*length(Mihili_inds*[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
%     plot(1.5+length(Mihili_inds),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
%     
%     %Mr T
%     allvals = vertcat(FRC{MrT_inds,i});
%     avvals = mean(allvals); 
%     [l_vals, u_vals] = boot_bounds(1000,@mean,allvals,2.5,97.5);
% 
%     plot([1 1]*(4.5+size(FileName,1)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
%     plot(4.5+size(FileName,1),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
%     
% end

%% 
% Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
% MrT_inds = find(strcmp(FileName(:,1),'MrT'));
% 
% figure; hold on; 
% 
% for daynum = 1:size(FileName,1)
%     if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
%     if size(TM{daynum},1) == 3; cfp = {'b','g','r'}; else cfp = {'b','g'}; end
%     for i = 1:size(TM{daynum},1)
%         plot(daynum + x_off, TM{daynum}{i},'.','Color',cfp{i},'MarkerSize',10);
%         plot(daynum*[1 1]+x_off, [TL{daynum}{i} TU{daynum}{i}],'Color',cfp{i},'LineWidth',1);
%     end
% end
% xlim([-2 daynum+3+3]);
% 
% cfp = {'b','g','r'};
% for i = 1:2
%     %Mihili
%     allvals = vertcat(TC{Mihili_inds,i});
%     avvals = mean(allvals); 
%     [l_vals, u_vals] = boot_bounds(1000,@mean,allvals,2.5,97.5);
%  
%     plot([1 1]*(1.5+length(Mihili_inds)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
%     plot(1.5+length(Mihili_inds),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
%     
%     %Mr T
%     allvals = vertcat(TC{MrT_inds,i});
%     avvals = mean(allvals);
%     [l_vals, u_vals] = boot_bounds(1000,@mean,allvals,2.5,97.5);
% 
%     plot([1 1]*(4.5+size(FileName,1)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
%     plot(4.5+size(FileName,1),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
%     
% end
