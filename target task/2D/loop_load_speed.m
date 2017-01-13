%%
% % Set up list of all days for inclusion
% FileName = {'Mihili','07112013',  [1 2]        ;...
%             'Mihili','07152013',  [1,2]     ;...
%             'Mihili','07192013',  [1,2]   ;...
%             'Mihili','08062013',  [1,2]     ;...
%             'Mihili','08122013',  [1,2]     ;...
%             'Mihili','08152013',  [1,2]        ;...
%             'MrT'   ,'05042013',  [1,2]        ;...
%             'MrT'   ,'05052013',  [1,2]        ;...
%             'MrT'   ,'05062013',  [1,2]        };
        
if ~exist('FileName','var')
FileName = {'Mihili','07112013',  2,  [1, 2]    ;...
    
            'Mihili','07152013',  2,  [1, 2]    ;...
            'Mihili','07152013',  3,  [1, 2]    ;...
            
            'Mihili','07192013',  2,  [1, 2]    ;...
            'Mihili','07192013',  3,  [1, 2]    ;...
            'Mihili','07192013',  5,  [1, 2]    ;...
            
            'Mihili','08062013',  2,  [1, 2]    ;...
            'Mihili','08062013',  3,  [1, 2]    ;...
            
            'Mihili','08122013',  2,  [1, 2]    ;...
            'Mihili','08122013',  3,  [1, 2]    ;...
            
            'Mihili','08152013',  2,  [1, 2]    ;...
            
            'MrT'   ,'05042013',  2,  [1, 2]    ;...
            'MrT'   ,'05052013',  2,  [1, 2]    ;...
            'MrT'   ,'05062013',  2,  [1, 3]    };
end
% FileName = {'Mihili','07192013',  [1,2]   };
           
[SL, SU, SM, FRC, COS,SPD] = deal(cell(size(FileName,1),1)); %initialize
colop = {'k','b','g','r'};
figure; hold on;
for daynum = 1:size(FileName,1)
    
    nameoffile = sprintf('C:\\Users\\limblab\\Desktop\\Mihili_bdf\\full alldays\\full_alldays_%s_%s', FileName{daynum,1},FileName{daynum,2});
    
    % Load file
    load(nameoffile);

    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    alldays(2).tt(alldays(2).tt(:,3)>100,:) = [];

    llist = flipud(unique(alldays(2).tt(:,3))); 
    alldays(2).tt(~ismember(alldays(2).tt(:,3),llist(FileName{daynum,4})),:) = [];
    
    speed_script;
    
    SPD{daynum} = topspeed;

%     SL{daynum} = V(2,2:end);
%     SU{daynum} = V(3,2:end);
%     SM{daynum} = V(1,2:end);
%     COS{daynum} = V([2 3 1],1);
%     % Fill variables
%     for li = 1:length(likes)
%         plot(daynum*[1 1] + li*0.1,[V(2,li) V(3,li)],colop{li},'LineWidth',4);
%         plot(daynum + li*0.1,V(1,li),'.','Color',colop{li},'MarkerSize',10);
%     end
%     
%     plot_cross(V(1,
    
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
end
%loop_load_behavior; close;
%%
spd_diffs = zeros(length(SPD),3);
figure; hold on; 
av_spd = zeros(length(SPD),length(LI{1}));
for i = 1:length(SPD)
    
    daynum = i;
    
    for j = 1:length(LI{daynum})
        
        is = LI{daynum}{j};
        av_spd(daynum,j) = mean(SPD{daynum}(is));
    end
    
    likemat = [SPD{daynum}(LI{daynum}{1}), ones(size(LI{daynum}{1})) ; ...
               SPD{daynum}(LI{daynum}{2}), 2*ones(size(LI{daynum}{2}))];
           
    difffunc = @(x) nanmean(x(x(:,end)==2,1)) - nanmean(x(x(:,end)==1,1));
    
    [l,h,brnd] = boot_bounds(1000,difffunc,likemat,2.5,97.5);
    
    spd_diffs(daynum,:) = [l h nanmean(brnd)];
    
    plot(daynum.*[1 1],[l h],'k'); plot(daynum,spd_diffs(daynum,3),'k.');
end
%%
figure; hold on; 
for i = 1:length(G)
    
    plot(i.*[1 1],-spd_diffs(G(i),1:2),'k');
    plot(i,-spd_diffs(G(i),3),'k.');
end
    
%% 
% Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
% MrT_inds = find(strcmp(FileName(:,1),'MrT'));
% 
% figure; subplot(2,1,1); hold on; 
% for daynum = 1:length(Mihili_inds)
%     di = Mihili_inds(daynum);
%     if size(SM{daynum},2) == 3; cfp = {'b','g','r'}; else cfp = {'b','g'}; end
%     for i = 1:size(SM{di},2)
%         plot([SL{di}(i) SU{di}(i)],[1 1]*FM{di}(i),'Color',cfp{i});
%         plot([1 1]*SM{di}(i),[FL{di}(i) FU{di}(i)],'Color',cfp{i});
%         plot(SM{di}(i),FM{di}(i),'.','Color',cfp{i});
%     end
% end
%     
% subplot(2,1,2); hold on; 
% for daynum = 1:length(MrT_inds)
%     di = MrT_inds(daynum);
%     if size(SM{daynum},2) == 3; cfp = {'b','g','r'}; else cfp = {'b','g'}; end
%     for i = 1:size(SM{di},2)
%         plot([SL{di}(i) SU{di}(i)],[1 1]*FM{di}(i),'Color',cfp{i});
%         plot([1 1]*SM{di}(i),[FL{di}(i) FU{di}(i)],'Color',cfp{i});
%         plot(SM{di}(i),FM{di}(i),'.','Color',cfp{i});
%     end
% end

