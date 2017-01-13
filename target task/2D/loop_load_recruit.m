
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
        
        
 %FileName = {'Mihili','08152013',  2        };

        
[AG,CR,UR, AGU, AGL, GA, P] = deal(cell(size(FileName,1),1)); %initialize
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2});
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5);
    
    %MAKE SURE FR_change{z} = nanmean(firing_perc{z}{bin},2);
    Activity_script2; close; % Obtain spike counts, etc. 
    PPC;  close; % Run population recruitment code
    recruitment_interp; close;
    
    % Fill variables
    AG{daynum} = av_gain;
    CR{daynum} = co_recruit;
    UR{daynum} = un_recruit;
    AGL{daynum} = bL; 
    AGU{daynum} = bU;
    GA{daynum} = GAIN;
    P{daynum} = p;
    
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
end

%%
% figure; hold on;
% for daynum = 1:size(FileName,1);
%     ming = min(horzcat(AG{daynum}{:}));
%     maxg = max(horzcat(AG{daynum}{:}));
%     
%     %plot(vertcat(AG{daynum}{:})');
%     %plot((vertcat(AG{daynum}{:})'-ming)./(maxg-ming));
%     
%     gains = vertcat(AG{daynum}{:});
%     gains_h = 0.5*(gains(:,1:314)+fliplr(gains(:,315:end)));
%     gains_s = [gains_h fliplr(gains_h)];
%     
%     dgain = diff(gains);
%     dgain_s = 0.5*(gains(:,1:314)+fliplr(gains(:,315:end)));
%     
%     bind_gain = bin_array(gains_h,size(gains_h,1),4);
%     a_da = [-bind_gain(:,end) , bind_gain(:,1)];
%     
%     cop = {'b','g','r'};
%     for lik = 1:size(bind_gain,1)   
%         plot(a_da(lik,:),0.3*(lik-1)+[daynum daynum],cop{lik},'LineWidth',4);
%     end    
%     plot([0 0],[-20 20],'k');
%     
% end
%%
num_degrees = 90;

ind_align = [(315-round(100*pi*num_degrees/360)):(314+round(100*pi*num_degrees/360))];
ind_nonalign = [1:round(100*pi*num_degrees/360) 628:-1:(629-round(100*pi*num_degrees/360))];

[CROSS, DOT] = deal(cell(size(FileName,1),1));
cfp = {'b','r','g'};
for daynum = 1:size(FileName,1)
    for lik = 1:size(GA{daynum},1)
        
%         comb_align = GA{daynum}{lik}(:,ind_align); comb_align = comb_align(:);
%         comb_nonalign = GA{daynum}{lik}(:,ind_nonalign); comb_nonalign = comb_nonalign(:);
        comb_align = GA{daynum}{lik}(:,ind_align); comb_align = mean(comb_align,2);
        comb_nonalign = GA{daynum}{lik}(:,ind_nonalign); comb_nonalign = mean(comb_nonalign,2);
        
        align_m = mean(comb_align);
        [align_l,align_u] = boot_bounds(1000,@mean,comb_align,2.5,97.5);
        nonalign_m = mean(comb_nonalign);
        [nonalign_l,nonalign_u] = boot_bounds(1000,@mean,comb_nonalign,2.5,97.5);
        
        CROSS{daynum}(lik,:) = [align_l align_u, nonalign_l, nonalign_u];
        DOT{daynum}(lik,:) = [align_m, nonalign_m];
    end 
end
%%
figure; hold on;
for daynum = 1:size(FileName,1)
    if strcmp(FileName{daynum,1},'Mihili'); markr = '-'; else markr = ':';end
    plot(DOT{daynum}(:,1),DOT{daynum}(:,2),markr,'Color',[0.5 0.5 0.5],'LineWidth',2);
end

for daynum = 1:size(FileName,1)
    if size(GA{daynum},1)==3; cfp = {'b','g','r'}; else cfp = {'b','g'}; end
    for lik = 1:size(GA{daynum},1)
        plot(CROSS{daynum}(lik,1:2),[1 1]*DOT{daynum}(lik,2),cfp{lik},'LineWidth',4);
        plot([1 1]*DOT{daynum}(lik,1),CROSS{daynum}(lik,3:4),cfp{lik},'LineWidth',4);
    end  
end
extremey = ylim; extremex = xlim;
plot([-1 max([extremey extremex])],[-1 max([extremey extremex])],'k--');
plot([-1 max([extremey extremex])],[0 0],'k');
plot([0 0],[-1 max([extremey extremex])],'k');

    