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
        
FileName = {'Mihili','05062014',  2        ;...
            'Mihili','11212013',  2        };
        
% FileName = {'Mihili','07192013',  2};
        
[FL, FU, FM, FRC] = deal(cell(size(FileName,1),1)); %initialize
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2}); 
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    [~,GOOD_list] = ISI_exam(alldays(2).PMd_units,1.7,0.1);
    Activity_script2; close; % Obtain spike counts, etc. 
   
    % Fill variables
    FL{daynum} = lowerbound{1};
    FU{daynum} = upperbound{1};
    FM{daynum} = midbound{1};
    for i = 1:length(like_ind{1})
       % FRC{daynum,i} = FR_change{1}(like_ind{1}{i});
        FRC{daynum,i} = firing_diffs{1}{1}(like_ind{1}{i},:);
    end
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
end
%% 
Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
MrT_inds = find(strcmp(FileName(:,1),'MrT'));

figure; hold on; 

for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if size(FM{daynum},1) == 3; cfp = {'b','g','r'}; else cfp = {'b','g'}; end
    for i = 1:size(FM{daynum},1)
        plot(daynum + x_off, FM{daynum}(i),'.','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off, [FL{daynum}(i) FU{daynum}(i)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);
%%
cfp = {'b','g','r'};
for i = 1:2
    %Mihili
    FRCmeans = cellfun(@(x) nanmean(x,1)',FRC,'UniformOutput',0);
    FRClengths = cell2mat(cellfun(@(x) length(x)>1,FRCmeans,'UniformOutput',0));
    allvals = vertcat(FRCmeans{Mihili_inds,i});
    avvals = mean(allvals); 
    [l_vals, u_vals] = boot_bounds(10000,@mean,allvals,2.5,97.5);
 
    plot([1 1]*(1.5+length(Mihili_inds)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
    plot(1.5+length(Mihili_inds),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
    
    %Mr T
    allvals = vertcat(FRCmeans{MrT_inds,i});
    avvals = mean(allvals); 
    [l_vals, u_vals] = boot_bounds(10000,@mean,allvals,2.5,97.5);

    plot([1 1]*(4.5+size(FileName,1)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
    plot(4.5+size(FileName,1),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
    
end
