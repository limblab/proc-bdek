%%
% Set up list of all days for inclusion
% FileName = {'Mihili','07112013',  [1 2]        ;...
%             'Mihili','07152013', [1,2,3]     ;...
%             'Mihili','07192013', [1,2,3,4]   ;...
%             'Mihili','08062013', [1,2,3]     ;...
%             'Mihili','08122013', [1,2,3]     ;...
%             'Mihili','08152013',  [1 2]        ;...
%             'MrT'   ,'05042013',  [1 2]        ;...
%             'MrT'   ,'05052013',  [1 2]        ;...
%             'MrT'   ,'05062013',  [1 2]        };
%    
% FileName = {'Mihili','05062014',  [1 2]        ;...
%             'Mihili','06122014',  [1 2]        };
        
FileName = {'Mihili','07192013',  [1 2]};
  
        
        
[FL, FU, FM, FRC] = deal(cell(size(FileName,1),1)); %initialize
for daynum = 1:size(FileName,1)
    
    % Load file
    [alldays, priors, monkey, MO, DA, YE] = load_processed(FileName{daynum,1},FileName{daynum,2}); 
    priors{1}.val = 1; priors{2}.val = 2; priors{3}.val = 3; % % set arbitrary prior values
    
    alldays(2).tt = vertcat(alldays(FileName{daynum,3}).tt); %Concatenate to include all prior blocks
    alldays(2).slices = alldays(2).tt(:,1:5); 
    
    [~,GOOD_list] = ISI_exam(alldays(2).PMd_units,1.7,0.1);
    rawcount_script; % Obtain spike counts, etc. 
    spike_count_raw{1}{1} = horzcat(spike_count_raw{1}{:}); spike_count_raw{1}(2:end) = [];
    Fano_scriptBIN; close; close;
   
    % Fill variables
    FL{daynum} = BF(:,1);
    FU{daynum} = BF(:,2);
    FM{daynum} = BF(:,3);
    for i = 1:length(LHvals)
        FRC{daynum,i} = FANOS{i};
    end
    fprintf('%d/%d\n',daynum,size(FileName,1)); 
end

%% 
Mihili_inds = find(strcmp(FileName(:,1),'Mihili'));
MrT_inds = find(strcmp(FileName(:,1),'MrT'));

figure; hold on; 

for daynum = 1:size(FileName,1)
    if ismember(daynum,Mihili_inds); x_off = 0; else x_off = 3; end
    if size(FM{daynum},1) == 4; cfp = {'k','b','g','r'}; else cfp = {'k','b','g'}; end
    for i = 1:size(FM{daynum},1)
        plot(daynum + x_off, FM{daynum}(i),'.','Color',cfp{i},'MarkerSize',10);
        plot(daynum*[1 1]+x_off, [FL{daynum}(i) FU{daynum}(i)],'Color',cfp{i},'LineWidth',1);
    end
end
xlim([-2 daynum+3+3]);

cfp = {'k','b','g','r'};
for i = 1:3
    %Mihili
    thesevals = cellfun(@(x) reshape(x,[],1), FRC(Mihili_inds,i),'UniformOutput',0);
    combvals = vertcat(thesevals{:});
    allvals = horzcat(combvals{:});
    avvals = nanmean(allvals); 
    [l_vals, u_vals] = boot_bounds(10000,@nanmean,allvals,2.5,97.5);
 
    plot([1 1]*(1.5+length(Mihili_inds)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
    plot(1.5+length(Mihili_inds),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
    
%     %Mr T
%     thesevals = cellfun(@(x) reshape(x,[],1), FRC(MrT_inds,i),'UniformOutput',0);
%     combvals = vertcat(thesevals{:});
%     allvals = horzcat(combvals{:});
%     avvals = nanmean(allvals); 
%     [l_vals, u_vals] = boot_bounds(10000,@nanmean,allvals,2.5,97.5);
% 
%     plot([1 1]*(4.5+size(FileName,1)),[l_vals, u_vals],'Color',cfp{i},'LineWidth',2);
%     plot(4.5+size(FileName,1),avvals,'o','Color',cfp{i},'MarkerSize',5,'MarkerFaceColor','w');
    
end
