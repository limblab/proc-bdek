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
  
FileName = {'Mihili','07112013',  2        ;...
            'Mihili','07152013',  2        ;...
            'Mihili','07192013',  2        ;...
            'Mihili','08062013',  2        ;...
            'Mihili','08122013',  2        ;...
            'Mihili','08152013',  2        ;...
            'MrT'   ,'05042013',  2        ;...
            'MrT'   ,'05052013',  2        ;...
            'MrT'   ,'05062013',  2        };
  
% FileName = {'Mihili','07192013',  2};


        
[RR2full,RR2full2,PR2,PAS] = deal(cell(size(FileName,1),1)); %initialize
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
     
    if exist(sprintf('PMd_PD90_%s_%s.mat',FileName{daynum,1},FileName{daynum,2}),'file')
        load(sprintf('PMd_PD90_%s_%s.mat',FileName{daynum,1},FileName{daynum,2}));
    else
        tuning_types_PD;
        save(sprintf('PMd_PD90_%s_%s',FileName{daynum,1},FileName{daynum,2}),'best_PDS');
    end
    
    
    pretarget_plot_PDOD; close all;
    VA = voif_all{1};
    
    GLM_pseudoR2;
    
    RR2full{daynum} = 100*(pseudoR2.full-pseudoR2.part)./pseudoR2.full;
    RR2full2{daynum} = 100*(pseudoR2.full2-pseudoR2.part)./pseudoR2.full2;
    PR2{daynum} = pseudoR2;
    PAS{daynum} = pass;
    % Fill variables
    clearvars -except RR2full RR2full2 PR2 FileName PAS;
    
end
