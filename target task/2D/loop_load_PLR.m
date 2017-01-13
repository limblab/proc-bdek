%% Files to Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
BRAIN_AREA = 'M1';
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
            
            'Mihili','07122013',  2,  [1, 2]    ;...
    
            'Mihili','08012013',  2,  [1, 2]    ;...
            'Mihili','08012013',  3,  [1, 2]    ;...
            'Mihili','08012013',  4,  [1, 2]    ;...
            
            'Mihili','08222013',  2,  [1, 2]    ;...
            'Mihili','08222013',  3,  [1, 2]    ;...
            
            'Mihili','09042013',  2,  [1, 2]    ;...
            
            'Mihili','09052013',  2,  [1, 2]    ;...
            
            'Mihili','09062013',  2,  [1, 2]    ;...
            
            'Mihili','09262013',  2,  [1, 2]    ;...

            'Mihili','10022013',  2,  [1, 2]    ;...


            'Mihili','10072013',  2,  [1, 2]    ;...

};

% FileName = {'Mihili','08062013',  2,  [1, 2]    };...
% FileName = {'MrT'   ,'05042013',  2,  [1, 2]    ;...
%             'MrT'   ,'05052013',  2,  [1, 2]    ;...
%             'MrT'   ,'05062013',  2,  [1, 3]    ;...
%             'MrT'   ,'07082013',  2,  [1, 2]    ...
% };        

session_limit = 10e10;

%% Behavior

[split_indices,PLRb,KRATS] = deal(cell(size(FileName,1),1)); %initialize
[PLR,LI,PLRb]=deal(cell(100,1));
[std_priors,std_likes,std_posts] = deal(nan(10000,2));
counter = 1;
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
    
    llist = flipud(unique(alldays(2).tt(:,3))); 
    alldays(2).tt(~ismember(alldays(2).tt(:,3),llist(FileName{daynum,4})),:) = [];

    split_indices{daynum} = 1:session_limit:size(alldays(2).tt,1);
    split_indices{daynum}(end) = size(alldays(2).tt,1);
    if(length(split_indices{daynum})==1); split_indices{daynum} = [1 size(alldays(2).tt,1)]; end
    
    [~,~,fitkrats,fitpris] = behavior_fit_circ(alldays(2).tt(:,[2 9 10 3]));
    
    ALLDAYS = alldays;
    if strcmp(FileName{daynum,1},'MrT'); numslcs = 10; else numslcs = 5; end
    for section = 1:(length(split_indices{daynum})-1)

        indis = split_indices{daynum}(section):split_indices{daynum}(section+1);
        alldays(2).tt = ALLDAYS(2).tt(indis,:);
                
        liklist = flipud(unique(alldays(2).tt(:,3)));
        for lks = 1:length(liklist)
            LI{counter}{lks} = find(alldays(2).tt(:,3)==liklist(lks));
        end
        [PLR{counter},PLRb{counter}] = behavior_fit_circ(alldays(2).tt(:,[2 9 10 3]),fitkrats,fitpris);

        std_priors(counter,:) = 1./sqrt([PLR{counter}{1}(1) PLR{counter}{2}(1)]);
        std_likes(counter,:)  = 1./sqrt([PLR{counter}{1}(2) PLR{counter}{2}(2)]);
        std_posts(counter,:)  = 1./sqrt([PLR{counter}{1}(3) PLR{counter}{2}(3)]);
        
        KRATS{counter} = fitkrats;
        counter = counter + 1;
    end
    clc; fprintf('%d/%d\n',daynum,size(FileName,1)); 

    clearvars -except FL FU FM FileName RB VA BPDS BRAIN_AREA G counter PLR ...
        split_indices std_priors std_likes std_posts LI PLRb session_limit loop_ranges KRATS AV EF NN;
    
end
PLR(cellfun(@isempty,PLR)) = [];
PLRb(cellfun(@isempty,PLRb)) = [];
LI(cellfun(@isempty,LI)) = [];
std_priors(isnan(sum(std_priors,2)),:) = [];
std_likes(isnan(sum(std_likes,2)),:) = [];
std_posts(isnan(sum(std_posts,2)),:) = [];

G = 1:length(VA);

dRES_prior = diff(std_priors,[],2);
dRES_likes = diff(std_likes,[],2);
dRES_posts = diff(std_posts,[],2);

dRES_bnd_prior = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,1))-1./sqrt(x{1}(:,1)),[2.5 97.5]),PLRb,'UniformOutput',0));
dRES_bnd_likes = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,2))-1./sqrt(x{1}(:,2)),[2.5 97.5]),PLRb,'UniformOutput',0));
dRES_bnd_posts = cell2mat(cellfun(@(x) prctile(1./sqrt(x{2}(:,3))-1./sqrt(x{1}(:,3)),[2.5 97.5]),PLRb,'UniformOutput',0));

dRES = dRES_posts;
dRES_bnd = dRES_bnd_posts;
