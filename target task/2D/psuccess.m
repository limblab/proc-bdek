%% Create directional bins
prediction_day_indices = 2;
dbins = [0 pi];
%dbins = 0:pi/16:pi;
Psuccess = @(e,t) mean(abs(circ_dist(e,t)) < (15/180*pi));
PofS = @(adays,inds) Psuccess(adays.tt(inds,2),adays.tt(inds,10));

LHvals = flipud(unique(alldays(prediction_day_indices).tt(:,3)));
centdirs = alldays(prediction_day_indices).tt(:,9);
movedirs = alldays(prediction_day_indices).tt(:,10);
colp = {'b','g','r'};
    
[cbininds] = deal(cell(length(LHvals),length(dbins)-1));
[psuc,numt] = deal(zeros(length(LHvals),length(dbins)-1));

dfromp = abs(circ_dist(centdirs,circ_mean(alldays(prediction_day_indices).tt(:,2))));
for likei = 1:length(LHvals)
    likecondinds = find(alldays(prediction_day_indices).tt(:,3)==LHvals(likei));
    for i = 1:(length(dbins)-1) 
        centinds = find(dfromp > dbins(i) & dfromp < dbins(i+1));
        cbininds{likei,i} = centinds(ismember(centinds,likecondinds)); 
    end
    
    psuc(likei,:) = cellfun(@(is) PofS(alldays(prediction_day_indices),is), cbininds(likei,:));
    numt(likei,:) = cellfun(@length, cbininds(likei,:));
    
end

%%
% TOI = [cbininds{3,1}; cbininds{3,2}; cbininds{2,3}; cbininds{2,4}];
TOI = [vertcat(cbininds{3,1:2}); vertcat(cbininds{2,3:end})];
alldays = ALLDAYS;
alldays(prediction_day_indices).tt = ALLDAYS(prediction_day_indices).tt(TOI,:);

prediction_day_indices = 2;
pretarget_plot;
%Activity_script2;
alldays = ALLDAYS;
