units = load('C:\Users\limblab\Desktop\Mihili_bdf\Unit_info\Mihili_PMd_Units_07192013.mat');
Wcell = struct2cell(units);
shapes = cellfun(@(x) x(4:51),Wcell,'UniformOutput',0);
waves = vertcat(shapes{:});
normwave = (waves-repmat(min(waves,[],2),1,48))./(repmat(max(waves,[],2)-min(waves,[],2),1,48));
for i =1:size(waves,1)
    dur(i) = find(normwave(i,:)==max(normwave(i,:)),1,'first')-find(normwave(i,:)==min(normwave(i,:)),1,'last'); 
end
t = kmeans(dur,2);
%%
inrange = @(x,range) x>range(1) & x<range(2);
time_range = [alldays(1).tt(1,1) alldays(1).tt(end,1)];

trains = cell(length(alldays(1).PMd_units),1);
for i = 1:length(alldays(1).PMd_units)
    trains{i} = train2bins(alldays(1).PMd_units{i}(inrange(alldays(1).PMd_units{i},time_range)),.001);
    clc; fprintf('Neuron: %d/%d\n',i,length(alldays(1).PMd_units));
    
end

%%
%trainmat = vertcat(trains{:});
%     trainscat = cellfun(@(x) x(1:min(cellfun(@length,trains))),trains,'UniformOutput',0);
trainmat2 = cellfun(@(x) padarray(x,[0 (max(cellfun(@length,trains))-length(x))],0,'post'),trains,'UniformOutput',0);
%     trainmat = vertcat(trainscat{:});
trainarray = vertcat(trainmat2{:});
%%
C = combnk(1:size(trainarray,1),2);
corrs = zeros(size(C,1),201);
[CMAT,LAGMAT] = deal(NaN(size(trainarray,1)));
for i = 1:size(C,1)
    clc; fprintf('Combination %d/%d\n',i,size(C,1));
    [corrs(i,:),lags] = xcorr(trainarray(C(i,1),:),trainarray(C(i,2),:),100,'coeff');
    CMAT(C(i,1),C(i,2)) = max(corrs(i,:));
    LAGMAT(C(i,1),C(i,2)) = lags(find(corrs(i,:)==max(corrs(i,:)),1,'first'));
end