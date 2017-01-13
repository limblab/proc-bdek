function [rotated_traces,tracex,tracey,paddedx,paddedy] = average_trace_path(traces,inds,doplot)

if nargin == 1
    inds{1} = 1:length(traces);
    doplot = 0;
elseif nargin == 2
    doplot = 0;
end
    
rotated_traces = cell(length(traces),1);

maxlength = max(cellfun(@length,traces));

tracex = zeros(length(traces),maxlength);
tracey = tracex;
paddedx = tracex;
paddedy = tracey;
paddedt = zeros(length(traces),500);

for i = 1:length(traces)
    
    trace = traces{i};
    startpoint = trace(1,:);
    shiftedtrace = trace - repmat(startpoint,size(trace,1),1);
    
    orig_angle = atan2(shiftedtrace(end,2),shiftedtrace(end,1));
    
    R = [cos(-orig_angle) -sin(-orig_angle); sin(-orig_angle), cos(-orig_angle)];
    
    rotated_traces{i} = (R*shiftedtrace')';
    
    tracex(i,1:length(trace)) = rotated_traces{i}(:,1)';
    tracey(i,1:length(trace)) = rotated_traces{i}(:,2)';
    
    paddedt(i,:) = interp1(1:length(trace),tracey(i,1:length(trace)),linspace(tracex(i,1),tracex(i,length(trace)),500));
    
    paddedx(i,:) = interp1(1:length(trace),rotated_traces{i}(:,1)',(1:maxlength)*(length(trace)/maxlength));
    paddedy(i,:) = interp1(1:length(trace),rotated_traces{i}(:,2)',(1:maxlength)*(length(trace)/maxlength));
end

if doplot ~= 0
    figure; hold on;
    plotcolors = {'k','b','r'};
    [avx,avy,avt] = deal(cell(length(inds),1));
    for i = 1:length(inds)
        curinds = inds{i};
        
        paddedy(:,end) = 0;
        paddedy(:,1) = 0;
        paddedx(:,1) = 0;
        
        avx{i}.M = nanmean(paddedx(curinds,:),1);
        avy{i}.M = nanmean(paddedy(curinds,:),1);
        
%         avt{i}.M = nanmean(paddedt(curinds,:),1);
%         [avt{i}.L, avt{i}.H] = boot_bounds(1000,@nanmean,paddedt(curinds,:),2.5,97.5);
        
%         [avx{i}.L, avx{i}.H] = boot_bounds(1000,@nanmean,paddedx(curinds,:),2.5,97.5);
        [avy{i}.L, avy{i}.H] = boot_bounds(1000,@nanmean,paddedy(curinds,:),2.5,97.5);
    

        plot(1:length(avy{i}.M),avy{i}.M','Color',plotcolors{i});
        patch([1:length(avy{i}.M) fliplr(1:length(avy{i}.M))],[avy{i}.H' fliplr(avy{i}.L')],plotcolors{i},'FaceAlpha',1);
    end
end
        
    


