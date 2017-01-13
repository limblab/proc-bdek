INPUT = 'average';
%INPUT = 'specified'; %Must define XINPUT

thets = 0.01:0.01:2*pi;

pds = linspace(0,2*pi,100); pds = pds(1:(end-1));
bases = zeros(length(pds),length(thets));
for i = 1:length(pds)
    bases(i,:) = 0.5 + 0.5*cos(thets-pds(i));
end

normal_weights = repmat(bases(:,end/2),1,length(thets));
recruit = sum(normal_weights.*bases);

%%
tune_change = @(x,basemod) ...
    sum(repmat(x(:,end/2)*basemod(2) + basemod(1),1,length(thets)).*...
                                            (x*basemod(2) + basemod(1)));
                                                                             
%%     
fit_TC = zeros(length(av_gain),628);
if strcmp(INPUT,'average')
    
    colforlike = {'b','g','r'};
    figure; hold on; plot(thets,0.5*cos(thets-pi)+0.5,'k');
    p = zeros(length(av_gain),2);
    for i = 1:length(av_gain)  

        gain_recruit = recruit.*av_gain{i};

        p(i,:) = fminsearch(@(bm) sum((tune_change(bases,bm) - gain_recruit).^2),[0,1]); 
        
        fit_TC(i,:) = p(i,2)*0.5*cos(thets-pi)+0.5+p(i,1);
        
        plot(thets,fit_TC(i,:),colforlike{i});
    end
    
else
    
    gain_recruit = recruit.*XINPUT;
    
    p = fminsearch(@(bm) sum((tune_change(bases,bm) - gain_recruit).^2),[0,1]); 
end
    
    
    
    