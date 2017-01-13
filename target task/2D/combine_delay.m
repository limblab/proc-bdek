delayF = cell(length(likes),1);
for lik = 1:length(likes)
    
    firedelay = fire_trial{lik}(10:20);
    firedelay_cat = horzcat(firedelay{:});
    dpd_cat = repmat(DS{lik},1,length(firedelay));
    
    for trial = 1:size(dpd_cat,1)
        for dir = 1:length(spatial_cents)

            dirinds = find(abs(circ_dist(dpd_cat(trial,:),spatial_cents(dir)))<(space_size/2));
            delayF{lik}(trial,dir) = mean(firedelay_cat(trial,dirinds));
        end 
    end
end
%%
x = .5*cos(spatial_cents)' + .5;
[r,rfull,lfrLIK,lfrBIN,lfrIND] = deal(cell(length(likes),1));
%x = nanmean(delayF{1})';
for lik = 1:length(likes)
    L = size(DS{lik},1);
    maxlag = (length(x)-1)/2;
    
    likinds = find(alldays(4).tt(:,12)==likes(lik));
    lfrLIK{lik} = lfr(likinds);
    
    for i = 1:L

        ds = circ_dist(lfrLIK{lik}(i),spatial_cents);
       
        lfrIND{lik}(i,:) = find(abs(ds)==min(abs(ds)));
        lfrBIN{lik}(i,:) = spatial_cents(abs(ds)==min(abs(ds)));

        y = delayF{lik}(i,:);

        scaled = (y - nanmean(y)) / nanstd(y);
        scaled(isnan(y)) = 0;

        cor = xcorr(x,scaled,maxlag);
        r{lik}(i) = spatial_cents(cor==max(cor));
        rfull{lik}(i,:) = cor;
    end
end
%shifts = linspace(-(length(x)-1)/2,(length(x)-1)/2,length(x));
% for i = 1:length(x)
%     
%     xshift = circshift(x,shifts(i));
%     
%     rc(i) = circ_corrcc(y,xshift);
% end
%-% re-center



