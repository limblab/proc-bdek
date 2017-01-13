%%
genfunc = @(muL,muP,sigL,sigP) (1./(sigL^2+sigP^2)).*(sigP^2.*(muL + randn(length(muL),1).*sigL)+ ...
                                          sigL^2.*(muP + randn(length(muP),1).*sigP));
ls = flipud(unique(alldays(2).tt(:,3)));
[stdres,dat_slope,KL,KP,KR,sigR,sigP,sigL] = deal(zeros(length(ls),1));
[doi,rs,g] = deal(cell(length(ls),1));
% figure; hold on; 
[~,mrs] = circ_polyfit(alldays(1).tt(:,2),alldays(1).tt(:,10));
motor_s = circ_std(mrs);
clrs = {'b','r','g'};
for i =  1:length(ls)
    is = find(alldays(2).tt(:,3)==ls(i));
    
    
    [doi{i},rs{i}] = circ_polyfit(alldays(2).tt(is,9),alldays(2).tt(is,10));
    stdres(i) = circ_std(rs{i});% - motor_s;
    dat_slope(i) = doi{i}(1);

%     KL(i) = kapres(i)./(dat_slope(i).^2 - dat_slope(i) + 1);

%       sigL(i) = stdres(i)/(dat_slope(i) + sqrt(dat_slope(i)*(1-dat_slope(i))));
      sigL(i) = sqrt(var(rs{i})/dat_slope(i));
      sigP(i) = sqrt(var(rs{i})/(1-dat_slope(i)));
      
%     KL(i) = ((dat_slope(i).^2 - dat_slope(i) + 1)./stdres(i))^2;
%     KP(i) = (1-dat_slope(i))*KL(i)./dat_slope(i);
%     
%     kf = zeros(length(is),1);
%     for j = 1:length(is)
%         [~,kf(j)] = VM_prod(alldays(2).tt(is(j),9),KL(i),circ_mean(alldays(2).tt(:,2)),KP(i));
%     end
      sigR(i) = std(rs{i});

    muLdata = alldays(2).tt(is,9);
    muPdata = repmat(circ_mean(alldays(2).tt(:,2)),size(muLdata,1),1);

    
%     plot(alldays(2).tt(is,9),g{i},'.','Color',clrs{i});
end
                      