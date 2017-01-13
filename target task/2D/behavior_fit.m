%%
genfunc = @(muL,muP,kL,kP) (1./(kL+kP)).*(kL.*(muL + vonmisrand_vec(length(muL),kL))+ ...
                                          kP.*(muP + vonmisrand_vec(length(muP),kP)));
ls = flipud(unique(alldays(2).tt(:,3)));
[kapres,dat_slope,KL,KP,KR,sigR,sigP,sigL,truesigP,truesigL,truesigR] = deal(zeros(length(ls),1));
[doi,rs,g] = deal(cell(length(ls),1));
% figure; hold on; 
[~,mrs] = circ_polyfit(alldays(1).tt(:,2),alldays(1).tt(:,10));
motor_s = circ_kappa(mrs);
clrs = {'b','r','g'};
for i =  1:length(ls)
    is = find(alldays(2).tt(:,3)==ls(i));
    
    [doi{i},rs{i}] = circ_polyfit(alldays(2).tt(is,9),alldays(2).tt(is,10));
%     kapres(i) = 1/(1/sqrt(circ_kappa(rs{i})) - 1/sqrt(motor_s))^2;
    kapres(i) = circ_kappa(rs{i});
    dat_slope(i) = doi{i}(1);

    KL(i) = dat_slope(i)*kapres(i);
    KP(i) = (1-dat_slope(i))*kapres(i);
    KR(i) = kapres(i);

    sigL(i) = 1/sqrt(KL(i));
    sigP(i) = 1/sqrt(KP(i));
    sigR(i) = 1/sqrt(kapres(i));

    muLdata = alldays(2).tt(is,9);
    muPdata = repmat(circ_mean(alldays(2).tt(:,2)),size(muLdata,1),1);
    
    g{i} = genfunc(muLdata,muPdata,KL(i),KP(i));
    
    truesigP(i) = 1/sqrt(circ_kappa(alldays(2).tt(is,2)));
    truesigL(i) = 1/sqrt(unique(alldays(2).tt(is,3)).*5);
    [~,kr] = VM_prod(0,1/(truesigL(i)^2),0,1/(truesigP(i)^2));
    truesigR(i) = 1/sqrt(kr);
    
%     plot(alldays(2).tt(is,9),g{i},'.','Color',clrs{i});
end
                      