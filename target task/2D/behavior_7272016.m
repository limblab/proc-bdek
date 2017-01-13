C = mod(round(alldays(1).tt(:,2)/(pi/4)), 8)+1;
[Kbd,Mbd] = deal(zeros(length(unique(C))+1,1));
for tL = 1:length(unique(C)) 
    Kbd(tL) = circ_kappa(circ_dist(alldays(1).tt(C==tL,10),alldays(1).tt(C==tL,2)));
    Mbd(tL) = circ_mean(circ_dist(alldays(1).tt(C==tL,10),alldays(1).tt(C==tL,2)));
end
Kbd(end) = Kbd(1);
Mbd(end) = Mbd(1);

Kbydir_func = @(x) interp1(0:pi/4:2*pi,Kbd,x);
Mbydir_func = @(x) interp1(0:pi/4:2*pi,Mbd,x);

%%
[v_mean,v_kap,p_mean,p_kap] = deal(zeros(size(alldays(2).tt,1),1));
for i = 1:size(alldays(2).tt,1)
    clc; fprintf('%d/%d\n',i,size(alldays(2).tt,1));
    [~,~,~,mdis] = boot_bounds(1000,@circ_mean,alldays(2).slices(i,:)',2.5,97.5);
    
    v_mean(i) = circ_mean(mdis);
    v_kap(i) = circ_kappa(mdis);
    
    if i > 1
        [~,~,~,pdis] = boot_bounds(1000,@circ_mean,alldays(2).tt(1:i,2),2.5,97.5);
        p_mean(i) = circ_mean(pdis); 
        p_kap(i) = circ_kappa(alldays(2).tt(1:i,2));
    else
        p_mean(i) = alldays(2).tt(i,2);
        p_kap(i) = 0.01;
    end
end

%%
[pofr,pofp,calc_pkap,est_postkap] = deal(zeros(size(alldays(2).tt,1),1));
for i = 1:size(alldays(2).tt,1)

    actreach = alldays(2).tt(i,10);
    actprior = circ_mean(alldays(2).tt(1:i,2));
    rloc = v_mean(i) + Mbydir_func(v_mean(i));
    rkap = 1/(1/v_kap(i) + 1/Kbydir_func(v_mean(i)));

    ploc = circ_mean(alldays(2).tt(1:i,2)) + Mbydir_func(v_mean(i));
    pkap = p_kap(i);
    
    pofr(i) = sum(circ_vmpdf((actreach-0.0349):0.0001:(actreach+0.0349),rloc,rkap))*0.0001;
    pofp(i) = sum(circ_vmpdf((actreach-0.0349):0.0001:(actreach+0.0349),ploc,pkap))*0.0001;
    
    
    calc_pkap(i) = rkap*circ_dist(actreach,rloc)./circ_dist(ploc,actreach);
    est_postkap(i) = calc_pkap(i) + rkap;
end
    


    