lines = alldays(2).slices;
% k_fun = @(ws) (ws(1) + ws(2))/(ws(1)^2/k_l + ws(2)^2/k_p) * (ws(1)^2 + ws(2)^2 + 2*ws(1)*ws(2)*cos(mu_p-mu_l))^0.5;
options = optimset('Display','Off','Algorithm','active-set');
[ww, ww2, ww3] = deal(zeros(size(lines,1),2));
[mu, mu2,mu3,wwval,ww2val,ww3val] = deal(zeros(size(lines,1),1));
for i = 1:size(lines,1)
    kappas(i) = circ_kappa(lines(i,:));
    
    move = alldays(2).tt(i,10);
    mu_l = mod(alldays(2).tt(i,9),2*pi);
    k_l = kappas(i);
    mu_p = circ_mean(alldays(2).tt(1:i,2));
    k_p = circ_kappa(alldays(2).tt(1:i,2));
    mu_prev = [0; alldays(2).tt(2:end,10)];
    
   % mu_p = circ_mean(mu_p,mu_prev(i));
    
    mean_fun = @(ws) mu_l + atan2(sin(mu_p - mu_l), (ws(1)./ws(2)) + cos(mu_p - mu_l));

    k_fun = @(ws) -((ws(1) + ws(2))/(ws(1).^2./k_l + ws(2).^2./k_p) * (ws(1).^2 + ws(2).^2 + 2.*ws(1).*ws(2).*cos(mu_p-mu_l)).^0.5);
    b_fun = @(ws) abs(circ_dist((mu_l + atan2(sin(mu_p - mu_l), (ws(1)./ws(2)) + cos(mu_p - mu_l))),move));
    simple_fun = @(ws) abs(circ_dist(circ_mean(ws.*mu_l, (1-ws).*mu_p), move));
    
	%[ww(i,:),wwval(i)] = fminsearch(k_fun,[k_l k_p]);
    %[ww2(i,:),ww2val(i)]= fmincon(b_fun,[1; 1],[-1 0; 0 -1],[0;0],[],[],[-1000; -1000],[1000; 1000],[],options);
    %[ww3(i,:),ww3val(i)]= fmincon(simple_fun,[1; 1],[-1 0; 0 -1],[0;0],[],[],[-1000; -1000],[1000; 1000],[],options);
    [ww3(i,1),ww3val(i)]= fminbnd(simple_fun,-1,1); ww3(i,2) = 1-ww3(i,1);
    
    %mu(i) = mu_l + atan2(sin(mu_p - mu_l),(ww(1)/ww(2)) + cos(mu_p - mu_l));
    %mu2(i) = mod(mean_fun(ww2(i,:)),2*pi);
    mu3(i) = mod(circ_mean(ww3(i).*mu_l, (1-ww3(i)).*mu_p),2*pi);
end
