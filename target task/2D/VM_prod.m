function [mu,kappa,fulldist] = VM_prod(xmu,xkappa,ymu,ykappa)

thets = 0:0.001:2*pi;

adist = circ_vmpdf(thets,xmu,xkappa);
bdist = circ_vmpdf(thets,ymu,ykappa);

fulldist = adist.*bdist;

fulldist = (fulldist./diff(thets(1:2)))./sum(fulldist);
mu = thets(fulldist==max(fulldist));

eq = @(k) exp(k*cos(thets-mu))./(2*pi*besseli(0,k));
diffeq = @(k) sum((eq(k)'-fulldist).^2);

kappa = fminbnd(diffeq,0.001,500);
