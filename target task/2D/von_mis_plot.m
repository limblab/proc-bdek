function[xs, dist] =  von_mis_plot(mu,kappa,do_plot)

mu = mu*pi/180;
xs = 0:0.01:2*pi;

dist = exp(kappa.*cos(mu-xs))./(2.*pi.*besseli(0,kappa));

if do_plot == 1
    figure; plot(xs,dist);
end