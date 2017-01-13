function[tune_curve] = dens_pref_dir(alldays,adj_tt,allday_ind,unit,t1,t2)

[ignore,Lout,Hout] = raster_UNT_circ(alldays(allday_ind).bdfM,alldays(allday_ind).PMd_units{unit},adj_tt{allday_ind},t1,t2);

raster = Lout;

totspikes = nansum(raster,2);
directions = adj_tt{allday_ind}(:,10);

kappa = 10;
xs = 0:0.01:2*pi;
tune_curve = zeros(1,length(xs));

for i = 1:length(directions)
    mu = directions(i);
    dist = exp(kappa.*cos(mu-xs))./(2.*pi.*besseli(0,kappa));
    
    weighted = dist.*totspikes(i);
    
    tune_curve = tune_curve + weighted; 
    
end

    
    