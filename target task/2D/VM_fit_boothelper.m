function[outmat] = VM_fit_boothelper(fires_dfp)

spatial_cents = -pi:pi/16:pi;
space_size = diff(spatial_cents(1:2));

fires = fires_dfp(:,1:(end/2));
dfp = fires_dfp(:,(end/2+1):end);

fires_vec = reshape(fires,[],1);
dfp_vec = reshape(dfp,[],1);

FIRE = nan(length(spatial_cents),1);
for dir = 1:length(spatial_cents)

    dirinds = abs(circ_dist(dfp_vec,spatial_cents(dir)))<(space_size/2);
    
    FIRE(dir,:) = nanmean(fires_vec(dirinds),1);
end 

%[pfits,~,curve,base,gain,fwhm] = VM_fit2(spatial_cents',FIRE,[-pi:0.01:pi]);
%[pfits,~,curve,base,gain,fwhm] = VM_fit_constrained(spatial_cents',FIRE,[-pi:0.01:pi]);
[pfits,~,curve,base,gain,fwhm] = VM_fit3(spatial_cents',FIRE,[-pi:0.01:pi]);

%outmat = [curve pfits base gain fwhm];
outmat = [curve pfits 0 base gain fwhm];

end