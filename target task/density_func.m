function[max_dens_pos,cont_dens,x_range] = density_func(positions)

%% Set up filter
width = 100;
N = 5001;

densfilt = fspecial('gaussian',[N 1],width);

dx = 0.01;
positions = round(positions./dx);

x_range = -100/dx:100/dx;

bins = hist(positions,x_range);

%% Convolve Gaussian
cont_dens = conv(bins,densfilt);

cont_dens(1:(N-1)/2) = [];
cont_dens(end-(N-1)/2+1:end)=[];

max_dens_pos = dx*x_range(cont_dens==max(cont_dens));

end