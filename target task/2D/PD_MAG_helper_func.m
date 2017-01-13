function[PDMAG,b] = PD_MAG_helper_func(xs_ys_spds)

xsys = xs_ys_spds(:,1:2);

xs = xsys(:,1);
ys = xsys(:,2);

if size(xs_ys_spds,2)==3
    spds = xs_ys_spds(:,3);
    b = glmfit([cos(xs) sin(xs) spds spds.*cos(xs) spds.*sin(xs)],ys,'poisson');
    PD = atan2(b(3),b(2));
    MAG = sqrt(b(2).^2+b(3).^2);
elseif size(xs_ys_spds,2)==2
    b = glmfit([cos(xs) sin(xs)],ys,'poisson');
    PD = atan2(b(3),b(2));
    MAG = sqrt(b(2).^2+b(3).^2);
else
    error('input must have either 2 or 3 columns'); 
end

PDMAG = [PD MAG];
  
end