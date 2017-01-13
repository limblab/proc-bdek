VM = @(x,mu,k) exp(k.*cos(x-mu))./(2*pi*besseli(0,k));

dt= 0.01;
xs = 0:dt:2*pi;


mu1 = pi/2;
mu2 = 3*pi/2;

k1 = 2;
k2 = 2;

D1 = VM(xs,mu1,k1);
D2 = VM(xs,mu2,k2);

figure; hold on; 
plot(xs,D1,'b'); plot(xs,D2,'r'); 

muApp = mu1 + atan2(sin(mu2-mu1),k1/k2 + cos(mu2-mu1));
kApp = sqrt(k1^2 + k2^2 + 2*k1*k2*cos(mu2-mu1));

DA = VM(xs,muApp,kApp);

plot(xs,DA,'k'); drawnow;


%% Estimate distribution

PP = zeros(1000,length(xs));
for i = 1:1000
  m = mu1 + vonmisrand(k1);
  d = mu2 + vonmisrand(k2);
  p_m = VM(m,mu1,k1);
  p_d = VM(d,mu2,k2);
      
  for j = 1:length(xs)
    z = xs(j);

    p_mlz = VM(z,mu1,k1);
    p_dlz = VM(z,mu2,k1);
    
    PP(i,j) = p_mlz*p_dlz/(p_m*p_d);
  end
end



        
        
        
