recruit_move = 0.5*cos(linspace(-pi,pi,1000)) + 1;

point_rec = zeros(1,1000); point_rec(500) = 1;
dist_rec  = circ_vmpdf(linspace(-pi,pi,1000),0,50);
pulse_rec = zeros(1,1000); pulse_rec(400:600) = 1;


c1 = conv(recruit_move,point_rec,'same');
c1n = c1./max(c1);

c2 = conv(recruit_move,dist_rec,'same');
c2n = c2./max(c2);
