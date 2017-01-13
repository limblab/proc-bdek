%%
figure; hold on;

for i = 1:length(a)

    theta = b(i);
    count = a(i);

    inout = [500 500+.5*count];

    p = ang2poly(theta,0.04,inout(1),inout(2));
    patch(p(:,1),p(:,2),p(:,3),'EdgeColor','k','FaceColor',[0.7 0.5 0.5]);
    axis([-1500 1500 -1500 1500]);

end