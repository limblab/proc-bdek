function plot_cross(x_mean,x_LH,y_mean,y_LH,color,width)

plot([x_LH(1),x_LH(2)],y_mean*[1 1],'-','Color',color,'LineWidth',width); hold on; 
plot(x_mean*[1 1],[y_LH(1),y_LH(2)],'-','Color',color,'LineWidth',width); 
end