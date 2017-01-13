blck = find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks));
trlblck = find(DAT.TrialBlocks{blck}.Trial == trl2plot);
thets = linspace(0,2*pi,100); figure; hold on; 
r1 = (str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)));
r2 = (str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)));
patch([r1*cos(thets) r2*cos(fliplr(thets))],[r1*sin(thets) r2*sin(fliplr(thets))],[.25 .25 1],'EdgeColor',[.25 .25 1]); 
axis square; axis equal; axis off;
if size(DAT.TrialBlocks{blck}.VisualCueAngles,2)>1
    plot([r1*cos(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:));r2*cos(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:))],...
         [r1*sin(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:));r2*sin(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:))],'r','LineWidth',5);
else
    a1 = (DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,' ')-1)));
    a2 = (DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,' ')-1)));
    patch([r1*cos(a1:pi/50:a2) r2*cos(fliplr(a1:pi/50:a2))],[r1*sin(a1:pi/50:a2) r2*sin(fliplr(a1:pi/50:a2))],'r','EdgeColor','r'); 
end


% plot(downsample(DAT.kinematics.XPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),downsample(DAT.kinematics.YPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),'k','LineWidth',5)

plot(downsample(DAT.kinematics.XPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),downsample(DAT.kinematics.YPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),'k','LineWidth',5)


plot(downsample(DAT.kinematics.XPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),downsample(DAT.kinematics.YPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks))}.Trial == trl2plot))),10),'k','LineWidth',5)
