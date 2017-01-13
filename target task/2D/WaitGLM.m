function WaitGLM(BDF,comp_tt,trains)

for i = 1:length(comp_tt)
    
    %% Get indices 
    targ_on = find(BDF.pos(:,1)>comp_tt(i,5),1,'first');
    go_cue = find(BDF.pos(:,1)>comp_tt(i,6),1,'first');
    trial_end = find(BDF.pos(:,1)>comp_tt(i,7),1,'first');
    
    %% Kinematics
    speed = sum(BDF.vel(go_cue:trial_end,2:3).^2,2);
        maxspeed = max(speed);
        maxind = find(speed==maxspeed)+go_cue-1;
        maxdir = atan2(BDF.vel(maxind,3),BDF.vel(maxind,2));
    cos_ang = cos(maxdir);
    sin_ang = sin(maxdir);
    vx = maxspeed*cos_ang;
    vy = maxspeed*sin_ang;
    guessx = BDF.pos(trial_end,2);
    guessy = BDF.pos(trial_end,3);
    
    %% Trial params
    
    y(i) 
end