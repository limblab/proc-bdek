aligned_slices = alldayso(2).slices - repmat(alldayso(2).tt(:,10),1,5);

hinds = alldayso(2).tt(:,3) == 1;

lefters = sum(circ_dist(aligned_slices,0) < 0,2) > 3;
righters = sum(circ_dist(aligned_slices,0) > 0,2) > 3;

left_inds = find(lefters + hinds == 2);
right_inds = find(righters + hinds == 2);

allleft = reshape(aligned_slices(left_inds,:),1,numel(aligned_slices(left_inds,:)));
allright = reshape(aligned_slices(right_inds,:),1,numel(aligned_slices(right_inds,:)));

figure; hold on;
plot(repmat(allright,2,1),[zeros(1,length(allright)) ; ones(1,length(allright))],'Color',[1 0 0]);
plot(repmat(allleft,2,1),[zeros(1,length(allleft)) ; ones(1,length(allleft))],'Color',[.5 .5 .5]);