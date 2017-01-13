function [PDS,PD_range,MAG_average] = get_PDS_fulldat(fulldat,allowed_range,varargin)

if nargin == 1
    allowed_range = pi/2;
end

[PD_range,PD_ok,MAG_average,MAG_ok] = deal(zeros(size(fulldat.PDS{end}{end},1),length(fulldat.PDS{end})));
MAG_ok(:) = 1;
for i = 1:length(fulldat.PDS{end})
    
    
    pds = fulldat.PDS{end}{i};
    pddist = fulldat.PDS_dist{end}{i};
    magdist = fulldat.MAGS_dist{end}{i};
    pdbounds = fulldat.PDS_bounds{end}{i};
    
    PD_range(:,i) = pdbounds(:,2)-pdbounds(:,1);
    
    for j = 1:size(pddist,1)
        PD_ok(j,i) = 1-circ_mtest(pddist(j,:),0);
    end

    MAG_ok(max(magdist,[],2) > 500,i) = 0;
    MAG_ok(mean(magdist,2) < 0.0001,i) = 0;
    MAG_ok(mode(magdist,2) < 0.0001,i) = 0;

    MAG_average(:,i) = mean(magdist,2);
end
PD_ok(PD_range > allowed_range) = NaN;
% PD_ok(PD_range > pi) = NaN;
PD_ok(PD_ok==0) = NaN;
MAG_ok(MAG_ok==0) = NaN; 

PD_range = PD_range.*PD_ok.*MAG_ok;
MAG_average = MAG_average.*MAG_ok.*PD_ok;

PDS = horzcat(fulldat.PDS{end}{:}).*MAG_ok.*PD_ok;