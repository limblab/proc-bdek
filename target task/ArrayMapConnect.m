function[fig_handle] = ArrayMapConnect(monkey,implant,unit_ids,from_to_neurs) 
%% Get Mapfile and assign metrics
map = ArrayMap(monkey,implant); % Get mapfile
% Function to parse unit information
info_parse = @(id) [floor(id) round(10*(id - floor(id)))];

%% Populate map with metrics and flip to prepare for plotting
%Initialize arrays
[from_map,flipped_map] = deal(cell(size(map,1),size(map,2),max(round(10*(unit_ids - floor(unit_ids))))));
    %deal(zeros(size(map,1),size(map,2),max(round(10*(unit_ids - floor(unit_ids))))));
[allmap,flipped_allmap] = deal(zeros(size(map,1),size(map,2),max(round(10*(unit_ids - floor(unit_ids))))));

orig_from2neurs = from_to_neurs;
    
unit_ids = floor(unit_ids);
diffids = diff(unit_ids);

samechan = find(diffids==0) + 1;

unit_ids = unit_ids + 0.1;
unit_ids(samechan) = unit_ids(samechan) + 0.1;

all2chans = floor(from_to_neurs);
nonzerochans = all2chans(all2chans~=0);
uniquenonzerochans = unique(nonzerochans);
counts  = histc(nonzerochans,uniquenonzerochans);
dups = uniquenonzerochans(counts > 1);
for i = 1:length(dups)
    dupchan = dups(i);
    
    chansofdup = find(all2chans == dupchan);
    valsofdup = from_to_neurs(chansofdup);
    
    numuns = unique(valsofdup);
    
    for qq = 1:length(valsofdup)
        all2chans(chansofdup(qq)) = all2chans(chansofdup(qq)) + .1*(find(valsofdup(qq)==numuns)-1);
    end 
end
    
from_to_neurs(from_to_neurs ~= 0) = all2chans(from_to_neurs ~= 0) + 0.1;


for i = 1:length(unit_ids) %Loop through units and place metrics on map
    
    chan_unit = info_parse(unit_ids(i)); % Parse channel/unit information
    channel = chan_unit(1); unitnum = chan_unit(2);
    
    [locx,locy] = find(map==channel); % Find channel location on array
    
    allmap(locx,locy,unitnum) = 1;
    
    if sum(from_to_neurs(i,:)) > 0
        for jj = 1:size(from_to_neurs(i,:),2)
            
            if from_to_neurs(i,jj) > 0
            
            tchan_unit = info_parse(from_to_neurs(i,jj)); % Parse channel/unit information
            tchannel = tchan_unit(1); tunit = tchan_unit(2); 
    
            [tlocx,tlocy] = find(map==tchannel); % Find channel location on array
        
            from_map{locx,locy,unitnum} = [tlocx-.2*(tunit-1),tlocy+.2*(tunit-1)];
            end
        end
    end
end

%% Do plotting

for i = 1:size(from_map,3) % Flip map to prepare for plotting
    flipped_map(:,:,i) = flipud(from_map(:,:,i))';
    flipped_allmap(:,:,i) = flipud(allmap(:,:,i))';
end

fig_handle = figure; hold on; 

plot(1,10,'w.','MarkerSize',0.00001); 
plot(10,1,'w.','MarkerSize',0.00001); 
plot(1,1,'w.','MarkerSize',0.00001); 
plot(10,10,'w.','MarkerSize',0.00001); 
axis square

for k = 1:size(from_map,3) % Loop through multi-unit channels
    for j = 1:size(from_map,2) % First dimension
        for i = 1:size(from_map,1) % Second dimension
            if flipped_allmap(i,j,k) > 0
                plot(i+.2*(k-1),j+.2*(k-1),'k.');
            end
            if ~isempty(flipped_map{i,j,k})
                                
                arrow([i+.2*(k-1),j+.2*(k-1)],[flipped_map{i,j,k}(2), 11-flipped_map{i,j,k}(1)],'Width',3);
                plot(i+.2*(k-1),j+.2*(k-1),'bo'); plot(flipped_map{i,j,k}(2), 11-flipped_map{i,j,k}(1),'bo')
%             elseif flipped_map 
%                 plot(i,j,'r.');
            end
        end
    end
end

axis off

