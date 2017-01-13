function[traces] = extract_move_traces(bdf,tt,col_start,col_end)

times = bdf.pos(:,1);
traces = cell(size(tt,1),1);

for i = 1:size(tt,1)    
    traces{i} = bdf.pos(times > tt(i,col_start) & times < tt(i,col_end),2:3);
end