function[M] = combine_Msessions(MALL)

max_blocks = max(cellfun(@(x) length(x.tt),MALL));

epsizes = cell2mat(cellfun(@(x) x.LENS,MALL,'Uni',0)');
minepsize = min(epsizes,[],1);
removetime = cell(length(MALL),1);
for i = 1:length(MALL)
    toolarge = find(epsizes(i,:)-minepsize ~= 0);
    if ~isempty(toolarge)
        removetime{i} = MALL{i}.CLENS(toolarge+1);
        for j = 1:max_blocks
            if j<=length(MALL{i}.tt)
                MALL{i}.pmd{j}.d(:,removetime{i},:) = [];
                MALL{i}.m1{j}.d(:,removetime{i},:) = [];
            end
        end
    end
end

M = [];

for i = 1:max_blocks %% i = TRIAL BLOCK
    curPMD = [];
    curM1 = [];
    curtargs = [];
    TT = cell(length(MALL),1);
    for j = 1:length(MALL) %% j = SESSION NUMBER
        if i<=length(MALL{j}.tt)
           TT{j,1} = MALL{j}.tt{i};

           curPMD = cat(3,curPMD,MALL{j}.pmd{i}.d);
           curM1 = cat(3,curM1,MALL{j}.m1{i}.d);
           curtargs = cat(1,curtargs,MALL{j}.targids{i});
        end
    end
    

    M.tt{i} = cell2mat(TT);
    M.pmd{i}.d = curPMD;
    M.m1{i}.d = curM1;
    M.targids{i} = curtargs;
end
M.CLENS = [1 cumsum(minepsize)];
M.LENS = minepsize;
M.meta = MALL{1}.meta;


