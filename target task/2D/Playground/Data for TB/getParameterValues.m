function[Params] = getParameterValues(Directory,Monkey,Date,Task)
%% find meta info
switch Task
    case 'UNT2D'
        tname = 'UT2D';
    case 'DC'
        tname = 'DC';
    case '
end

GeneralParams = {'X Offset','Y Offset', 'Cursor Rotation','Reward Pulse','Reward Jackpot',...
                 'Reward Jackpot Chance','Curl Angle','Static X Force','Static Y Force'};           
GeneralParams = cellfun(@(x) x(~isspace(x)),GeneralParams,'Uni',0);

listing = dir(Directory);
allfiles = {listing.name}';
matchmonkey = cellfun(@(x) length(strfind(x,sprintf('%s',Monkey))),allfiles);
matchdate = cellfun(@(x) length(strfind(x,sprintf('%s',Date(Date ~= '-')))),allfiles);

fileidx = find(matchmonkey & matchdate);%1,'first');
Params = [];
for i = 1:length(fileidx)
    fname = allfiles{fileidx(i)}(1:end-4);
    f = fopen([Directory fname '.xml']);
    c = textscan(f,'%s'); c = horzcat(c{1}{:});
    
    %% Find general params
    for params = 1:length(GeneralParams)
        paramloc = strfind(c,GeneralParams{params});
        catstr = c(paramloc:end);
        valbounds1 = strfind(catstr,'<value>'); valbounds1 = valbounds1(1) + 7;
        valbounds2 = strfind(catstr,'</value>'); valbounds2 = valbounds2(1) - 1;
        value = catstr(valbounds1:valbounds2);
        Params.(sprintf('file_%d',i)).Name = fname;
        Params.(sprintf('file_%d',i)).General.(GeneralParams{params}) = str2double(value);
    end
    
    %% Find Task specific 
    taskparamlocs = strfind(c,tname);
    for params = 1:length(taskparamlocs)
        catstr = c(taskparamlocs(params):end);
        
        paramname_1 = length(tname)+1;
        paramname_2 = strfind(catstr,'</key>'); paramname_2 = paramname_2(1)-1;
        
        valbounds1 = strfind(catstr,'<value>'); valbounds1 = valbounds1(1) + 7;
        valbounds2 = strfind(catstr,'</value>'); valbounds2 = valbounds2(1) - 1;
        
        Params.(sprintf('file_%d',i)).(Task).(catstr(paramname_1:paramname_2)) = str2double(catstr(valbounds1:valbounds2));
    end
    fclose(f);
end