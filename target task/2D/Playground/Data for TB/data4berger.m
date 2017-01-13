task = 'CF';
monkey = 'Chewie';
date = '10-05-2016';
switch task
    case 'UNT2D' % Uncertainty
        DAT = [];   

        DAT.meta.monkey = monkey;
        DAT.meta.date = date;
        DAT.meta.task = 'Uncertainty task with Center-Out (multiple target distributions)';
        paramlist = {'UT2DOuterTargetSize' ,  'targetsize'         ,'degrees';...
                     'UT2DMovementLength'  ,  'ringradius'         ,'cm'     ;...
                     'UT2DTargetDiameter'  ,  'centertargetradius' ,'cm'     ;...
                     'UT2DOuterTargetDepth',  'ringdepth'          ,'cm'};

        % find meta info
        dirloc = 'Z:\Animal-Miscellany\Mihili 12A3\Behavior Parameters\';
        
        listing = dir(dirloc);
        allfiles = {listing.name}';
        matchmonkey = cellfun(@(x) length(strfind(x,sprintf('%s',DAT.meta.monkey))),allfiles);
        matchdate = cellfun(@(x) length(strfind(x,sprintf('%s',DAT.meta.date(DAT.meta.date ~= '-')))),allfiles);

        fileidx = find(matchmonkey & matchdate,1,'first');

        f = fopen([dirloc allfiles{fileidx}]);
        c = textscan(f,'%s'); c = c{1};
        ccomb = horzcat(c{:});
        for params = 1:size(paramlist,1)
            paramloc = strfind(ccomb,paramlist{params,1});
            catstr = ccomb(paramloc:end);
            valbounds1 = strfind(catstr,'<value>'); valbounds1 = valbounds1(1) + 7;
            valbounds2 = strfind(catstr,'</value>'); valbounds2 = valbounds2(1) - 1;
            value = catstr(valbounds1:valbounds2);

            DAT.meta.(paramlist{params,2}) = [value ' ' paramlist{params,3}];
        end

        % Kinematics
        DAT.kinematics = table(alldays(1).kin.pos(:,1),alldays(1).kin.pos(:,2),alldays(1).kin.pos(:,3),...
                                'VariableNames',{'Time','XPosition','YPosition'});

        % Trials
        triallengths = [0 cellfun(@(x) size(x,1),{alldays.tt})]; 
        for i = 1:length(alldays)          

            TargetAngle = alldays(i).tt(:,2);
            CenterTargetOnTime = alldays(i).tt(:,4);
            OuterCueOnTime = alldays(i).tt(:,5);
            GoCueTime = alldays(i).tt(:,6);
            EndOfTrialTime = alldays(i).tt(:,7);
            Resultmat = alldays(i).tt(:,8)==32;
            Result = cell(size(Resultmat,1),1);
            Result(Resultmat==1) = {'S'}; Result(Resultmat==0) = {'F'};
            Result = cell2mat(Result);
            VisualCueAngles = alldays(i).slices;
        %     if(size(VisualCueAngles,2)==1); VisualCueAngles = num2cell(VisualCueAngles); end
            Trial = (1:length(TargetAngle))' + sum(triallengths(1:i));

            DAT.TrialBlocks{i,:} = ...
                table(Trial,TargetAngle,VisualCueAngles,CenterTargetOnTime,OuterCueOnTime,GoCueTime,...
                      EndOfTrialTime,Result);     
        end

        % Set up PMd
        for i = 1:length(alldays(1).PMd_units)
            if isfield(alldays(1),'PMd_locations')
                DAT.PMd{i,1}.bank = alldays(1).PMd_locations{i}(1);
                DAT.PMd{i,1}.pin = alldays(1).PMd_locations{i}(2:end);
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i};
            else
                DAT.PMd{i,1}.channel = floor(alldays(1).PMd_units{i}(1));
                DAT.PMd{i,1}.unit = round(10*(alldays(1).PMd_units{i}(1) - floor(alldays(1).PMd_units{i}(1))));
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i}(2:end);
            end
        end

        % Set up M1
        for i = 1:(length(alldays(1).M1_units))
            if isfield(alldays(1),'M1_locations')
                DAT.M1{i,1}.bank = alldays(1).M1_locations{i}(1);
                DAT.M1{i,1}.pin = alldays(1).M1_locations{i}(2:end);
                DAT.M1{i,1}.timestamps = alldays(1).M1_units{i};
            else
                if alldays(1).M1_units{i}(1) < 100
                    DAT.M1{i,1}.channel = floor(alldays(1).M1_units{i}(1));
                    DAT.M1{i,1}.unit = round(10*(alldays(1).M1_units{i}(1) - floor(alldays(1).M1_units{i}(1))));
                    DAT.M1{i,1}.timestamps = alldays(1).M1_units{i}(2:end);
                end
            end
        end
        clearvars -except alldays DAT
        %
        %d = 2;
        % 
        % blck = find(cellfun(@(x) sum(x.Trial==trl2plot),DAT.TrialBlocks));
        % trlblck = find(DAT.TrialBlocks{blck}.Trial == trl2plot);
        % thets = linspace(0,2*pi,100); figure; hold on; 
        % r1 = str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1));
        % r2 = str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1));
        % patch([r1*cos(thets) r2*cos(fliplr(thets))],[r1*sin(thets) r2*sin(fliplr(thets))],[.25 .25 1],'EdgeColor',[.25 .25 1]); 
        % axis square; axis equal; axis off;
        % if size(DAT.TrialBlocks{blck}.VisualCueAngles,2)>1
        %     plot([r1*cos(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:));r2*cos(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:))],...
        %          [r1*sin(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:));r2*sin(DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:))],'r','LineWidth',5);
        % else
        %     a1 = DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,' ')-1));
        %     a2 = DAT.TrialBlocks{blck}.VisualCueAngles(trlblck,:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,' ')-1));
        %     patch([r1*cos(a1:pi/50:a2) r2*cos(fliplr(a1:pi/50:a2))],[r1*sin(a1:pi/50:a2) r2*sin(fliplr(a1:pi/50:a2))],'r','EdgeColor','r'); 
        % end


        %d = 2;
        %
        % plotfunc = @(DAT,trial) eval(sprintf('figure; hold on; title([''Trial: '' num2str(%d)],''FontSize'',18); patch([(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*cos(linspace(0,2*pi,100)) (str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*cos(fliplr(linspace(0,2*pi,100)))],[(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*sin(linspace(0,2*pi,100)) (str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*sin(fliplr(linspace(0,2*pi,100)))],[.25 .25 1],''EdgeColor'',[.25 .25 1]); axis square; axis equal; axis off; if size(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles,2)>1; plot([(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*cos(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:));(str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*cos(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:))],[(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*sin(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:));(str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*sin(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:))],''r'',''LineWidth'',5); else patch([(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*cos((DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1))):pi/50:(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1)))) (str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*cos(fliplr((DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1))):pi/50:(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1)))))],[(str2double(DAT.meta.ringradius(1))-.5*str2double(DAT.meta.ringdepth(1)))*sin((DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1))):pi/50:(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1)))) (str2double(DAT.meta.ringradius(1))+.5*str2double(DAT.meta.ringdepth(1)))*sin(fliplr((DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)-.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1))):pi/50:(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.VisualCueAngles(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d),:)+.5*pi/180*str2double(DAT.meta.targetsize(1:strfind(DAT.meta.targetsize,'' '')-1)))))],''r'',''EdgeColor'',''r''); end; plot(downsample(DAT.kinematics.XPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d))),10),downsample(DAT.kinematics.YPosition(DAT.kinematics.Time>DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.OuterCueOnTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d)) & DAT.kinematics.Time<DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.EndOfTrialTime(find(DAT.TrialBlocks{find(cellfun(@(x) sum(x.Trial==%d),DAT.TrialBlocks))}.Trial == %d))),10),''k'',''LineWidth'',5);',repmat(trial,1,50)));
        
    case 'UCK' % Cisek
        DAT = [];

        DAT.meta.monkey = monkey;
        DAT.meta.date = date;
        DAT.meta.task = 'Two target task with delay';
        paramlist = {'UCKMoveLength' ,  'reachlength'         ,'cm'     ;...
                     'UCKCenterSize' ,  'centertargetradius'  ,'cm'     ;...
                     'UCKOuterSize'  ,  'outertargetradius'   ,'cm'     ;...
                     'UCKCueOneR'    ,  'cue1colorR'          ,''       ;...
                     'UCKCueOneG'    ,  'cue1colorG'          ,''       ;...
                     'UCKCueOneB'    ,  'cue1colorB'          ,''       ;...
                     'UCKCueTwoR'    ,  'cue2colorR'          ,''       ;...
                     'UCKCueTwoG'    ,  'cue2colorG'          ,''       ;...
                     'UCKCueTwoB'    ,  'cue2colorB'          ,''       ;...
                     'UCKCueOffset'  ,  'correctcueradius'    ,'cm'     ;};

        % find meta info
        if strcmpi(DAT.meta.monkey,'Mihili')
            dirloc = 'Z:\Animal-Miscellany\Mihili 12A3\Behavior Parameters\';
        else
            dirloc = 'Z:\Animal-Miscellany\Chewie 8I2\Behavior Parameters\UCK\';
        end
        listing = dir(dirloc);
        allfiles = {listing.name}';
        matchmonkey = cellfun(@(x) length(strfind(x,sprintf('%s',DAT.meta.monkey))),allfiles);
        matchdate = cellfun(@(x) length(strfind(x,sprintf('%s',DAT.meta.date(DAT.meta.date ~= '-')))),allfiles);

        fileidx = find(matchmonkey & matchdate,1,'first');

        f = fopen([dirloc allfiles{fileidx}]);
        c = textscan(f,'%s'); c = c{1};
        ccomb = horzcat(c{:});
        for params = 1:size(paramlist,1)
            paramloc = strfind(ccomb,paramlist{params,1});
            catstr = ccomb(paramloc:end);
            valbounds1 = strfind(catstr,'<value>'); valbounds1 = valbounds1(1) + 7;
            valbounds2 = strfind(catstr,'</value>'); valbounds2 = valbounds2(1) - 1;
            value = catstr(valbounds1:valbounds2);

            DAT.meta.(paramlist{params,2}) = [value ' ' paramlist{params,3}];
        end
        DAT.meta.Cue1Color = ['RGB = [',DAT.meta.cue1colorR, DAT.meta.cue1colorG,DAT.meta.cue1colorB,']'];
        DAT.meta.Cue2Color = ['RGB = [',DAT.meta.cue2colorR, DAT.meta.cue2colorG,DAT.meta.cue2colorB,']'];
        DAT.meta = rmfield(DAT.meta,{'cue1colorR','cue1colorG','cue1colorB','cue2colorR',...
                                     'cue2colorG','cue2colorB'});
        % Kinematics
        DAT.kinematics = table(alldays(1).kin.pos(:,1),alldays(1).kin.pos(:,2),alldays(1).kin.pos(:,3),...
                                'VariableNames',{'Time','XPosition','YPosition'});

        % Trials
        tcorrect = @(t1,t2,tc) (abs(circ_dist(t1,tc))<eps)+2*(abs(circ_dist(t2,tc))<.01);
        for i = 1:length(alldays)
            if isfield(alldays(i),'numtarg') 
                if ~isempty(strfind(alldays(i).numtarg,'catch'))
                    alldays(i).tt(:,8) = -1;
                end
            end
        end

        TT = vertcat(alldays.tt);  
        if strcmp(alldays(1).labels{1},'1: number()') % We're in new format and need to rearrange columns
            TTN = zeros(size(TT,1),18);
            
            TTN(:,2:3) = TT(:,11:12)*pi/180;
            TTN(:,4:9) = TT(:,4:9);
            TTN(:,10) = nanmean(TT(:,[10 3]),2);
            TTN(:,13) = TT(:,15);
            TTN(:,14) = (TT(:,18)==2)*.5 + (TT(:,18)==1)*1000;
            rd = TT(:,19)*pi/180;
            TTN(:,16) = 32*(abs(circ_dist(rd,TTN(:,13)))<(pi/4));
            
            TT = TTN;
        end
        [~,torder] = sortrows(TT(:,9));
        TT = TT(torder,:);
        
        Trial = (1:size(TT,1))';
        Cue1Angle = TT(:,2); 
        Cue2Angle = TT(:,3);
        CuedTarget = tcorrect(TT(:,2),TT(:,3),TT(:,13));
        numtargets = (TT(:,14)>1) + 2*(TT(:,14)<1);
        Cue1Angle(CuedTarget==2 & numtargets==1) = NaN;
        Cue2Angle(CuedTarget==1 & numtargets==1) = NaN;
        CenterTargetOnTime = TT(:,4);
        CenterTargetHoldTime = TT(:,5);
        OuterCueOnTime = TT(:,6);
        OuterCueOffTime = TT(:,7);
        CorrectTargetOnTime = TT(:,8);
        GoCueTime = TT(:,9);
        EndOfTrialTime = TT(:,10);
        Resultmat = TT(:,16)==32;
        Result = cell(size(Resultmat,1),1);
        Result(Resultmat==1) = {'S'}; Result(Resultmat==0) = {'F'};
        Result = cell2mat(Result);


        DAT.Trials = ...
            table(Trial,Cue1Angle,Cue2Angle,CuedTarget,CenterTargetOnTime,CenterTargetHoldTime,...
                  OuterCueOnTime,OuterCueOffTime,CorrectTargetOnTime,GoCueTime,EndOfTrialTime,Result);     


        % Set up PMd
        for i = 1:length(alldays(1).PMd_units)
            if isfield(alldays(1),'PMd_locations')
                DAT.PMd{i,1}.bank = alldays(1).PMd_locations{i}(1);
                DAT.PMd{i,1}.pin = alldays(1).PMd_locations{i}(2:end);
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i};
            else
                DAT.PMd{i,1}.channel = floor(alldays(1).PMd_units{i}(1));
                DAT.PMd{i,1}.unit = round(10*(alldays(1).PMd_units{i}(1) - floor(alldays(1).PMd_units{i}(1))));
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i}(2:end);
            end
        end

        % Set up M1
        for i = 1:(length(alldays(1).M1_units))
            if isfield(alldays(1),'M1_locations')
                DAT.M1{i,1}.bank = alldays(1).M1_locations{i}(1);
                DAT.M1{i,1}.pin = alldays(1).M1_locations{i}(2:end);
                DAT.M1{i,1}.timestamps = alldays(1).M1_units{i};
            else
                if alldays(1).M1_units{i}(1) < 100
                    DAT.M1{i,1}.channel = floor(alldays(1).M1_units{i}(1));
                    DAT.M1{i,1}.unit = round(10*(alldays(1).M1_units{i}(1) - floor(alldays(1).M1_units{i}(1))));
                    DAT.M1{i,1}.timestamps = alldays(1).M1_units{i}(2:end);
                end
            end
        end
        clearvars -except alldays DAT
        
    case 'CF' % Curl Field
        DAT = [];   

        DAT.meta.monkey = monkey;
        DAT.meta.date = date;
        DAT.meta.task = 'Center Out reaches with velocity-dependent curl field';

%         DAT.meta.;
        DAT.meta.reachlength = num2str(alldays(1).reachlength);
        DAT.meta.targetshape = 'square';
        DAT.meta.targetwidth = num2str(alldays(1).targetsize);
        fangle = input('please input force angle parameter value\n>>');
        fmag = input('please input force magnitude value\n>>');
        DAT.meta.forceangle = [num2str((pi-fangle)/pi*180) ' degrees (from inst. movement direction)'];
        DAT.meta.forcemagnitude = [num2str(fmag) ' N*s/cm'];

        % Kinematics
        DAT.kinematics = table(alldays(1).kin.pos(:,1),alldays(1).kin.pos(:,2),alldays(1).kin.pos(:,3),...
                                'VariableNames',{'Time','XPosition','YPosition'});

        % Trials
        TrialTypes = {'Baseline','Force Field','Washout'};
        triallengths = [0 cellfun(@(x) size(x,1),{alldays.tt})]; 
        for i = 1:length(alldays)          

            TargetAngle = mod(alldays(i).tt(:,8),2*pi);
            OuterCueOnTime = alldays(i).tt(:,5);
            GoCueTime = alldays(i).tt(:,6);
            EndOfTrialTime = alldays(i).tt(:,3);
            Resultmat = alldays(i).tt(:,4)==32;
            Result = cell(size(Resultmat,1),1);
            Result(Resultmat==1) = {'S'}; Result(Resultmat==0) = {'F'};
            Result = cell2mat(Result);
        %     if(size(VisualCueAngles,2)==1); VisualCueAngles = num2cell(VisualCueAngles); end
            Trial = (1:length(TargetAngle))' + sum(triallengths(1:i));
            
            DAT.TrialBlocks{i,:}.description = TrialTypes{i};
            DAT.TrialBlocks{i,:}.trials = ...
                table(Trial,TargetAngle,OuterCueOnTime,GoCueTime,...
                      EndOfTrialTime,Result);  
        end

        % Set up PMd
        for i = 1:length(alldays(1).PMd_units)
            if isfield(alldays(1),'PMd_locations')
                DAT.PMd{i,1}.bank = alldays(1).PMd_locations{i}(1);
                DAT.PMd{i,1}.pin = alldays(1).PMd_locations{i}(2:end);
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i};
            else
                DAT.PMd{i,1}.channel = floor(alldays(1).PMd_units{i}(1));
                DAT.PMd{i,1}.unit = round(10*(alldays(1).PMd_units{i}(1) - floor(alldays(1).PMd_units{i}(1))));
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i}(2:end);
            end
        end

        % Set up M1
        for i = 1:(length(alldays(1).M1_units))
            if isfield(alldays(1),'M1_locations')
                DAT.M1{i,1}.bank = alldays(1).M1_locations{i}(1);
                DAT.M1{i,1}.pin = alldays(1).M1_locations{i}(2:end);
                DAT.M1{i,1}.timestamps = alldays(1).M1_units{i};
            else
                if alldays(1).M1_units{i}(1) < 100
                    DAT.M1{i,1}.channel = floor(alldays(1).M1_units{i}(1));
                    DAT.M1{i,1}.unit = round(10*(alldays(1).M1_units{i}(1) - floor(alldays(1).M1_units{i}(1))));
                    DAT.M1{i,1}.timestamps = alldays(1).M1_units{i}(2:end);
                end
            end
        end
        clearvars -except alldays DAT 
        
    case 'VR' % Visual Rotation
        DAT = [];   

        DAT.meta.monkey = monkey;
        DAT.meta.date = date;
        DAT.meta.task = 'Center Out reaches with visuomotor rotation';

        %         DAT.meta.;
        DAT.meta.reachlength = num2str(alldays(1).reachlength);
        DAT.meta.targetshape = 'square';
        DAT.meta.targetwidth = num2str(alldays(1).targetsize);
        rangle = input('please input visual rotation parameter value\n>>');
        DAT.meta.rotation = [num2str(rangle/pi*180) ' degrees'];
        DAT.meta.note = 'Kinematic data reflects the location of the hand, not the cursor';

        % Kinematics
        DAT.kinematics = table(alldays(1).kin.pos(:,1),alldays(1).kin.pos(:,2),alldays(1).kin.pos(:,3),...
                                'VariableNames',{'Time','XPosition','YPosition'});

        % Trials
        TrialTypes = {'Baseline','Visual Rotation','Washout'};
        triallengths = [0 cellfun(@(x) size(x,1),{alldays.tt})]; 
        for i = 1:length(alldays)          

            TargetAngle = mod(alldays(i).tt(:,8),2*pi);
            OuterCueOnTime = alldays(i).tt(:,5);
            GoCueTime = alldays(i).tt(:,6);
            EndOfTrialTime = alldays(i).tt(:,3);
            Resultmat = alldays(i).tt(:,4)==32;
            Result = cell(size(Resultmat,1),1);
            Result(Resultmat==1) = {'S'}; Result(Resultmat==0) = {'F'};
            Result = cell2mat(Result);
        %     if(size(VisualCueAngles,2)==1); VisualCueAngles = num2cell(VisualCueAngles); end
            Trial = (1:length(TargetAngle))' + sum(triallengths(1:i));

            DAT.TrialBlocks{i,:}.description = TrialTypes{i};
            DAT.TrialBlocks{i,:}.trials = ...
                table(Trial,TargetAngle,OuterCueOnTime,GoCueTime,...
                      EndOfTrialTime,Result);  
        end

        % Set up PMd
        for i = 1:length(alldays(1).PMd_units)
            if isfield(alldays(1),'PMd_locations')
                DAT.PMd{i,1}.bank = alldays(1).PMd_locations{i}(1);
                DAT.PMd{i,1}.pin = alldays(1).PMd_locations{i}(2:end);
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i};
            else
                DAT.PMd{i,1}.channel = floor(alldays(1).PMd_units{i}(1));
                DAT.PMd{i,1}.unit = round(10*(alldays(1).PMd_units{i}(1) - floor(alldays(1).PMd_units{i}(1))));
                DAT.PMd{i,1}.timestamps = alldays(1).PMd_units{i}(2:end);
            end
        end

        % Set up M1
        for i = 1:(length(alldays(1).M1_units))
            if isfield(alldays(1),'M1_locations')
                DAT.M1{i,1}.bank = alldays(1).M1_locations{i}(1);
                DAT.M1{i,1}.pin = alldays(1).M1_locations{i}(2:end);
                DAT.M1{i,1}.timestamps = alldays(1).M1_units{i};
            else
                if alldays(1).M1_units{i}(1) < 100
                    DAT.M1{i,1}.channel = floor(alldays(1).M1_units{i}(1));
                    DAT.M1{i,1}.unit = round(10*(alldays(1).M1_units{i}(1) - floor(alldays(1).M1_units{i}(1))));
                    DAT.M1{i,1}.timestamps = alldays(1).M1_units{i}(2:end);
                end
            end
        end
        clearvars -except alldays DAT 
end
