%% smooth and center position data
gocol = 6;%9;
endcol = 3;%3;

% still_1 = 5;%5;
% still_2 = 5;%6;


TS = alldays(1).kin.pos(:,1);
dt = diff(alldays(1).kin.pos(1:2,1));
% POS = [smooth(alldays(1).kin.pos(:,2),50) smooth(alldays(1).kin.pos(:,3),50)];
POS = alldays(1).kin.pos(:,2:3);
% poscent = zeros(size(alldays(1).tt,1),2);
% 
% combtt = vertcat(alldays.tt);
% for tri = 1:size(combtt)
%     i0 = find(TS > combtt(tri,still_1),1,'first');
%     i1 = find(TS < combtt(tri,still_2),1,'last');
%     poscent(tri,:) = nanmean(POS(i0:i1,:),1);
% end
% POS = POS - repmat(nanmean(poscent,1),size(POS,1),1);

dfunc = @(x) [nan; nan; (3/2)*x(3:end) - 2*x(2:(end-1)) + (1/2)*x(1:(end-2))];
SPD = sqrt(dfunc(POS(:,1)).^2 + dfunc(POS(:,2)).^2)./dt;
RT = cell(length(alldays),1);
for tb = 1:length(alldays)
    TT = alldays(tb).tt;
    [react_time,topspeed, time_topspeed] = deal(nan(size(TT,1),1));

    for tri = 1:size(TT,1)

        igo = find(TS < TT(tri,gocol),1,'last');
        iend = find(TS < TT(tri,endcol),1,'last');
        tis = 1:(iend-igo+1);
        SPDsnip = SPD(igo:iend);
        [pkspd,pkloc] = max(SPDsnip);
        befpk = SPDsnip(pkloc:-1:1);
        befts = pkloc:-1:1;
        lastunder25 = befts(befpk < (0.25*pkspd));

        if ~isempty(lastunder25)
            react_time(tri) = lastunder25(1);
        else
            react_time(tri) = NaN;
        end

        clc; fprintf('%d/%d\n',tri,length(TT));
    end
    RT{tb} = react_time;
    alldays(tb).RT = alldays(tb).tt(:,gocol)+RT{tb}.*dt;
end
clear POS SPD TS