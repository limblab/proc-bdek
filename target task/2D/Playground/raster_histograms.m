dircol = 19;
align_col = 6;
window_bounds = [0 750];

[~,~,~,fulldatPDS] = tuning_types_PD_func(alldays,'PMd',1,dircol,{align_col},{[0 750]});
PDS = get_PDS_fulldat(fulldatPDS,pi/2);

PDnum = mod(round((PDS+2*pi)./(pi/4))-1,8)+1;
Tnum{1} = mod(round((alldays(1).tt(:,dircol)+2*pi)./(pi/4))-1,8)+1;
Tnum{2} = mod(round((alldays(2).tt(:,dircol)+2*pi)./(pi/4))-1,8)+1;
%%
alldays(1).tt(isnan(alldays(1).tt(:,20)),20) = alldays(1).tt(isnan(alldays(1).tt(:,20)),9);
alldays(2).tt(isnan(alldays(2).tt(:,20)),20) = alldays(2).tt(isnan(alldays(2).tt(:,20)),9);
rast = cell(1,2);
getPD = @(pd,targs) pd == (mod(targs-1,8)+1);
getOD = @(pd,targs) pd == (mod(targs+4-1,8)+1);
getORTH = @(pd,targs) pd == (mod(targs+2-1,8)+1) | pd == (mod(targs+6-1,8)+1);
for i = 1:length(PDS)
    clc; fprintf('getting raster %d/%d\n',i,length(PDS));
    if ~isnan(PDS(i))
        
        for j = 1:2
            pdind = find(getPD(PDnum(i),Tnum{j}));
            odind = find(getOD(PDnum(i),Tnum{j}));
            orthind = find(getORTH(PDnum(i),Tnum{j}));

            r1 = raster_get(alldays(1).PMd_units{i},alldays(j).tt,[-100 750]/1000,6);
            rast{j}.pd{i,1} = bin_array(r1(pdind,:),length(pdind),size(r1,2)./10,'sum');
            rast{j}.od{i,1} = bin_array(r1(odind,:),length(odind),size(r1,2)./10,'sum');
            rast{j}.orth{i,1} = bin_array(r1(orthind,:),length(orthind),size(r1,2)./10,'sum');

            r1 = raster_get(alldays(1).PMd_units{i},alldays(j).tt,[0 250]/1000,7);
            rast{j}.pd{i,2} = bin_array(r1(pdind,:),length(pdind),size(r1,2)./10,'sum');
            rast{j}.od{i,2} = bin_array(r1(odind,:),length(odind),size(r1,2)./10,'sum');
            rast{j}.orth{i,2} = bin_array(r1(orthind,:),length(orthind),size(r1,2)./10,'sum');

            r1 = raster_get(alldays(1).PMd_units{i},alldays(j).tt,[0 250]/1000,8);
            rast{j}.pd{i,3} = bin_array(r1(pdind,:),length(pdind),size(r1,2)./10,'sum');
            rast{j}.od{i,3} = bin_array(r1(odind,:),length(odind),size(r1,2)./10,'sum');
            rast{j}.orth{i,3} = bin_array(r1(orthind,:),length(orthind),size(r1,2)./10,'sum');

            r1 = raster_get(alldays(1).PMd_units{i},alldays(j).tt,[-100 200]/1000,20);
            rast{j}.pd{i,4} = bin_array(r1(pdind,:),length(pdind),size(r1,2)./10,'sum');
            rast{j}.od{i,4} = bin_array(r1(odind,:),length(odind),size(r1,2)./10,'sum');
            rast{j}.orth{i,4} = bin_array(r1(orthind,:),length(orthind),size(r1,2)./10,'sum');

            clear r1 r2
        end
    end
end
