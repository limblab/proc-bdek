dirbin = 10;
timebin = 15;

Lrast(isnan(Lrast))=0;
Hrast(isnan(Hrast))=0;

Ltt = comp_tt(comp_tt(:,3)==max(comp_tt(:,3)),:);
Htt = comp_tt(comp_tt(:,3)==min(comp_tt(:,3)),:);

[dirsL, ordered_inds_L] = sortrows(Ltt(:,end));
[dirsH, ordered_inds_H] = sortrows(Htt(:,end));

ord_Ltt = Ltt(ordered_inds_L,:);
ord_Htt = Htt(ordered_inds_H,:);

ord_Lrast = Lrast(ordered_inds_L,:);
ord_Hrast = Hrast(ordered_inds_H,:);

% Create direction bins
lowedge = min([dirsL; dirsH]);
highedge = max([dirsL; dirsH]);

dir_binedge = linspace(lowedge,highedge,dirbin+1);
dir_bincent = (dir_binedge(1:end-1) + dir_binedge(2:end))./2;

% Time bin the arrays
tb_Lrast = bin_array(ord_Lrast,size(ord_Lrast,1),timebin);
tb_Hrast = bin_array(ord_Hrast,size(ord_Hrast,1),timebin);

% Loop through and bin directions
fin_Lrast = zeros(dirbin,timebin);
fin_Hrast = zeros(dirbin,timebin);
for i = 1:dirbin
    fin_Lrast(i,:) = mean(tb_Lrast((ord_Ltt(:,end) > dir_binedge(i) & ord_Ltt(:,end) <= dir_binedge(i+1)),:),1);
    fin_Hrast(i,:) = mean(tb_Hrast((ord_Htt(:,end) > dir_binedge(i) & ord_Htt(:,end) <= dir_binedge(i+1)),:),1);
end

fin_Lrast(isnan(fin_Lrast))=0;
fin_Hrast(isnan(fin_Hrast))=0;


