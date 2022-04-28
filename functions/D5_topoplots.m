function  D5_topoplots(grand_mean_matrix,locs)

[~,latency]=max(grand_mean_matrix(30,:));
time_window=[latency-51, latency+51];
potentials=mean(grand_mean_matrix(:,time_window(1):time_window(2)),2);
topoplot(potentials,locs);
colorbar;