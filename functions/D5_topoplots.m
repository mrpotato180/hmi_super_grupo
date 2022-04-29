function  D5_topoplots(grand_mean_matrix,locs)

[~,latency]=max(grand_mean_matrix(30,2*512:2.35*512));


time_window=[2*512+latency-51, 2*512+latency+51];
disp(time_window/512)
potentials=mean(grand_mean_matrix(:,time_window(1):time_window(2)),2);
topoplot(potentials,locs);
colorbar;