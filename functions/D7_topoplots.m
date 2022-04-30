function  D7_topoplots(grand_mean_matrix,locs)

potentials=zeros(61,1);
for ch=1:61
    potentials(ch,1)=min(grand_mean_matrix(ch,round(2.4*512):round(2.6*512)));
end
topoplot(squeeze(potentials),locs,'maplimits',[-60 0]);
colormap turbo;
colorbar;