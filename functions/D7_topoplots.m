function  D7_topoplots(grand_mean_matrix,locs)

potentials=zeros(61,1);
for ch=1:61
    potentials(ch,1)=min(grand_mean_matrix(ch,2.4*512:2.6*512));
end

topoplot(potentials,locs);
colorbar;
