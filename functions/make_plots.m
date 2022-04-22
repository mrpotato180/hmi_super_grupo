function make_plots(ms,mschannel,gm,gmchannel,subject_list)
fm=512;
t=-2.5+1/fm:1/fm:6;
figure;
tiledlayout(4,2,"Padding","compact","TileSpacing","tight");
for i=1:length(ms(1,1,:))
    nexttile
    
    
    plot(t,ms(mschannel,:,i))
    xlim([-2.5 6])
    xline(-2)
    xline(0)
    title(subject_list(i))
end

figure;
plot(t,gm(gmchannel,:))
xlim([-2.5 6])
xline(-2)
xline(0)
title('Grand mean')
xlabel('Time (s)')
ylabel('Intensity (*10^-6 V)')
end
