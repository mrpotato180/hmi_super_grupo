function [elim_list, zero_chan] = channel_rejection(path,ica,loc)
%loads the run and returns bad channel numbers
sr=load(fullfile(path));
channels=sr.EEG.data(1:61,:);
fm=512;%Hz
filt_chan=zeros(size(channels));
if ~ica
    [b,a]=butter(4,[0.3/(fm/2) 70/(fm/2)]);
    for i= 1:61
        filt_chan(i,:)=filtfilt(b,a,channels(i,:));
    end
else
    filt_chan=eeglab_ica(path, loc);
    filt_chan=filt_chan(4:64,:);
end


zero_chan=zeros(size(filt_chan));
for i=1:61
    zero_chan(i,:)=detrend(filt_chan(i,:),'constant');
end


kurt_chan=zeros(size(1:61));
for i=1:61
    kurt_chan(i)=kurtosis(zero_chan(i,:));
end

sort_kurt=sort(kurt_chan);
corr_kurt=sort_kurt(7:55);
mean_kurt=mean(corr_kurt);
sdv_kurt=std(corr_kurt);
high_thresh=mean_kurt+5*sdv_kurt;
low_thresh=mean_kurt-5*sdv_kurt;

elim_list=false(size(1:61));
for i=1:61
   if kurt_chan(i)>high_thresh || kurt_chan(i)<low_thresh
    elim_list(i) = true;
    disp(i)
   end
end

end

