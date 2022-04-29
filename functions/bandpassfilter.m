function filtersignal = bandpassfilter(EEGsignals,lowthresh,highthresh)
fm=512;%Hz
filtersignal=zeros(size(EEGsignals));
[b,a]=butter(4,[lowthresh/(fm/2) highthresh/(fm/2)]);
for i= 1:length(EEGsignals(:,1))
    filtersignal(i,:)=filtfilt(b,a,EEGsignals(i,:));
end