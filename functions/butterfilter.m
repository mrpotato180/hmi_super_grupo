function filtersignal = butterfilter(EEGsignals,highthresh)
fm=512;%Hz
filtersignal=zeros(size(EEGsignals));
[b,a]=butter(4,[0.3/(fm/2) highthresh/(fm/2)]);
for i= 1:61
    filtersignal(i,:)=filtfilt(b,a,EEGsignals(i,:));
end