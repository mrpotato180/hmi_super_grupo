
function eegout=eeglab_ica(path,loc)
%a function that automatically rejects components calculated by an acsobiro
%ica decomposition
load(fullfile(path),"-mat",'EEG')
a=cat(1,EEG.data(62:64,:),EEG.data(1:61,:)); %add the EOG channels to the beginning of the variable
a=butterfilter(a,70);%filter the signal, it was found it is better to filter before ica, not after

%this was done because when performing the ica rejection by hand, it was
%noticed that this produces better results
save("temp.mat","a","-mat") 
%we had to save it to a temporary file 
% because variable loading in eeglab works only with global variables
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; close; %we open eeglab and close the gui

EEG = pop_importdata('dataformat','matlab','nbchan',0,'data','temp.mat','srate',512,'pnts',0,'xmin',0,'chanlocs',loc);
%we have now imported the data and channels into eeglab with the correct
%sampling rate
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
EEG = eeg_checkset( EEG );

%below is the typical eeglab ica procedure for any algorithm
EEG.icawinv = acsobiro(EEG.data);
if isempty(EEG.icaweights)
    EEG.icaweights = pinv(EEG.icawinv);
end
if isempty(EEG.icasphere)
    EEG.icasphere  = eye(size(EEG.icaweights,2));
end
if isempty(EEG.icawinv)
    EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere); 
end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%finally, we label the components by type
EEG = pop_iclabel(EEG, 'default');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
comp_rej_lst=[];

for k=1:EEG.nbchan
    if EEG.etc.ic_classification.ICLabel.classifications(k,3)>0.85
    %if EEG.etc.ic_classification.ICLabel.classifications(k,3)>0.5
        comp_rej_lst=[comp_rej_lst,k];
    end
    for i=[2,4,5,6,7]
        if EEG.etc.ic_classification.ICLabel.classifications(k,i)>0.85
            comp_rej_lst=[comp_rej_lst,k];
        end
    end
end
EEG = pop_subcomp( EEG, comp_rej_lst, 0);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

eegout=EEG.data;

