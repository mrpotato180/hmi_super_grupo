
%NOTE: This function is used to select the bad channels. In order to
% achieve this, it is using the "channel_rejection" function.

function badchannels_subject=d1_bad_channels(subject_list,eeglab_ICA_bool,locs,datapath,icapath,resfolder)
badchannels_subject=zeros(8,10,61,'logical');
tstart=tic;
if ~eeglab_ICA_bool
    load_path=strcat(resfolder,filesep,"bad_channels.mat");
    if exist(load_path,"file")
    msgbox('The old bad channel file was taken.ICA is off.');
    return
    end
end
if eeglab_ICA_bool
    load_path=strcat(resfolder,filesep,"bad_channels_with_ica.mat");
    if exist(load_path,"file")
    msgbox('The old bad channel file was taken. ICA is on.');
    return
    end
    
end
close(findall(0,'type','figure','tag','TMWWaitbar'))
    for i=1:length(subject_list)
        curr_folder=append(datapath,'\',subject_list(i));
        if eeglab_ICA_bool
            if ~exist(strcat(icapath,filesep,subject_list(i)), 'dir')
                mkdir(strcat(icapath,filesep,subject_list(i)));
            end
        end
        
        for j=1:10
            
            if isempty(findall(0,'type','figure','tag','TMWWaitbar'))
                waiter=waitbar(0,'Calculating...','Name','Progress');
            end
            waitbar(((j-1)+((i-1)*10))/80,waiter,strcat('Preparing run ',string(j), ' from subject  ',string(subject_list(i))),'Name','Progress');
            
            if j<10
                run_code=strcat('0',string(j));
            else
                run_code=string(j);
            end
            curr_subject=subject_list(i);
            if length(char(subject_list(i))) == 2
            subject_code=strcat('0',extractAfter(curr_subject,1));
            else
            subject_code=strcat('1',extractAfter(curr_subject,2));
            end
            filename = strcat('ME_S',string(subject_code),'_r',string(run_code));
            path=append(curr_folder,filesep,filename);
            if eeglab_ICA_bool
                [badchannels_subject(i,j,:), filtered_data]= channel_rejection(path,eeglab_ICA_bool,locs);
            else
                [badchannels_subject(i,j,:),filtered_data]= channel_rejection(path,eeglab_ICA_bool, 'a');
            end
            if eeglab_ICA_bool
                save_path=strcat(icapath,filesep,subject_list(i),filesep,'filtered_data_',string(j),'.mat');
                save(save_path,"filtered_data")
            end
          
        end
   end
        
       
    close(findall(0,'type','figure','tag','TMWWaitbar'))
    
    %C=num2cell(filtered_data_subject,[1,]);
if ~eeglab_ICA_bool
        save_path=strcat(resfolder,filesep,"bad_channels.mat");
        save(save_path,"badchannels_subject")

else
        save_path=strcat(resfolder,filesep,"bad_channels_with_ica.mat");
        save(save_path,"badchannels_subject")
       
        clearvars -global EEG ALLCOM ALLEEG CURRENTSET CURRENTSTUDY PLUGINLIST STUDY LASTCOM

end
tEnd = toc(tstart);
disp(['Elapsed time is ' num2str(floor(tEnd/60))  ' minutes and ' num2str(rem(tEnd,60)) ' seconds']);
