
%NOTE: This function is used to select the bad channels. In order to
% achieve this, it is using the "channel_rejection" function.

function badchannels_subject=d1_bad_channels(subject_list,eeglab_ICA_bool,locs,datapath,icapath)
badchannels_subject=zeros(8,10,61,'int8');
tstart=tic;
if ~eeglab_ICA_bool
    if exist("bad_channels.mat","file")
    msgbox('The old bad channel file was taken.ICA is off.');
    return
    end
end
if eeglab_ICA_bool
    if exist("icapath","file") % esto necesita atención
    msgbox('The old bad channel file was taken. ICA is on.');
    return
    end
    
end
close(findall(0,'type','figure','tag','TMWWaitbar'))
    for i=1:length(subject_list)
        curr_folder=append(datapath,'\',subject_list(i));
        filenames=dir(curr_folder);
        if eeglab_ICA_bool
            if ~exist(strcat(icapath,'\',subject_list(i)), 'dir')
                mkdir(strcat(icapath,'\',subject_list(i)));
            end
        end
        disp(subject_list(i))
        
        for j=1:10
            
            if isempty(findall(0,'type','figure','tag','TMWWaitbar'))
                waiter=waitbar(0,'Calculating...','Name','Progress');
            end
            waitbar(((j-1)+((i-1)*10))/80,waiter,strcat('Preparing run ',string(j), ' from subject  ',string(subject_list(i))),'Name','Progress');
            filename = filenames(j+2).name;
            path=append(curr_folder,'\',filename);
            if eeglab_ICA_bool
                [badchannels_subject(i,j,:), filtered_data]= channel_rejection(path,eeglab_ICA_bool,locs);
            else
                [badchannels_subject(i,j,:),filtered_data]= channel_rejection(path,eeglab_ICA_bool, 'a');
            end
            if eeglab_ICA_bool
                save_path=strcat(icapath,'\',subject_list(i),'\filtered_data_',string(j),'.mat');
                save(save_path,"filtered_data")
            end
          
        end
   end
        
       
    close(findall(0,'type','figure','tag','TMWWaitbar'))
    
    %C=num2cell(filtered_data_subject,[1,]);
if ~eeglab_ICA_bool

        save("bad_channels.mat","badchannels_subject")

    else
        save("bad_channels_with_ica.mat","badchannels_subject")
        clearvars -global globalvars{:} globalvars

end
tEnd = toc(tstart);
disp(['Elapsed time is ' num2str(floor(tEnd/60))  ' minutes and ' num2str(rem(tEnd,60)) ' seconds']);
