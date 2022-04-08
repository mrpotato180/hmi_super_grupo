function d2_bad_trials(subject_list,movement_code,eeglab_ica_bool,auto_detect_movement,datapath,icapath)
cont=0;
cont2=1;
count3=1;
matrix=zeros(61,4352,60,8);
events_matrix_1536=zeros(6,3,10,8);
fstart=tic;
%dimensiones : channels / samples / trial / subject
if eeglab_ica_bool
    if exist("bad_channels_with_ica.mat","file")
        bad_channels=load("bad_channels_with_ica.mat","badchannels_subject");
        bad_channels=bad_channels.badchannels_subject;
    else
        waitfor(msgbox('D1 has not been executed. Please execute D1 first.'))
        return
    end
else
    if exist("bad_channels.mat","file")
        bad_channels=load("bad_channels.mat","badchannels_subject");
        bad_channels=bad_channels.badchannels_subject;
    else
        waitfor(msgbox('D1 has not been executed. Please execute D1 first.'))
        return
    end
end
    
if movement_code==1536
    if eeglab_ica_bool
        if exist("bad_trials_1536_with_ica.mat","file")
            return
        else
            boxmsg="Calculating bad trials for Elbow flexion";
        end
    else
        if exist("bad_trials_1536.mat","file")
            return
        else
            boxmsg="Calculating bad trials for Elbow flexion";
        end
    end
end
if movement_code==1541
    if eeglab_ica_bool
        if exist("bad_trials_1541_with_ica.mat","file")
            return
        else
            boxmsg="Calculating bad trials for Hand Opening";
        end
    else
        if exist("bad_trials_1541.mat","file")
            return
        else
            boxmsg="Calculating bad trials for Hand Opening";
        end
    end
end

close(findall(0,'type','figure','tag','TMWWaitbar'))

for i=1:length(subject_list)
    if eeglab_ica_bool
        icapath2=append(icapath,'\',subject_list(i));
        icanames=dir(icapath2);
        
    end
    curr_folder=append(datapath,'\',subject_list(i));
    
    filenames=dir(curr_folder);
    for j=1:10
        if isempty(findall(0,'type','figure','tag','TMWWaitbar'))
                waiter=waitbar(0,'Calculating...','Name',boxmsg);
        end
        waitbar(((j-1)+((i-1)*10))/80,waiter,strcat('Preparing run ',string(j), ' from subject  ',string(subject_list(i))),'Name',boxmsg);
        filename = filenames(j+2).name;
        path=append(curr_folder,'\',filename);
        sr=load(path);
        Events=sr.EEG.events(:,:);
        if ~eeglab_ica_bool
            channels=sr.EEG.data(1:61,:);
            % Band-pass filtering (0.3-70Hz) all the EEG data (channels 1 to 61)
            filteredchannels=butterfilter(channels);
        else
            icaname = icanames(j+2).name;
            icafile=append(icapath2,'\',icaname);
            sr=load(icafile);
            filteredchannels=sr.filtered_data(:,:);
        end
    
         % Signal segmentation- acquiring the different trials
        
        for w = 1:length(Events(:,1))
            if Events(w,1) == movement_code
                temporal_variable=filteredchannels(:,Events(w,2)-(2.5*512):Events(w,2)+4351-(2.5*512));
                matrix(:,:,cont+1,cont2)=temporal_variable;
                events_matrix_1536(count3,:,j,i)=Events(w,:);
                count3=count3+1;
                cont=cont+1;
            end
        end
        count3=1;
    end
    cont=0;
    cont2=cont2+1;
    
end

close(findall(0,'type','figure','tag','TMWWaitbar'))


if auto_detect_movement
    [events_matrix_1536,~,~] = Optional2_dunction(subject_list,events_matrix_1536,movement_code,datapath);

end
clearvars channels sr
%%
% For each trial and EEG channel, subtract the mean value of each channel (zero-mean). 
matrix_detrend=zeros(61,4352,60,8);

if ~eeglab_ica_bool
for z = 1:length(matrix(1,1,1,:))
    for x = 1:length(matrix(:,1,1,1))
        for y = 1:length(matrix(1,1,:,1))
        temp_mean=mean(matrix(x,:,y,z));
        matrix_detrend(x,:,y,z)=(matrix(x,:,y,z))-temp_mean;
        end
    end
end
else
    matrix_detrend=matrix;
end
clearvars matrix
%%
% Here we mark the bad channels from the previous section (d1) as nan in the signal. 
for q = 1:length(bad_channels(:,1,1)) %subjects
    for i = 1:length(bad_channels(1,:,1)) %runs
        a=bad_channels(q,i,:);
        for j = a%channels
            if j ~= 80
                matrix_detrend(j,:, 1+(6*(i-1)):i*6,q)=NaN;
            end
        end
    end
end
clearvars bad_channels
%%
% Calculating kurtosis
kurt_trial_channel=ones(61,60,8);
for h=1:length(matrix_detrend(1,1,1,:))
    for z = 1: length(matrix_detrend(1,1,:,1))
        for x = 1: length(matrix_detrend(:,1,1,1))
            kurt_trial_channel(x,z,h)=kurtosis(matrix_detrend(x,:,z,h));
        end
    end
end
%%
% Bad trial obtained if is an outlier based on kurtosis
count=0;
conditions_matrix=ones(61,60,8);
for h=1:8
    for x=1:61
        countNaNs=0;
        sort_kurt=sort(kurt_trial_channel(x,:,h));
        for i = sort_kurt
            if isnan(i)
                countNaNs=countNaNs+1;
            end
        end
        corr_kurt=sort_kurt(5+countNaNs:55);
        mean_kurt=mean(corr_kurt);
        sdv_kurt=std(corr_kurt);
        high_thresh=mean_kurt+5*sdv_kurt;
        low_thresh=mean_kurt-5*sdv_kurt;
        for i=1:60
            if kurt_trial_channel(x,i,h)>high_thresh || kurt_trial_channel(x,i,h)<low_thresh
                conditions_matrix(x,i,h)=0; 
            end
            % Mark bad trial if its amplitude exceeds ±150µV
            if max(abs(matrix_detrend(x,:,i,h)),[],'omitnan')>150
                conditions_matrix(x,i,h)=0;
                count=count+1;
            end 
        end
    end
end
clearvars matrix_detrend
%%
%  mark a trial as bad if the number of free-artifact channels is lower than
% the 75% of channels
for i = 1:length(conditions_matrix(1,1,:))
    for j = 1: length(conditions_matrix(1,:,1))
        if nnz(conditions_matrix(:,j,i))<45
            conditions_matrix(:,j,i)=0;
        end
    end
end
%disp("here we still have 95% good trials")
clearvars kurt_trial_channel
%%
% Mark trials as bad if:
    % -The physical movement starts before 100ms from the appearance of the stimulus
    % -The movement starts after 2s from the appearance of the stimulus.

count4=0;
for p = 1:8
    for i= 1:6
        for j=1:10
            if (events_matrix_1536(i,3,j,p)==0) 
                conditions_matrix(:,(i-1)*10+j,p)=0;
                count4=count4+1;
            elseif((events_matrix_1536(i,3,j,p)-events_matrix_1536(i,2,j,p))<(512*0.1))
                conditions_matrix(:,(i-1)*10+j,p)=0;
                count4=count4+1;
            elseif((events_matrix_1536(i,3,j,p)-events_matrix_1536(i,2,j,p))>(512*2))
                conditions_matrix(:,(i-1)*10+j,p)=0;
                count4=count4+1;
            end
        end
    end
end
%% 


if movement_code==1536
    if eeglab_ica_bool && ~auto_detect_movement
    save('bad_trials_1536_with_ica.mat', "conditions_matrix")
    elseif eeglab_ica_bool && auto_detect_movement
    save('bad_trials_1536_with_ica_with_auto.mat', "conditions_matrix")
    elseif ~eeglab_ica_bool && auto_detect_movement
    save('bad_trials_1536_with_auto.mat', "conditions_matrix")
    else
    save('bad_trials_1536.mat', "conditions_matrix")
    end
end
if movement_code==1541
    if eeglab_ica_bool && ~auto_detect_movement
    save('bad_trials_1541_with_ica.mat', "conditions_matrix")
    elseif eeglab_ica_bool && auto_detect_movement
    save('bad_trials_1541_with_ica_with_auto.mat', "conditions_matrix")
    elseif ~eeglab_ica_bool && auto_detect_movement
    save('bad_trials_1541_with_auto.mat', "conditions_matrix")
    else
    save('bad_trials_1541.mat', "conditions_matrix")
    end
end
tEnd=toc(fstart);
disp(['Elapsed time is ' num2str(floor(tEnd/60))  ' minutes and ' num2str(rem(tEnd,60)) ' seconds']);
end