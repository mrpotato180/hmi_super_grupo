movement_code=1541;
subject_list=["S1", "S2", "S3", "S6", "S7", "S9", "S11", "S13"];
conditions=[false,false];

if movement_code==1541
    if conditions(1) == true
        basefolder=dir(icapath);
        if conditions(2)==true
            condit=load("bad_trials_1541_with_ica_with_auto.mat", "conditions_matrix");
            events=load("events_matrix_1541_with_ica_with_auto.mat","events_matrix_1536");
        else
            condit=load("bad_trials_1541_with_ica.mat", "conditions_matrix");
            events=load("events_matrix_1541_with_ica.mat","events_matrix_1536");
        end
    else
        basefolder=dir(datapath);
        
        if conditions(2)==true
            condit=load("bad_trials_1541_with_auto.mat", "conditions_matrix");
            events=load("events_matrix_1541_with_auto.mat","events_matrix_1536");
        else
            condit=load("bad_trials_1541.mat", "conditions_matrix");
            events=load("events_matrix_1541.mat","events_matrix_1536");
        end
    end
else
    if conditions(1) == true
        basefolder=dir(icapath);
        if conditions(2)==true
            condit=load("bad_trials_1536_with_ica_with_auto.mat", "conditions_matrix");
            events=load("events_matrix_1536_with_ica_with_auto.mat","events_matrix_1536");
        else
            condit=load("bad_trials_1536_with_ica.mat", "conditions_matrix");
            events=load("events_matrix_1536_with_ica.mat","events_matrix_1536");
        end
    else
        basefolder=dir(datapath);
        
        if conditions(2)==true
            condit=load("bad_trials_1536_with_auto.mat", "conditions_matrix");
            events=load("events_matrix_1536_with_auto.mat","events_matrix_1536");
        else
            condit=load("bad_trials_1536.mat", "conditions_matrix");
            events=load("events_matrix_1536.mat","events_matrix_1536");
        end
    end
end
condit=condit.conditions_matrix;
events=events.events_matrix_1536;
    for i=1:length(subject_list)
        curr_folder=dir(strcat(basefolder(1).folder,'\',subject_list(i)));
        for j=1:10 %runs
            
                
            sr=load(strcat(basefolder(1).folder,'\',subject_list(i),'\',curr_folder(j+2).name), "-mat");
            
            try
                sr=sr.EEG.data(1:61,:); %if no ica is chosen
                %sr=butterfilter(sr,70);
                for ch=1:61
                sr(ch,:)=detrend(sr(ch,:),'constant');
                end
            catch
                sr=sr.filtered_data(:,:); %if ica is chosen
            end
        end
    end

conditions_matrix=condit; %61,60,8
signals_new=butterfilter(sr,10);
out_matrix=nan(60,61,2560,8);
for i=1:8
    for j=1:10
        for tr=1:6
            if conditions_matrix(:,(j-1)*6+tr,i)==1
                out_matrix((j-1)*6+tr,:,:,i)=signals_new(:,events(tr,3,j,i)-(2.5*512):events(tr,3,j,i)+(2.5*512)-1);
            end
                out_matrix((j-1)*6+tr,:,:,i)=out_matrix((j-1)*6+tr,:,:,i)-mean(out_matrix((j-1)*6+tr,:,:,i),2,'omitnan');
        end
    end
end
mean_subject_matrix=zeros(61,2560,8);
   for i=1:length(subject_list)
       mean_subject_matrix(:,:,i)=mean(out_matrix(:,:,:,i),1,'omitnan');
   end
    grand_mean_matrix=squeeze(mean(mean_subject_matrix,3));

%% D5


locs='C:\Users\Alex\Documents\UPC\Q2\HMI\EEGLAB demo-20220304\loc61eeg.loc';


[peakvalue,latency]=max(grand_mean_matrix(30,:));
time_window=[latency-51, latency+51];
potentials=mean(grand_mean_matrix(:,time_window(1):time_window(2)),2);
figure
topoplot(potentials,locs)
colorbar


