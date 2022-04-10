function [grand_mean_matrix,mean_subject_matrix] = mean_all_files(movement_code,conditions,subject_list,datapath,icapath)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
out_matrix=nan(61,60,4352,8);

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
                sr=sr.EEG.data(1:61,:);
                sr=butterfilter(sr);
                for ch=1:61
                sr(ch,:)=detrend(sr(ch,:),'constant');
                end
            catch
                sr=sr.filtered_data(:,:);
            end
            
            for tr=1:6
                if condit(:,(j-1)*6+tr,i)==1
                    out_matrix(:,(j-1)*6+tr,:,i)=sr(:,events(tr,2,j,i)-(2.5*512):events(tr,2,j,i)+(6*512)-1);
                end
            end
                
            
        end
    end

   mean_subject_matrix=zeros(61,4352,8);
   for i=1:length(subject_list)
       mean_subject_matrix(:,:,i)=mean(out_matrix(:,:,:,i),2,'omitnan');
   end
    grand_mean_matrix=mean(mean_subject_matrix,3);
end