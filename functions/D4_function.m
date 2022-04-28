function [grand_mean_matrix,mean_subject_matrix] = D4_function(movement_code,conditions,subject_list,datapath,icapath)


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
conditions_matrix=condit; %61,60,8
events=events.events_matrix_1536;
out_matrix=nan(60,61,2560,8);

for i=1:length(subject_list)

    for j=1:10 %runs
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
        if conditions(1)
            icaname=strcat('filtered_data_',string(j));
        end
        filename = strcat('ME_S',string(subject_code),'_r',string(run_code));

        if ~conditions(1) % if no ica
            sr=load(strcat(fullfile(basefolder(1).folder),filesep,subject_list(i),filesep,filename), "-mat");


            sr=sr.EEG.data(1:61,:); %if no ica is chosen
            %sr=butterfilter(sr,70);
            %                 for ch=1:61
            %                 sr(ch,:)=detrend(sr(ch,:),'constant');
            %                 end
        else
            sr=load(strcat(fullfile(basefolder(1).folder),filesep,subject_list(i),filesep,icaname), "-mat");
            sr=sr.filtered_data(:,:); %if ica is chosen
        end



        signals_new=butterfilter(sr,10);

        for tr=1:6
            for ch=1:61
                if conditions_matrix(ch,(j-1)*6+tr,i)==1

                    if events(tr,3,j,i) ~= 0
                        out_matrix((j-1)*6+tr,ch,:,i)=signals_new(ch,events(tr,3,j,i)-(2.5*512):events(tr,3,j,i)+(2.5*512)-1);

                    end
                end
            end

        end
    end
end
%%
johnny_walker=out_matrix;
for i=1:8
    for tr=1:60
        for ch=1:61
                johnny_walker(tr,ch,:,i)=out_matrix(tr,ch,:,i)-mean(out_matrix(tr,:,:,i),2,'omitnan');

        end
    end
end
mean_subject_matrix=zeros(61,2560,8);
   for i=1:length(subject_list)
       mean_subject_matrix(:,:,i)=mean(johnny_walker(:,:,:,i),1,'omitnan');
   end
    grand_mean_matrix=squeeze(mean(mean_subject_matrix,3));
end
