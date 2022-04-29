
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
conditions_matrix=condit; %61,60,8
events=events.events_matrix_1536;
out_matrix_mu=nan(60,61,4096,8);
out_matrix_beta=nan(60,61,4096,8);

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



        mu_signals_new=bandpassfilter(sr,8,13);
        beta_signals_new=bandpassfilter(sr,14,30);

        for tr=1:6
            for ch=1:61
                if conditions_matrix(ch,(j-1)*6+tr,i)==1

                    if events(tr,3,j,i) ~= 0
                        out_matrix_mu((j-1)*6+tr,ch,:,i)=mu_signals_new(ch,events(tr,3,j,i)-(2.5*512):events(tr,3,j,i)+(5.5*512)-1);
                        out_matrix_beta((j-1)*6+tr,ch,:,i)=beta_signals_new(ch,events(tr,3,j,i)-(2.5*512):events(tr,3,j,i)+(5.5*512)-1);

                    end
                end
            end

        end
    end
end
%%
musignals=out_matrix_mu;
betasignals=out_matrix_beta;
for i=1:8
    for tr=1:60
        for ch=1:61
                musignals(tr,ch,:,i)=out_matrix_mu(tr,ch,:,i)-mean(out_matrix_mu(tr,:,:,i),2,'omitnan');
                betasignals(tr,ch,:,i)=out_matrix_beta(tr,ch,:,i)-mean(out_matrix_beta(tr,:,:,i),2,'omitnan');

        end
    end
end
clearvars out_matrix
%% will the real d6 please stand up

muenergy=musignals.^2;
betaenergy=betasignals.^2;

mu_av_energy=squeeze(mean(muenergy,1,'omitnan'));
beta_av_energy=squeeze(mean(betaenergy,1,'omitnan'));

mu_mov_av=movmean(mu_av_energy,round(0.15*512),2);
beta_mov_av=movmean(beta_av_energy,round(0.15*512),2);

mu_ref=squeeze(mean(mu_mov_av(:,1*512:1.5*512,:),2,'omitnan'));
beta_ref=squeeze(mean(beta_mov_av(:,1*512:1.5*512,:),2,'omitnan'));

ERD_mu=nan(61,4096,8);
ERD_beta=nan(61,4096,8);

for subject=1:8
    for ch=1:61

        ERD_mu(ch,:,subject)=100*(mu_mov_av(ch,:,subject)-mu_ref(ch,subject))/mu_ref(ch,subject);
        ERD_beta(ch,:,subject)=100*(beta_mov_av(ch,:,subject)-beta_ref(ch,subject))/beta_ref(ch,subject);
    end
end


grand_mean_mu=squeeze(mean(ERD_mu,3));
grand_mean_beta=squeeze(mean(ERD_beta,3));


%% D7
grand_mean_matrix=grand_mean_mu;
eeglab;close;
locs='C:\Users\Alex\Documents\UPC\Q2\HMI\EEGLAB demo-20220304\loc61eeg.loc';
potentials=zeros(61,1);
for ch=1:61
    potentials(ch,1)=min(grand_mean_matrix(ch,2.4*512:2.6*512));
end
figure
topoplot(squeeze(potentials),locs,'maplimits',[-60 0]);

colormap turbo;
colorbar;

  
