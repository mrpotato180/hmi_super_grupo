function [grand_mean_mu,grand_mean_beta,ERD_mu,ERD_beta,grand_mean_mu_ERP,grand_mean_beta_ERP,ERP_mu,ERP_beta,ERD_mu_pj,ERD_beta_pj,ERP_mu_pj,ERP_beta_pj] = D6_function(movement_code,conditions,subject_list,datapath,icapath, ERP, project)

if movement_code==1541
    if conditions(1) == true
        basefolder=dir(icapath);
        if conditions(2)==true
            condit=load("bad_trials_1541_with_ica_with_auto.mat", "conditions_matrix");
            events=load("events_matrix_1541_with_auto.mat","events_matrix_1536");
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
            events=load("events_matrix_1536_with_auto.mat","events_matrix_1536");
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
out_matrix_mu_ERP=nan(60,61,6144,8);
out_matrix_beta_ERP=nan(60,61,6144,8);

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
                        if ERP
                            if events(tr,4,j,i) ~= 0
                                out_matrix_mu_ERP((j-1)*6+tr,ch,:,i)=mu_signals_new(ch,events(tr,4,j,i)-(8.5*512):events(tr,4,j,i)+(3.5*512)-1);
                                out_matrix_beta_ERP((j-1)*6+tr,ch,:,i)=beta_signals_new(ch,events(tr,4,j,i)-(8.5*512):events(tr,4,j,i)+(3.5*512)-1);
                            end
                        end
                    end
                end
            end
        end
    end
end
clearvars sr mu_signals_new beta_signals_new events
%%

for i=1:8
    for tr=1:60
        for ch=1:61
                out_matrix_mu(tr,ch,:,i)=out_matrix_mu(tr,ch,:,i)-mean(out_matrix_mu(tr,:,:,i),2,'omitnan');
                out_matrix_beta(tr,ch,:,i)=out_matrix_beta(tr,ch,:,i)-mean(out_matrix_beta(tr,:,:,i),2,'omitnan');

        end
    end
end

%6145 length
if ERP
    musignals_ERP=out_matrix_mu_ERP;
    betasignals_ERP=out_matrix_beta_ERP;
    for i=1:8
        for tr=1:60
            for ch=1:61
                musignals_ERP(tr,ch,:,i)=out_matrix_mu_ERP(tr,ch,:,i)-mean(out_matrix_mu_ERP(tr,:,:,i),2,'omitnan');
                betasignals_ERP(tr,ch,:,i)=out_matrix_beta_ERP(tr,ch,:,i)-mean(out_matrix_beta_ERP(tr,:,:,i),2,'omitnan');

            end
        end
    end
end


%% will the real d6 please stand up

muenergy=out_matrix_mu.^2;
betaenergy=out_matrix_beta.^2;

clearvars out_matrix_mu out_matrix_beta

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

if ERP
    muenergy_ERP=musignals_ERP.^2;
    betaenergy_ERP=betasignals_ERP.^2;

    mu_av_energy_ERP=squeeze(mean(muenergy_ERP,1,'omitnan'));
    beta_av_energy_ERP=squeeze(mean(betaenergy_ERP,1,'omitnan'));

    mu_mov_av_ERP=movmean(mu_av_energy_ERP,round(0.15*512),2);
    beta_mov_av_ERP=movmean(beta_av_energy_ERP,round(0.15*512),2);

    ERP_mu=nan(61,6144,8);
    ERP_beta=nan(61,6144,8);

    for subject=1:8
        for ch=1:61

            ERP_mu(ch,:,subject)=100*(mu_mov_av_ERP(ch,:,subject)-mu_ref(ch,subject))/mu_ref(ch,subject);
            ERP_beta(ch,:,subject)=100*(beta_mov_av_ERP(ch,:,subject)-beta_ref(ch,subject))/beta_ref(ch,subject);
        end
    end


    grand_mean_mu_ERP=squeeze(mean(ERP_mu,3));
    grand_mean_beta_ERP=squeeze(mean(ERP_beta,3));
end

%extra for the project of dades
if project
    ERD_mu_pj=nan(60,61,4096,8);
    ERD_beta_pj=nan(60,61,4096,8);
    
    mu_mov_av_pj=movmean(muenergy,round(0.15*512),3);
    beta_mov_av_pj=movmean(betaenergy,round(0.15*512),3);
    
    for subject=1:8
        for tr=1:60
            for ch=1:61
    
                ERD_mu_pj(tr,ch,:,subject)=100*(mu_mov_av_pj(tr,ch,:,subject)-mu_ref(ch,subject))/mu_ref(ch,subject);
                ERD_beta_pj(tr,ch,:,subject)=100*(beta_mov_av_pj(tr,ch,:,subject)-beta_ref(ch,subject))/beta_ref(ch,subject);
            end
        end
    end
    
    ERP_mu_pj=nan(60,61,6144,8);
    ERP_beta_pj=nan(60,61,6144,8);
    
    mu_mov_av_pj=movmean(muenergy_ERP,round(0.15*512),3);
    beta_mov_av_pj=movmean(betaenergy_ERP,round(0.15*512),3);
    
    for subject=1:8
        for tr=1:60
            for ch=1:61
    
                ERP_mu_pj(tr,ch,:,subject)=100*(mu_mov_av_pj(tr,ch,:,subject)-mu_ref(ch,subject))/mu_ref(ch,subject);
                ERP_beta_pj(tr,ch,:,subject)=100*(beta_mov_av_pj(tr,ch,:,subject)-beta_ref(ch,subject))/beta_ref(ch,subject);
            end
        end
    end
end

end

