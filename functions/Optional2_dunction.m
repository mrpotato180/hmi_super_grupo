function [samples_of_movement_1536,residous,chanels_93] = Optional2_dunction(subject_list,events_matrix_1536,event_code,datapath)

if event_code == 1536
    channel_number=93;
    thres_variable=8;
    thres_variable2=0.045;
elseif event_code == 1541
    channel_number=69;
    thres_variable=-0.028;
    thres_variable2=0.005;    
end



 
chanels_93=zeros(10,200000,8);

for i=1:length(subject_list)
    folder_path=datapath;
    curr_folder=append(folder_path,filesep,subject_list(i));
    filenames=dir(curr_folder);
    for j=1:10
        filename = filenames(j+2).name;
        
        path=append(curr_folder,filesep,filename);
        sr=load(fullfile(path));
        chanels_93(j,1:length(sr.EEG.data(1,:)),i)=sr.EEG.data(channel_number,:);
    end
end
%events_matrix_1536=zeros(6,3,10,8);

samples_of_movement_1536=zeros(6,4,10,8);
samples_of_movement_1536(:,1:2,:,:)=events_matrix_1536(:,1:2,:,:);
chanels_93_derivative=zeros(10,199999,8);

for t=1:8
    for y=1:10
        chanels_93_derivative(y,:,t)=diff(chanels_93(y,:,t),1,2);
    end
end

for x = 1:length(chanels_93(1,1,:)) %8
    for y = 1:10
        count2=0;
        if event_code == 1536
            for i = 1:length(events_matrix_1536(:,2,y,x)) %6
                if max(chanels_93(y,events_matrix_1536(i,2,y,x):events_matrix_1536(i,2,y,x)+2850,x)) - chanels_93(y,events_matrix_1536(i,2,y,x),x) > thres_variable
                    while abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count2,x))< thres_variable2
                        count2=count2+1;
                    end
                    samples_of_movement_1536(i,3,y,x)=events_matrix_1536(i,2,y,x)+count2;
                    [localmax,locallatency]=max(chanels_93(y,samples_of_movement_1536(i,3,y,x):samples_of_movement_1536(i,3,y,x)+1800,x));
                    count3=count2+locallatency;
                    count2=0;
                    while abs(chanels_93(y,events_matrix_1536(i,2,y,x)+count3,x))>0.5*localmax && count3 < 20000
                        count3=count3+1;
                    end
                    while (abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count3,x))> thres_variable2 || abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count3,x))== 0) && events_matrix_1536(i,2,y,x)+ count3 <199998  
                        count3=count3+1;
                    end
                    if count3 >= 15000
                        samples_of_movement_1536(i,4,y,x)=0;
                    else
                        samples_of_movement_1536(i,4,y,x)=events_matrix_1536(i,2,y,x)+count3;
                    end

                else
                    samples_of_movement_1536(i,3,y,x)=0;
                end
            end
        elseif event_code == 1541
            for i = 1:length(events_matrix_1536(:,2,y,x)) %6
                if min(chanels_93(y,events_matrix_1536(i,2,y,x):events_matrix_1536(i,2,y,x)+2850,x))- chanels_93(y,events_matrix_1536(i,2,y,x),x) < thres_variable
                    while abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count2,x))< thres_variable2
                        count2=count2+1;
                    end
                    samples_of_movement_1536(i,3,y,x)=events_matrix_1536(i,2,y,x)+count2;
                    [localmax,locallatency]=min(chanels_93(y,samples_of_movement_1536(i,3,y,x):samples_of_movement_1536(i,3,y,x)+1300,x));
                    count3=count2+locallatency;
                    count2=0;
                    while abs(chanels_93(y,events_matrix_1536(i,2,y,x)+count3,x))>0.5*abs(localmax-chanels_93(y,samples_of_movement_1536(i,3,y,x),x)) && count3 < 20000
                        count3=count3+1;
                    end
                    while (abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count3,x))> thres_variable2 || abs(chanels_93_derivative(y,events_matrix_1536(i,2,y,x)+count3,x))== 0) && events_matrix_1536(i,2,y,x)+ count3 <199998  
                        count3=count3+1;
                    end
                    if count3 >= 15000
                        samples_of_movement_1536(i,4,y,x)=0;
                    else
                        samples_of_movement_1536(i,4,y,x)=events_matrix_1536(i,2,y,x)+count3;
                    end
                else
                    samples_of_movement_1536(i,3,y,x)=0;
                end
            end
        end
    end
end

residous=zeros(6,10,8);
for j = 1:8
    for i = 1:10
        for x = 1:6
            residous(x,i,j)=samples_of_movement_1536(x,3,i,j)-events_matrix_1536(x,3,i,j);
        end
    end
end
%si se hace plot de los casos que salen raros, se observa congruencia con
%el grafico representado. 

%testing parameter
testingresidous=residous;
for j = 1:8
    for i = 1:10
        for x = 1:6
            if abs(residous(x,i,j))>1000
                testingresidous(x,i,j)=0;
            end
        end
    end
end
M4=mean(mean(mean(testingresidous)));
M5 = ['The average difference with respect to the provided events is: ',num2str(M4/512),'s'];
disp(M5)

