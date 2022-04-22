function create_table(channels,conditions,subject_list)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if conditions(1) ==false
    if conditions(2) == false
        to_load=["bad_trials_1536.mat","bad_trials_1541.mat"];
    else
        to_load=["bad_trials_1536_with_auto.mat","bad_trials_1541_with_auto.mat"];
    end
elseif conditions(1) ==true
    if conditions(2) == false
        to_load=["bad_trials_1536_with_ica.mat","bad_trials_1541_with_ica.mat"];
    else
        to_load=["bad_trials_1536_with_ica_with_auto.mat","bad_trials_1541_with_ica_with_auto.mat"];
    end
end
    %create a table of good trial percentage per user/movement

    bad_trials_1536=load(fullfile(to_load(1)),"conditions_matrix");
    bad_trials_1536=bad_trials_1536.conditions_matrix;

    bad_trials_1541=load(fullfile(to_load(2)),"conditions_matrix");
    bad_trials_1541=bad_trials_1541.conditions_matrix;
       


    T=[];
    for subject=1:length(bad_trials_1536(1,1,:))
            T(subject,1)=(nnz((bad_trials_1536(channels,:,subject))))/numel(bad_trials_1536(channels,:,subject));
            T(subject,2)=(nnz((bad_trials_1541(channels,:,subject))))/numel(bad_trials_1541(channels,:,subject));
    
    end
T=round(T,4)*100;
T=array2table(T,"RowNames",subject_list,"VariableNames",["Elbow Flexion", "Hand Opening"]);
disp(T)
end