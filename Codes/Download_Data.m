%% Download from Amazon S3
clearvars
clc
% load('Subjects.mat');
load('ids.mat');

% Subjects = Subjects(Subjects(:,2)==1,1);
Subjects = ids;
Bucket = 'hcp-openaccess';
err = [];

for i=1:length(Subjects)
    % Key: HCP_1200/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii
%     file = {['Data/' num2str(Subjects(i)) '/rmegpreproc/' num2str(Subjects(i)) '_MEG_3-Restin_rmegpreproc.mat'],...
%             ['Data/' num2str(Subjects(i)) '/rmegpreproc/' num2str(Subjects(i)) '_MEG_4-Restin_rmegpreproc.mat'],...
%             ['Data/' num2str(Subjects(i)) '/rmegpreproc/' num2str(Subjects(i)) '_MEG_5-Restin_rmegpreproc.mat']};
%     url = {['/MEG/Restin/rmegpreproc/' num2str(Subjects(i)) '_MEG_3-Restin_rmegpreproc.mat '],...
%             ['/MEG/Restin/rmegpreproc/' num2str(Subjects(i)) '_MEG_4-Restin_rmegpreproc.mat '],...
%             ['/MEG/Restin/rmegpreproc/' num2str(Subjects(i)) '_MEG_5-Restin_rmegpreproc.mat ']};
    file = {['Data/' num2str(Subjects(i)) '/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'],...
            ['Data/' num2str(Subjects(i)) '/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'],...
            ['Data/' num2str(Subjects(i)) '/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii'],...
            ['Data/' num2str(Subjects(i)) '/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii']};
    url = {['/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii '],...
           ['/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii '],...
           ['/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii '],...
           ['/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii ']};
    if(sum(isfile(file))~=4)
        tic
        for j=1:length(file)
            if(~isfile(file{j}))
                path = ['s3://hcp-openaccess/HCP_1200/' num2str(Subjects(i)) url{j} file{j}];
                [status,cmdout] = system(['aws s3 cp ' path]);
                if(status==0)
                    fprintf([num2str(i) '.' num2str(j) ' Data of subject ' num2str(Subjects(i)) ': Ok! \n']);
                else
                    fprintf([num2str(i) '.' num2str(j) ' Data of subject ' num2str(Subjects(i)) ': Error! \n']);
                    err = [err Subjects(i)];
                end
            else
                fprintf([num2str(i) '.' num2str(j) ' Data exist for subject ' num2str(Subjects(i)) ' \n']);
            end
        end
        toc
    else
        fprintf([num2str(i) '. Data exist for subject ' num2str(Subjects(i)) ' \n']);
    end
end