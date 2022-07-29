%% Seperate Subjects into two clusters
clearvars
clc
addpath(genpath('../matlabGiftiCifti'))

FileList_Language = dir('../Data');
% FileList_Language = FileList_Language(1:end-1);
clear files_language;
ids = [];
for i = 3:size(FileList_Language,1)
    a = FileList_Language(i).name;
    ids = [ids; a];
%     files_language((i+1)/4) = str2num(a(1:6));
end
ids = str2num(ids);
load('IQ.mat');
t = ismember(IQ_Distribution(:,1), ids);
% t(t==0) = 1;
IQ = IQ_Distribution(t==1,2:3);
rm_subs = (IQ(:,1)==0) + (IQ(:,2)==0);
rm_subs = rm_subs==2; % 2 subject is removed
IQ = IQ(~rm_subs,:);
ids = ids(~rm_subs);
type = 1; % 1: Unadj 2: Age Adjusted
IQ = IQ(:,type);

save('HighLowSubs.mat','IQ','ids','FileList_Language','rm_subs');
%% Extract Correlation Matrix
addpath(genpath('../matlabGiftiCifti'))
rs = [];
FileList_Language2 = [];
t=1;
ids2 = [];
for i = 3:size(FileList_Language,1)
    a = FileList_Language(i).name;
    if(isfile("../Data/"+a+"/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii"))
%         FileList_Language2(t) = FileList_Language(i);
        t = t + 1;
        ids2 = [ids2 str2num(a)];
    else
        i - 2
    end
end
%%
FileList_Language = FileList_Language([1:8 10:36 38:84 86:end]);
load('nifti_to_gifti.mat');
g1 = gifti('180areas.L.label.gii');
ROI_L = g1.cdata;
g2 = gifti('180areas.R.label.gii');
ROI_R = g2.cdata;
files = {'rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii', ...
         'rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii', ...
         'rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii', ...
         'rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii'};
for i = 3:size(FileList_Language,1)
    i-2
    a = FileList_Language(i).name;
    tic
    for k=1:4
        s = cifti_read("../Data/"+a+"/"+files{k});
        c = s.cdata;
        gL = zeros(32482,1200);
        gR = zeros(32482,1200);
        for j=1:size(leftcortex,1)
            gL(leftcortex(j,2)+1,:) = c(leftcortex(j,1)+1,:);
        end
        for j=1:size(rightcortex,1)
            gR(rightcortex(j,2)+1,:) = c(rightcortex(j,1)+1,:);
        end
        ROIs = zeros(360,1200);
        gL(gL == 0) = NaN;
        gR(gR == 0) = NaN;
        for j=1:180
            ROIs(j,:) = squeeze(mean(gL(ROI_L == (j + 180), :),1, 'omitnan')); 
        end
        for j=181:360
            ROIs(j,:) = squeeze(mean(gR(ROI_R == (j - 180), :),1, 'omitnan')); 
        end
        for j=1:6
            c2 = ROIs(:, (j-1)*200+1:j*200);
            [r, p] = corrcoef(c2');
            z=.5.*log((1+r)./(1-r));
            z = z';
            rs(i-2,k,j,:) = z(tril(true(size(r)),-1));
        end
    end
    toc
end
IQ = IQ(sum(ids==ids2,2)==1);
save('Features.mat');
%% Load data and IQ
clc
clearvars
load('Features.mat');

IQs = zeros(83,4,6);
for i=1:83
    IQs(i,:,:) = IQ(i);
end

IQs = reshape(IQs, [83*4*6 1]);
data = reshape(rs, [83*4*6 size(rs,4)]);

t0 = t - std(IQ(IQ(:,type)>0,1));
t1 = t + std(IQ(IQ(:,type)>0,1));

HIQ = IQs>median(IQ);
LIQ = IQs<=median(IQ);
IQ = IQs > median(IQ);

%% Train model SVM
[coeff,score,latent] = pca(data);
dim_vec = linspace(1,size(coeff,2),10);

ratio = 0.2;
rng('shuffle') 
Accs = [];
for i=1:length(dim_vec)
    i
    data_reduced = data * coeff(:,1:dim_vec(i));
    t = 1:size(data_reduced,1);
    sel_subs_test = randperm(size(data_reduced,1),ceil(size(data_reduced,1)*ratio));
    sel_subs_train = t(setdiff(1:end,sel_subs_test));
    d_SVM_test = data_reduced(sel_subs_test,:);
    d_SVM_train = data_reduced(sel_subs_train,:);
    SVMModel = fitcsvm(d_SVM_train, IQ(sel_subs_train), 'KernelFunction','rbf');
    CVSVMModel = crossval(SVMModel);
    Acc = 0;
    for j=1:CVSVMModel.KFold
        predicted = predict(CVSVMModel.Trained{j},d_SVM_test);
        Acc = Acc + (1-sum(abs(predicted-IQ(sel_subs_test)))/size(sel_subs_train,2))*100;
    end
    Accs(i) = Acc / CVSVMModel.KFold;
end

plot(dim_vec,Accs)
%% Train Model KNN
[coeff,score,latent] = pca(data);
dim_vec = linspace(1,size(coeff,2),10);
% dim_vec = linspace(1,50,10);

ratio = 0.2;
rng('shuffle') 
Accs = [];
for i=1:length(dim_vec)
    i
    data_reduced = data * coeff(:,1:dim_vec(i));
    t = 1:size(data_reduced,1);
    sel_subs_test = randperm(size(data_reduced,1),ceil(size(data_reduced,1)*ratio));
    sel_subs_train = t(setdiff(1:end,sel_subs_test));
    d_SVM_test = data_reduced(sel_subs_test,:);
    d_SVM_train = data_reduced(sel_subs_train,:);
    KNNMdl = fitcknn(d_SVM_train,IQ(sel_subs_train),...
                        'OptimizeHyperparameters','auto');
    CVSVMModel = crossval(KNNMdl);
    Acc = 0;
    for j=1:CVSVMModel.KFold
        predicted = predict(CVSVMModel.Trained{j},d_SVM_test);
        Acc = Acc + (1-sum(abs(predicted-IQ(sel_subs_test)))/size(sel_subs_train,2))*100;
    end
    Accs(i) = Acc / CVSVMModel.KFold;
end

plot(dim_vec,Accs)
%% Train Model KNN ==> 50
[coeff,score,latent] = pca(data);
% dim_vec = linspace(1,size(coeff,2),10);
dim_vec = [200];

ratio = 0.2;
rng('shuffle') 
Accs = [];
for i=1:length(dim_vec)
    data_reduced = data * coeff(:,1:dim_vec(i));
    t = 1:size(data_reduced,1);
    sel_subs_test = randperm(size(data_reduced,1),ceil(size(data_reduced,1)*ratio));
    sel_subs_train = t(setdiff(1:end,sel_subs_test));
    d_SVM_test = data_reduced(sel_subs_test,:);
    d_SVM_train = data_reduced(sel_subs_train,:);
    KNNMdl = fitcknn(d_SVM_train,IQ(sel_subs_train),'NumNeighbors',3);
    CVSVMModel = crossval(KNNMdl);
    Acc = 0;
    for j=1:CVSVMModel.KFold
        predicted = predict(CVSVMModel.Trained{j},d_SVM_test);
        Acc = Acc + (1-sum(abs(predicted-IQ(sel_subs_test)))/size(sel_subs_train,2))*100;
    end
    Accs(i) = Acc / CVSVMModel.KFold
end

plot(dim_vec,Accs)