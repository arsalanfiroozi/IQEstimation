%% Seperate Subjects into two clusters
clearvars
addpath(genpath('../../matlabGiftiCifti'))

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

% load('Subjects.mat');
% Subjects = Subjects(Subjects(:,2)==1,1);
% ids = Subjects;
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
%% Find Common Channels
clc
clearvars
load('HighLowSubs.mat');

load("./Data/"+FileList_Language(3).name+"/rmegpreproc/"+FileList_Language(3).name+"_MEG_3-Restin_rmegpreproc.mat");
comm_chan = data.label;
for i = 3:size(FileList_Language,1)
    a = FileList_Language(i).name;
    for j=3:5
        tic
        load("./Data/"+a+"/rmegpreproc/"+a+"_MEG_"+j+"-Restin_rmegpreproc.mat");
        toc
        labels = data.label;
        comm_chan = comm_chan(ismember(comm_chan,labels));
        j
    end
    i
end

save('CommonChannels.mat','comm_chan');
%% Extract Power measures based on Common Channels
clc 
clearvars
load('HighLowSubs.mat')
load('CommonChannels.mat')

% Central, Left Temporal, Right Temporal, Anterior, Posterior
sel_channels = {'A23' 'A33' 'A156' 'A181' 'A174' 'A142' 'A124' 'A158' 'A218' 'A242'};
comm_chan = {'A23' 'A33' 'A156' 'A181' 'A174' 'A142' 'A124' 'A158' 'A218' 'A242'};

Bounds = [[0 4];[4 7];[8 15];[16 31];[32 Inf];[8 12]];
Bounds_Label = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Mu'};

% Extract Power measures based on Common Channels
[comm_chan, ~] = sort(comm_chan);
sel_indices = [];
for i=1:length(sel_channels)
    ind = find(strcmp(comm_chan,sel_channels{i}));
    sel_indices = [sel_indices ind];
end

dataset = [];
dataset_coherence = [];
for i = 3:size(FileList_Language,1)
    tic
    a = FileList_Language(i).name;
    for j=3:5
        load("..\..\data\"+a+"\rmegpreproc\"+a+"_MEG_"+j+"-Restin_rmegpreproc.mat");
        d_trials = [];
        [labels, inds] = sort(data.label);
        for k=1:length(data.trial)
            t = data.trial{1,k};
            t = t(inds,:);
            d_trials = t;
            xd = d_trials(ismember(labels,comm_chan),:);
            Fs = data.fsample;
            for tt=1:size(xd,1)
                x = xd(tt,:);
                N = length(x);
                xdft = fft(x);
                xdft = xdft(1:N/2+1);
                psdx = (1/(Fs*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                freq = 0:Fs/length(x):Fs/2;
                % Power
                for r=1:length(Bounds)
                    sel_f = freq > Bounds(r,1) & freq < Bounds(r,2);
                    dataset(i-2,j-2,k,tt,r) = sum(psdx(sel_f));
                end
                % Normalized Power
                for r=1:length(Bounds)
                    sel_f = freq > Bounds(r,1) & freq < Bounds(r,2);
                    dataset(i-2,j-2,k,tt,r+length(Bounds)) = sum(psdx(sel_f))/sum(psdx);
                end
                % Relative Power
                t = 1;
                for k1=1:length(Bounds)
                    for k2=1:k1-1
                        sel_f1 = freq > Bounds(k1,1) & freq < Bounds(k1,2);
                        sel_f2 = freq > Bounds(k2,1) & freq < Bounds(k2,2);
                        dataset(i-2,j-2,k,tt,2*length(Bounds)+t) = sum(psdx(sel_f1))/sum(psdx(sel_f2));
                        t = t + 1;
                    end
                end        
            end
            % Coherency
            t = 1;
            xdn = xd(sel_indices,:);
            for k1=1:size(xdn,1)
                for k2=1:k1-1
                    x1 = xdn(k1,:);
                    x2 = xdn(k2,:);
                    [c, ft] = mscohere(x1,x2,[],[],[],Fs);
                    for r=1:length(Bounds)
                        sel_f = ft > Bounds(r,1) & ft < Bounds(r,2);
                        dataset_coherence(i-2,j-2,k,t,r) = sum(c(sel_f));
                    end
                    t = t + 1;
                end
            end
        end
        j
    end
    i
    toc
end

save('ExtractedFeatures.mat','dataset','dataset_coherence');
%% Normalize Model for Training
clc 
clearvars
load('HighLowSubs.mat')
load('CommonChannels.mat')
load('ExtractedFeatures.mat')
addpath(genpath('matlabGiftiCifti'))
N_Features = 27;
iqm = 2;

% Power Normalizing
dataset = dataset(~rm_subs, :, :, :, :);
dataset_r = reshape(dataset, [size(dataset,2)*size(dataset,1)*size(dataset,3) size(dataset,4) size(dataset,5)]);
dataset_norm = reshape(dataset, [size(dataset,2)*size(dataset,1)*size(dataset,3)*size(dataset,4) size(dataset,5)]);
% Normalizing Features
for i=1:N_Features
    dataset_r(:,:,i) = dataset_r(:,:,i) ./ max(dataset_norm(:,i));
end

% Coherency Normalizing
dataset_coherence = dataset_coherence(~rm_subs, :, :, :, :);
dataset_coherence_r = reshape(dataset_coherence, [size(dataset_coherence,2)*size(dataset_coherence,1)*size(dataset_coherence,3) size(dataset_coherence,4) 6]);
dataset_coherence_norm = reshape(dataset_coherence, [size(dataset_coherence,2)*size(dataset_coherence,1)*size(dataset_coherence,3)*size(dataset_coherence,4) 6]);
for i=1:6
    dataset_coherence_r(:,:,i) = dataset_coherence_r(:,:,i) ./ max(dataset_coherence_norm(:,i));
end

pvals = [];
d_SVM = [];
for i=1:size(dataset,4)
    for j=1:N_Features
        d = dataset_r(:,i,j);
        
        d_SVM = [d_SVM d];
%         [h,p] = ttest(d, IQt);
%         pvals(i,j) = p;
    end
end

pvals_coherence = [];
for i=1:size(dataset_coherence,4)
    for j=1:6
        d = dataset_coherence_r(:,i,j);
        
        d_SVM = [d_SVM d];
%         [h,p] = ttest(d, IQt);
%         pvals_coherence(i,j) = p;
    end
end
% sum(sum(pvals < 0.05)) + sum(sum(pvals_coherence < 0.05))

%% Regression
X = [ones(size(d_SVM(:,1))) d_SVM];
IQs = ones([86 3 149])*10;
for i=1:86
    IQs(i,:,:) = IQ(i);
end
IQs = reshape(IQs, [size(IQs,1)*size(IQs,2)*size(IQs,3) 1]);
ratio = 0.2;
rng('shuffle') 
t = 1:size(X,1);
sel_subs_test = randperm(size(X,1),ceil(size(X,1)*ratio));
sel_subs_train = t(setdiff(1:end,sel_subs_test));
x_NN_test = X(sel_subs_test,:);
x_NN_train = X(sel_subs_train,:);
y_NN_test = IQs(sel_subs_test);
y_NN_train = IQs(sel_subs_train);
%% Linear Regression Model
[b,bint,r,rint,stats] = regress(y_NN_train,x_NN_train);
y_test_predicted = x_NN_test * b;
[r, p] = corrcoef(y_test_predicted,y_NN_test)
[~, s] = sort(y_test_predicted);
scatter(y_NN_test(s),y_test_predicted(s));
xlabel('IQ')
ylabel('Predicted')
% [x, y] = sort(y_test_predicted);
% scatter(x, y_NN_test(y))
% err = sqrt(sum((y_test_predicted-y_NN_test).^2)./length(y_NN_test))
%% GLM
b = glmfit(x_NN_train(:,1:end-1),y_NN_train);
y_test_predicted = x_NN_test * b;
[r, p] = corrcoef(y_test_predicted,y_NN_test)
% err = sqrt(sum((y_test_predicted-y_NN_test).^2)./length(y_NN_test))
%% SVM
clc
Mdl = fitrsvm(x_NN_train(:,1:end-1),y_NN_train);
y_test_predicted = predict(Mdl,x_NN_test(:,1:end-1));
% tic
% CVSVMModel = crossval(Mdl);
% toc
% t = 0;
% for i=1:CVSVMModel.KFold
%     predicted = predict(CVSVMModel.Trained{i},x_NN_test(:,1:end-1));
% %     t = t + predict(Mdl,x_NN_test(:,1:end-1));
%     t = t + predicted;
% end
% y_test_predicted = t ./ CVSVMModel.KFold;
[r, p] = corrcoef(y_test_predicted,y_NN_test)
% err = sqrt(sum((y_test_predicted-y_NN_test).^2)./length(y_NN_test))
%% Gaussian Process Noise
Mdl = fitrgp(x_NN_train(:,1:end-1),y_NN_train);
y_test_predicted = predict(Mdl,x_NN_test(:,1:end-1));
% CVSVMModel = crossval(Mdl);
% t = 0;
% for i=1:CVSVMModel.KFold
%     predicted = predict(CVSVMModel.Trained{i},x_NN_test(:,1:end-1));
%     t = t + predicted;
% end
% y_test_predicted = t ./ CVSVMModel.KFold;
% err = sqrt(sum((y_test_predicted-y_NN_test).^2)./length(y_NN_test))
[r, p] = corrcoef(y_test_predicted,y_NN_test)

save('Gaussian.mat');