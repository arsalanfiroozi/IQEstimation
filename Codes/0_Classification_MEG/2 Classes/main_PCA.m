%% Seperate Subjects into two clusters
clearvars
addpath(genpath('../../matlabGiftiCifti'))

FileList_Language = dir('../../Data');
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
t = median(IQ(IQ(:,type)>0,1));
HIQ = IQ(:,type)>t;
sum(HIQ)
LIQ = IQ(:,type)<=t;
sum(LIQ)

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
% histogram(IQ,15,'FaceColor','blue');
box off
title('IQ Distribution of MEG Subjects');
histogram(IQ(LIQ),30,'FaceColor','blue');
hold on
box off
histogram(IQ(HIQ),30,'FaceColor','red');
legend({'Low IQ', 'High IQ'});
export_fig('3.png','-r600');

% save('HighLowSubs.mat','HIQ','LIQ','ids','FileList_Language','rm_subs','IQ');
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
%% Freq Bounds Spec & Extract Features -- Approximate Running Time: 6h
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
%% Normalize for Model Training
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
        IQt = [IQ(:,iqm) IQ(:,iqm) IQ(:,iqm)]; 
%         IQt = reshape(IQt, [size(IQt,1)*size(IQt,2) 1]);
        
        d_SVM = [d_SVM d];
%         [h,p] = ttest(d, IQt);
%         pvals(i,j) = p;
    end
end

pvals_coherence = [];
for i=1:size(dataset_coherence,4)
    for j=1:6
        d = dataset_coherence_r(:,i,j);
%         IQt = [IQ(:,iqm) IQ(:,iqm) IQ(:,iqm)]; 
%         IQt = reshape(IQt, [size(IQt,1)*size(IQt,2) 1]);
        
        d_SVM = [d_SVM d];
%         [h,p] = ttest(d, IQt);
%         pvals_coherence(i,j) = p;
    end
end
% sum(sum(pvals < 0.05)) + sum(sum(pvals_coherence < 0.05))
%% SVM
clc
[coeff,score,latent] = pca(d_SVM);
Accs = [];
dim_vec = linspace(1,size(coeff,2),10);
j = 2;
d_SVM_Reduced = d_SVM * coeff(:,1:dim_vec(j));
h = scatter3(d_SVM_Reduced(IQ_Labels==1,1),d_SVM_Reduced(IQ_Labels==1,2),d_SVM_Reduced(IQ_Labels==1,3),'filled')
h.SizeData = 10; hold on
h = scatter3(d_SVM_Reduced(IQ_Labels==0,1),d_SVM_Reduced(IQ_Labels==0,2),d_SVM_Reduced(IQ_Labels==0,3),'filled')
h.SizeData = 10;

% 20% Test, 80% Train
ratio = 0.2;
rng('shuffle') 
t = 1:size(d_SVM_Reduced,1);
sel_subs_test = randperm(size(d_SVM_Reduced,1),ceil(size(d_SVM_Reduced,1)*ratio));
sel_subs_train = t(setdiff(1:end,sel_subs_test));
d_SVM_test = d_SVM_Reduced(sel_subs_test,:);
d_SVM_train = d_SVM_Reduced(sel_subs_train,:);
% IQ_Labels = IQt>median(IQ(IQ(:,iqm)>0,1));
%         IQ_Labels = [HIQ HIQ HIQ];
%         IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2) 1]);
IQ_Labels = ones([86 3 149])*10;
for i=1:86
    IQ_Labels(i,:,:) = HIQ(i);
end
IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2)*size(IQ_Labels,3) 1]);

SVMModel = fitcsvm(d_SVM_train, IQ_Labels(sel_subs_train),'KernelFunction','polynomial');
2
%                         'OptimizeHyperparameters','auto');
%                         'OptimizeHyperparameters','auto',...
%                         'HyperparameterOptimizationOptions',struct('UseParallel',true));
CVSVMModel = crossval(SVMModel);
Acc = 0;
for i=1:CVSVMModel.KFold
    predicted = predict(CVSVMModel.Trained{i},d_SVM_test);
    Acc = Acc + (1-sum(abs(predicted-IQ_Labels(sel_subs_test)))/size(sel_subs_train,2))*100;
end
Acc = Acc / CVSVMModel.KFold
% 90.1073 Linear
% 93.3662 RBF
% 92.0798 polynomial
%% SVM
[coeff,score,latent] = pca(d_SVM);
Accs = [];
dim_vec = linspace(1,size(coeff,2),20);
for k=1:1
    k
    for j=1:length(dim_vec)
        tic
        j
        d_SVM_Reduced = d_SVM * coeff(:,1:dim_vec(j));

        % 20% Test, 80% Train
        ratio = 0.2;
        rng('shuffle') 
        t = 1:size(d_SVM_Reduced,1);
        sel_subs_test = randperm(size(d_SVM_Reduced,1),ceil(size(d_SVM_Reduced,1)*ratio));
        sel_subs_train = t(setdiff(1:end,sel_subs_test));
        d_SVM_test = d_SVM_Reduced(sel_subs_test,:);
        d_SVM_train = d_SVM_Reduced(sel_subs_train,:);
        % IQ_Labels = IQt>median(IQ(IQ(:,iqm)>0,1));
%         IQ_Labels = [HIQ HIQ HIQ];
%         IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2) 1]);
        IQ_Labels = ones([86 3 149])*10;
        for i=1:86
            IQ_Labels(i,:,:) = HIQ(i);
        end
        IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2)*size(IQ_Labels,3) 1]);
        
        SVMModel = fitcsvm(d_SVM_train, IQ_Labels(sel_subs_train),'KernelFunction','rbf');
%                         'OptimizeHyperparameters','auto');
        CVSVMModel = crossval(SVMModel);
        Acc = 0;
        for i=1:CVSVMModel.KFold
            predicted = predict(CVSVMModel.Trained{i},d_SVM_test);
            Acc = Acc + (1-sum(abs(predicted-IQ_Labels(sel_subs_test)))/size(sel_subs_train,2))*100;
        end
        Acc = Acc / CVSVMModel.KFold;
        Accs(k,j) = Acc;
        toc
    end
end
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(dim_vec, mean(Accs,1));
ylabel('Accuracy');
xlabel('#Dimensions after PCA');
box off
export_fig('SVM_PCA_Accs2.png','-r600');

save('SVMResult2.mat','SVMModel','Accs','sel_subs_test','d_SVM_train','CVSVMModel');
%% Neural Network
clc
[coeff,score,latent] = pca(d_SVM);
Accs = [];
for k=1:1
    for j=size(coeff,2)
        j
        d_NN_Reduced = d_SVM * coeff(:,1:j);
        
%         IQ_Labels = [HIQ HIQ HIQ];
%         IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2) 1]);
        IQ_Labels = ones([86 3 149])*10;
        for i=1:86
            IQ_Labels(i,:,:) = HIQ(i);
        end
        IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2)*size(IQ_Labels,3) 1]);
        
        ratio = 0.2;
        rng('shuffle') 
        t = 1:size(d_NN_Reduced,1);
        sel_subs_test = randperm(size(d_NN_Reduced,1),ceil(size(d_NN_Reduced,1)*ratio));
        sel_subs_train = t(setdiff(1:end,sel_subs_test));
        d_NN_test = d_NN_Reduced(sel_subs_test,:);
        d_NN_train = d_NN_Reduced(sel_subs_train,:);
        
        x = d_NN_train';
        t = IQ_Labels(sel_subs_train)';
        trainFcn = 'trainscg';

%         hiddenLayerSize = [j*(j+1)/2];
%         hiddenLayerSize = [251*125];
        hiddenLayerSize = [251*125 100];
        net = patternnet(hiddenLayerSize, trainFcn);

        net.divideParam.trainRatio = 64/100;
        net.divideParam.valRatio = 16/100;
        net.divideParam.testRatio = 0/100;
        
        net.trainParam.showWindow = false;

        [net,tr] = train(net,x,t);

        y = net(d_NN_test') > 0.5;
        e = gsubtract(IQ_Labels(sel_subs_test)',y);

        Accs(k,j) = sum(abs(e))/length(e)*100;
    end
end
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:size(coeff,2), mean(Accs(1,:),1));
ylabel('Accuracy');
xlabel('#Dimensions after PCA');
box off
export_fig('NN_PCA_Accs.png','-r600');

save('NNModel.mat');
%% KNN

[coeff,score,latent] = pca(d_SVM);
Accs = [];
dim_vec = linspace(1,size(coeff,2),1);
for k=1:1
    k
    for j = 1:length(dim_vec)
        j
        d_NN_Reduced = d_SVM * coeff(:,1:ceil(dim_vec(j)));
        
%         IQ_Labels = [HIQ HIQ HIQ];
%         IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2) 1]);
        IQ_Labels = ones([86 3 149])*10;
        for i=1:86
            IQ_Labels(i,:,:) = HIQ(i);
        end
        IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2)*size(IQ_Labels,3) 1]);

        ratio = 0.2;
        rng('shuffle') 
        t = 1:size(d_NN_Reduced,1);
        sel_subs_test = randperm(size(d_NN_Reduced,1),ceil(size(d_NN_Reduced,1)*ratio));
        sel_subs_train = t(setdiff(1:end,sel_subs_test));
        d_NN_test = d_NN_Reduced(sel_subs_test,:);
        d_NN_train = d_NN_Reduced(sel_subs_train,:);
        
%         KNNMdl = fitcknn(d_NN_train,IQ_Labels(sel_subs_train),'NumNeighbors',5,'Standardize',1);
        KNNMdl = fitcknn(d_NN_train,IQ_Labels(sel_subs_train),...
                        'OptimizeHyperparameters','auto');
%         SVMModel = fitcsvm(d_SVM_train, IQ_Labels(sel_subs_train), 'KernelFunction','rbf');
        CVSVMModel = crossval(KNNMdl);
        Acc = 0;
        for i=1:CVSVMModel.KFold
            predicted = predict(CVSVMModel.Trained{i},d_NN_test);
            Acc = Acc + (1-sum(abs(predicted-IQ_Labels(sel_subs_test)))/size(sel_subs_train,2))*100;
        end
        Acc = Acc / CVSVMModel.KFold;
        Accs(k,j) = Acc;
    end
end
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(dim_vec(1:7),mean(Accs(1,:),1));
ylabel('Accuracy');
xlabel('#Dimensions after PCA');
box off
export_fig('KNN_PCA_Accs3.png','-r600');

save('KNNModel3.mat');
%|   16 | Best   |     0.24677 |      178.91 |     0.24677 |     0.24682 |           61 |  correlation 
%% Bayesian
clc
[coeff,score,latent] = pca(d_SVM);
Accs = [];
dim_vec = linspace(1,size(coeff,2),1);
for k=1:1
    k
    for j=1:length(dim_vec)
        j
        d_NN_Reduced = d_SVM * coeff(:,1:dim_vec(j));
        
%         IQ_Labels = [HIQ HIQ HIQ];
%         IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2) 1]);
        
        IQ_Labels = ones([86 3 149])*10;
        for i=1:86
            IQ_Labels(i,:,:) = HIQ(i);
        end
        IQ_Labels = reshape(IQ_Labels, [size(IQ_Labels,1)*size(IQ_Labels,2)*size(IQ_Labels,3) 1]);
    
        ratio = 0.2;
        rng('shuffle') 
        t = 1:size(d_NN_Reduced,1);
        sel_subs_test = randperm(size(d_NN_Reduced,1),ceil(size(d_NN_Reduced,1)*ratio));
        sel_subs_train = t(setdiff(1:end,sel_subs_test));
        d_NN_test = d_NN_Reduced(sel_subs_test,:);
        d_NN_train = d_NN_Reduced(sel_subs_train,:);
        
%         KNNMdl = fitcknn(d_NN_train,IQ_Labels(sel_subs_train),'NumNeighbors',5,'Standardize',1);
%         KNNMdl = fitcknn(d_NN_train,IQ_Labels(sel_subs_train),...
%                         'OptimizeHyperparameters','auto');
%         SVMModel = fitcsvm(d_SVM_train, IQ_Labels(sel_subs_train), 'KernelFunction','rbf');
        BMdl = fitcnb(d_NN_train,IQ_Labels(sel_subs_train));
%                  'OptimizeHyperparameters','auto');
        CVSVMModel = crossval(BMdl);
        Acc = 0;
        for i=1:CVSVMModel.KFold
            predicted = predict(CVSVMModel.Trained{i},d_NN_test);
            Acc = Acc + (1-sum(abs(predicted-IQ_Labels(sel_subs_test)))/size(sel_subs_train,2))*100;
        end
        Acc = Acc / CVSVMModel.KFold;
        Accs(k,j) = Acc;
    end
end
close all

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(dim_vec,mean(Accs(1:end,:),1));
ylabel('Accuracy');
xlabel('#Dimensions after PCA');
box off
export_fig('Bayesian_PCA_Accs.png','-r600');

save('BayesianModel.mat');