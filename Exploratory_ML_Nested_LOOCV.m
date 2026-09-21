%% EXPLORATORY ML: leakage-free nested LOOCV on 150 broad-band features
% Input: 05_AllBands_Changes_dB_vs_PRE.xlsx
%
% Runs DURING_vs_PRE and POST_vs_PRE separately.
%
% Pipeline inside EACH LOOCV training fold:
%   1) min-max scaling learned from training data only
%   2) chi-square ranking -> top K1
%   3) mRMR -> top K2
%   4) classifier trained only on that fold
%
% Permutation testing repeats the ENTIRE nested feature-selection + LOOCV
% pipeline after permuting labels.
%
% Six classifiers are exploratory. Holm-adjusted permutation p-values are
% reported across the six models for each contrast.
%
% Requires Statistics and Machine Learning Toolbox.
% 1000 permutations can take a long time. For debugging only, temporarily
% reduce NPERM, then restore it to 1000 for final analysis.

clear; clc; close all;
rng(42);

CHANNELS = {'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4','P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4','TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};
BANDS = {'Delta','Theta','Alpha','Beta','Gamma'};
CONTRASTS = {'DURING_vs_PRE','POST_vs_PRE'};
MODELS = {'Naive Bayes','LDA','SVM','KNN','Decision Tree','Ensemble'};

K1 = 50;
K2 = 10;
NPERM = 1000;


filename=resolve_file('05_AllBands_Changes_dB_vs_PRE.xlsx');
Tall=readtable(filename,'VariableNamingRule','preserve');

expectedVars=cell(1,150);
k=1;
for ch=1:numel(CHANNELS)
    for b=1:numel(BANDS)
        expectedVars{k}=[CHANNELS{ch} '_' BANDS{b}];
        k=k+1;
    end
end

actualVars=Tall.Properties.VariableNames(6:end);
assert(isequal(actualVars,expectedVars), ...
    '150-feature/channel order does not match the verified extraction order.');
fprintf('150-feature/channel-order verification PASSED.\n');

for cc=1:numel(CONTRASTS)

    contrast = CONTRASTS{cc};

    T = Tall(strcmp(string(Tall.Contrast),contrast),:);

    % -------- Contrast-specific exclusions --------
    if strcmp(contrast,'DURING_vs_PRE')

        EXCLUDE_SUBJECTS = {'E_15','C_02', 'C_11', 'E_17'};

    elseif strcmp(contrast,'POST_vs_PRE')

        EXCLUDE_SUBJECTS = {'E_15'};

    else
        EXCLUDE_SUBJECTS = {};
    end

    % Apply exclusion
    T = apply_exclusions(T,EXCLUDE_SUBJECTS);

    fprintf('\nContrast: %s\n',contrast);

    if isempty(EXCLUDE_SUBJECTS)
        fprintf('Excluded subjects: none\n');
    else
        fprintf('Excluded subjects: %s\n', ...
            strjoin(EXCLUDE_SUBJECTS,', '));
    end


    Y=double(T.Label);
    X=table2array(T(:,expectedVars));

    assert(~any(isnan(X(:)) | isinf(X(:))), ...
        'Features contain NaN/Inf. Resolve missing data before ML.');
    assert(numel(unique(Y))==2 && all(ismember(unique(Y),[0 1])), ...
        'Labels must be 0=Control and 1=Experiment.');

    N=size(X,1);
    nModels=numel(MODELS);

    fprintf('\n============================================================\n');
    fprintf('ML CONTRAST: %s\n',contrast);
    fprintf('N=%d (Control=%d, Experiment=%d), predictors=%d\n', ...
        N,sum(Y==0),sum(Y==1),size(X,2));
    fprintf('K1=%d, K2=%d, permutations=%d\n',K1,K2,NPERM);
    fprintf('============================================================\n');

    % ---------- Actual nested LOOCV ----------
    folds=nested_feature_selection(X,Y,K1,K2);

    yhat=zeros(N,nModels);
    scores=zeros(N,nModels);

    for m=1:nModels
        for i=1:N
            [yhat(i,m),scores(i,m)] = train_predict_one( ...
                MODELS{m},folds.Xtr{i},folds.Ytr{i},folds.Xte{i});
        end
    end

    accuracy=nan(nModels,1);
    balancedAcc=nan(nModels,1);
    auc=nan(nModels,1);

    for m=1:nModels
        accuracy(m)=mean(yhat(:,m)==Y);
        sens=mean(yhat(Y==1,m)==1);
        spec=mean(yhat(Y==0,m)==0);
        balancedAcc(m)=(sens+spec)/2;
        [~,~,~,auc(m)]=perfcurve(Y,scores(:,m),1);
    end

    % ---------- Permutation null ----------
    aucNull=nan(NPERM,nModels);

    fprintf('Running permutations');
    for p=1:NPERM
        Yp=Y(randperm(N));

        % Rerun ALL preprocessing/selection inside each permuted LOOCV.
        foldsP=nested_feature_selection(X,Yp,K1,K2);

        scoreP=zeros(N,nModels);
        for m=1:nModels
            for i=1:N
                [~,scoreP(i,m)] = train_predict_one( ...
                    MODELS{m},foldsP.Xtr{i},foldsP.Ytr{i},foldsP.Xte{i});
            end
            [~,~,~,aucNull(p,m)]=perfcurve(Yp,scoreP(:,m),1);
        end

        if mod(p,max(1,round(NPERM/20)))==0
            fprintf('.');
        end
    end
    fprintf(' done.\n');

    pPerm=nan(nModels,1);
    for m=1:nModels
        pPerm(m)=(1+sum(aucNull(:,m)>=auc(m)))/(NPERM+1);
    end
    pHolm=holm_adjust(pPerm);

    R=table(string(MODELS(:)),100*accuracy,100*balancedAcc,auc,pPerm,pHolm, ...
        'VariableNames',{'Model','LOOCV_Accuracy_pct','BalancedAccuracy_pct', ...
        'AUC','Permutation_p','Permutation_p_Holm6'});
    R=sortrows(R,'AUC','descend');

    fprintf('\n=== Exploratory nested LOOCV: %s ===\n',contrast);
    disp(R);

    writetable(R,sprintf('ML_nested_LOOCV_%s_results.csv',contrast));

    % Feature-selection stability: how often each predictor was selected
    % across the N actual-data LOOCV folds.
    selectedCounts=zeros(150,1);
    for i=1:N
        idx=folds.SelectedIdx{i};
        selectedCounts(idx)=selectedCounts(idx)+1;
    end
    Stability=table(string(expectedVars(:)),selectedCounts,selectedCounts/N, ...
        'VariableNames',{'Feature','SelectedFolds','SelectionFrequency'});
    Stability=sortrows(Stability,'SelectionFrequency','descend');
    writetable(Stability,sprintf('ML_feature_stability_%s.csv',contrast));

    fprintf('Top 15 features by nested-fold selection frequency:\n');
    disp(Stability(1:15,:));
end


%% -------------------------- LOCAL FUNCTIONS -----------------------------
function folds=nested_feature_selection(X,Y,K1,K2)
    N=size(X,1);
    folds.Xtr=cell(N,1);
    folds.Xte=cell(N,1);
    folds.Ytr=cell(N,1);
    folds.SelectedIdx=cell(N,1);

    for i=1:N
        testMask=false(N,1);
        testMask(i)=true;
        trainMask=~testMask;

        Xtr0=X(trainMask,:);
        Xte0=X(testMask,:);
        Ytr=Y(trainMask);

        % Scaling learned ONLY from training subjects.
        mn=min(Xtr0,[],1);
        mx=max(Xtr0,[],1);
        rg=mx-mn;
        rg(rg==0)=1;

        Xtr=(Xtr0-mn)./rg;
        Xte=(Xte0-mn)./rg;

        % Avoid numerical excursions outside [0,1] in the held-out sample;
        % selection is still learned exclusively from Xtr/Ytr.
        Xte=max(0,min(1,Xte));

        k1u=min(K1,size(Xtr,2));
        idxChi=fscchi2(Xtr,Ytr);
        chiTop=idxChi(1:k1u);

        XtrChi=Xtr(:,chiTop);
        k2u=min(K2,size(XtrChi,2));
        idxMRMR=fscmrmr(XtrChi,Ytr);
        selected=chiTop(idxMRMR(1:k2u));

        folds.Xtr{i}=Xtr(:,selected);
        folds.Xte{i}=Xte(:,selected);
        folds.Ytr{i}=Ytr;
        folds.SelectedIdx{i}=selected;
    end
end

function [yh,s]=train_predict_one(modelName,Xtr,Ytr,Xte)
    switch modelName
        case 'Naive Bayes'
            try
                mdl=fitcnb(Xtr,Ytr,'ClassNames',[0 1]);
            catch
                mdl=fitcnb(Xtr,Ytr,'DistributionNames','kernel','ClassNames',[0 1]);
            end

        case 'LDA'
            % pseudoLinear is more stable than ordinary linear LDA for
            % small-N/high-dimensional training folds.
            mdl=fitcdiscr(Xtr,Ytr,'DiscrimType','pseudoLinear', ...
                'Prior','empirical','ClassNames',[0 1]);

        case 'SVM'
            % Data were already scaled inside the training fold.
            mdl=fitcsvm(Xtr,Ytr,'KernelFunction','linear', ...
                'Standardize',false,'ClassNames',[0 1]);

        case 'KNN'
            k=min(5,size(Xtr,1)-1);
            mdl=fitcknn(Xtr,Ytr,'NumNeighbors',max(k,1), ...
                'Standardize',false);

        case 'Decision Tree'
            mdl=fitctree(Xtr,Ytr,'ClassNames',[0 1]);

        case 'Ensemble'
            mdl=fitcensemble(Xtr,Ytr,'Method','Bag','ClassNames',[0 1]);

        otherwise
            error('Unknown model: %s',modelName);
    end

    [yh,sc]=predict(mdl,Xte);
    classes=mdl.ClassNames;
    posCol=find(classes==1,1);

    if isempty(posCol)
        error('Could not find positive class (1) in classifier.');
    end

    if size(sc,2)==1
        % Some binary learners may provide a single signed score.
        s=sc(1);
    else
        s=sc(1,posCol);
    end
end

function pAdj=holm_adjust(p)
    p=p(:);
    m=numel(p);
    [ps,ord]=sort(p);
    adjSorted=zeros(m,1);
    running=0;
    for i=1:m
        val=(m-i+1)*ps(i);
        running=max(running,val);
        adjSorted(i)=min(running,1);
    end
    pAdj=zeros(m,1);
    pAdj(ord)=adjSorted;
end

function filename=resolve_file(baseName)
    if isfile(baseName), filename=baseName; return; end
    [~,stem,ext]=fileparts(baseName);
    d=dir([stem '*' ext]);
    if isempty(d), error('Could not find %s.',baseName); end
    [~,ix]=max([d.datenum]); filename=d(ix).name;
    fprintf('Using file: %s\n',filename);
end

function T=apply_exclusions(T,subjects)
    if isempty(subjects), return; end
    T=T(~ismember(string(T.SubjectUID),string(subjects)),:);
end
