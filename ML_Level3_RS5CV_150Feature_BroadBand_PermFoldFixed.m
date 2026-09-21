%% LEVEL 3 - REPEATED STRATIFIED 5-FOLD CV
% 150 channel x broad-band features. Standalone script.
% Separate from your existing chi-square -> mRMR -> six-classifier analysis.
% PERMUTATION FIX: each permutation receives new 5-fold partitions stratified on Yp.
clear; clc; close all;

METHOD_TAG='WELCH';
CONTRASTS={'DURING_vs_PRE','POST_vs_PRE'};
EXCLUDE_DURING={'E_15','C_02','E_17'};
EXCLUDE_POST={'E_15'};
K=5; NREPEATS=50; NPERM=1000;
BASE_SEED=42000; PERM_SEED=91000; PERM_FOLD_SEED=150000;

CHANNELS={'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4','P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4','TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};
BANDS={'Delta','Theta','Alpha','Beta','Gamma'};
filename=resolve_file('05_AllBands_Changes_dB_vs_PRE.xlsx'); Tall=readtable(filename,'VariableNamingRule','preserve');
expectedVars=cell(1,150); k=1; for ch=1:numel(CHANNELS), for b=1:numel(BANDS), expectedVars{k}=[CHANNELS{ch} '_' BANDS{b}]; k=k+1; end, end
assert(isequal(Tall.Properties.VariableNames(6:end),expectedVars),'150-feature order check failed.');
fprintf('\n%s - RS5CV ML LEVEL 3 - PERMUTATION-FOLD FIXED\n',METHOD_TAG);
for cc=1:numel(CONTRASTS)
    contrast=CONTRASTS{cc}; T=Tall(strcmp(string(Tall.Contrast),contrast),:);
    if strcmp(contrast,'DURING_vs_PRE'), exclusions=EXCLUDE_DURING; else, exclusions=EXCLUDE_POST; end
    T=apply_exclusions(T,exclusions); T=sortrows(T,'SubjectUID'); Y=double(T.Label); X=table2array(T(:,expectedVars)); check_binary_data(X,Y);
    foldID=make_repeated_stratified_folds(Y,K,NREPEATS,BASE_SEED); Rep=repeated_cv_balanced_ridge(X,Y,foldID,K);
    meanAUC=mean(Rep.AUC); sdAUC=std(Rep.AUC); meanAcc=mean(Rep.Accuracy); sdAcc=std(Rep.Accuracy); meanBal=mean(Rep.BalancedAccuracy); sdBal=std(Rep.BalancedAccuracy); meanSens=mean(Rep.Sensitivity); meanSpec=mean(Rep.Specificity);
    N=numel(Y); nullMeanAUC=nan(NPERM,1);
    fprintf('\n%s | LEVEL 3 | %s\nN=%d (C=%d,E=%d), predictors=150\nExcluded: %s\nObserved mean CV AUC=%.4f\nRunning %d permutations',METHOD_TAG,contrast,N,sum(Y==0),sum(Y==1),show_exclusions(exclusions),meanAUC,NPERM);
    for p=1:NPERM
        % Permute labels deterministically.
        oldRng=rng; rng(PERM_SEED+p-1,'twister'); Yp=Y(randperm(N)); rng(oldRng);

        % IMPORTANT: regenerate repeated stratified folds from the PERMUTED labels.
        % This mirrors the observed analysis and avoids reusing folds stratified
        % according to the original labels.
        permFoldID=make_repeated_stratified_folds(Yp,K,NREPEATS, ...
            PERM_FOLD_SEED+(p-1)*NREPEATS);

        Mp=repeated_cv_balanced_ridge(X,Yp,permFoldID,K);
        nullMeanAUC(p)=mean(Mp.AUC);
        if mod(p,max(1,round(NPERM/20)))==0, fprintf('.'); end
    end
    fprintf(' done.\n'); pPerm=(1+sum(nullMeanAUC>=meanAUC))/(NPERM+1); nullCenter=mean(nullMeanAUC); nullSD=std(nullMeanAUC);
    if abs(nullCenter-0.5)>0.03, warning('Permutation null mean AUC=%.3f.',nullCenter); end
    R=table(string(METHOD_TAG),string(contrast),"Level3_150BroadBand",N,sum(Y==0),sum(Y==1),150,K,NREPEATS,NPERM,meanAUC,sdAUC,100*meanAcc,100*sdAcc,100*meanBal,100*sdBal,100*meanSens,100*meanSpec,pPerm,nullCenter,nullSD, ...
      'VariableNames',{'Method','Contrast','FeatureSpace','N','N_Control','N_Experiment','N_Predictors','K_Folds','N_Repeats','N_Permutations','Mean_CV_AUC','SD_CV_AUC','Mean_CV_Accuracy_pct','SD_CV_Accuracy_pct','Mean_CV_BalancedAccuracy_pct','SD_CV_BalancedAccuracy_pct','Mean_CV_Sensitivity_pct','Mean_CV_Specificity_pct','Permutation_p','PermutationNullMeanAUC','PermutationNullSD_AUC'});
    disp(R); out=sprintf('%s_RS5CV_PermFoldFixed_Level3_%s',METHOD_TAG,contrast); writetable(R,[out '_summary.csv']); writetable(Rep,[out '_repeat_metrics.csv']); SubjectOrder=table(string(T.SubjectUID),Y,'VariableNames',{'SubjectUID','Label'}); writetable(SubjectOrder,[out '_subject_order.csv']); save([out '_details.mat'],'foldID','Rep','nullMeanAUC','meanAUC','pPerm','SubjectOrder');
end

%% =========================== LOCAL FUNCTIONS ============================
function foldID = make_repeated_stratified_folds(Y,K,NREPEATS,BASE_SEED)
Y=double(Y(:)); N=numel(Y); foldID=zeros(N,NREPEATS); oldRng=rng;
for r=1:NREPEATS
    rng(BASE_SEED+r-1,'twister');
    cvp=cvpartition(Y,'KFold',K);
    for f=1:K
        foldID(test(cvp,f),r)=f;
    end
    if any(foldID(:,r)==0), error('Fold assignment failed in repeat %d.',r); end
end
rng(oldRng);
end

function M = repeated_cv_balanced_ridge(X,Y,foldID,K)
X=double(X); Y=double(Y(:)); N=size(X,1); R=size(foldID,2);
AUC=nan(R,1); Accuracy=nan(R,1); BalancedAccuracy=nan(R,1); Sensitivity=nan(R,1); Specificity=nan(R,1);
for r=1:R
    yhat=nan(N,1); scores=nan(N,1);
    for f=1:K
        te=(foldID(:,r)==f); tr=~te;
        Xtr0=X(tr,:); Xte0=X(te,:); Ytr=Y(tr);
        if numel(unique(Ytr))<2, error('Repeat %d fold %d has one training class.',r,f); end
        mu=mean(Xtr0,1); sd=std(Xtr0,0,1); sd(sd==0 | ~isfinite(sd))=1;
        Xtr=(Xtr0-mu)./sd; Xte=(Xte0-mu)./sd;
        n0=sum(Ytr==0); n1=sum(Ytr==1);
        W=zeros(size(Ytr)); W(Ytr==0)=numel(Ytr)/(2*n0); W(Ytr==1)=numel(Ytr)/(2*n1);
        lambda=1/size(Xtr,1);
        mdl=fitclinear(Xtr,Ytr,'Learner','logistic','Regularization','ridge', ...
            'Lambda',lambda,'Weights',W,'ClassNames',[0 1]);
        [yh,sc]=predict(mdl,Xte);
        cls=double(mdl.ClassNames(:)); posCol=find(cls==1,1);
        yhat(te)=double(yh); scores(te)=sc(:,posCol);
    end
    Accuracy(r)=mean(yhat==Y);
    sens=mean(yhat(Y==1)==1); spec=mean(yhat(Y==0)==0);
    Sensitivity(r)=sens; Specificity(r)=spec; BalancedAccuracy(r)=(sens+spec)/2;
    [~,~,~,AUC(r)]=perfcurve(Y,scores,1);
end
M=table((1:R)',AUC,Accuracy,BalancedAccuracy,Sensitivity,Specificity, ...
    'VariableNames',{'Repeat','AUC','Accuracy','BalancedAccuracy','Sensitivity','Specificity'});
end

function T=apply_exclusions(T,subjects)
if isempty(subjects), return; end
T=T(~ismember(string(T.SubjectUID),string(subjects)),:);
end

function filename=resolve_file(baseName)
if isfile(baseName), filename=baseName; return; end
[~,stem,ext]=fileparts(baseName); d=dir([stem '*' ext]);
if isempty(d), error('Could not find %s in current folder.',baseName); end
[~,ix]=max([d.datenum]); filename=d(ix).name; fprintf('Using file: %s\n',filename);
end

function txt=show_exclusions(x)
if isempty(x), txt='none'; else, txt=strjoin(x,', '); end
end

function check_binary_data(X,Y)
assert(~any(isnan(X(:)) | isinf(X(:))),'Features contain NaN/Inf.');
assert(numel(unique(Y))==2 && all(ismember(unique(Y),[0 1])),'Labels must be 0/1.');
end
