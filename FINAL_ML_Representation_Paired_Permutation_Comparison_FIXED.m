%% FINAL ML REPRESENTATION PAIRED PERMUTATION COMPARISON - ROBUST VERSION
% =========================================================================
% Compares the FINAL permutation-fold-fixed repeated stratified 5-fold CV
% results across:
%   Level 1 = 1 whole-scalp 34-43 Hz feature
%   Level 2 = 30 channel-wise 34-43 Hz features
%   Level 3 = 150 multiband features
%
% IMPORTANT:
% This script does NOT need obsFoldID to be stored in the *_details.mat files.
% It only requires:
%   - Rep
%   - nullMeanAUC
%   - meanAUC
%   - SubjectOrder
%
% If pPerm was not saved, it is recalculated from nullMeanAUC.
%
% Pairwise representation comparisons use the SAME permutation index across
% levels. This is valid only for the final scripts that used the same:
%   PERM_SEED = 91000
%   PERM_FOLD_SEED = 150000
%   NREPEATS = 50
%   K = 5
% and the same sorted SubjectUID order.
%
% Outputs:
%   WELCH_FINAL_ML_Representation_Summary.csv
%   WELCH_FINAL_ML_Pairwise_AUC_Comparisons.csv
%   WELCH_FINAL_ML_AllResults.mat
% =========================================================================

clear; clc; close all;

%% ----------------------------- SETTINGS --------------------------------
METHOD_TAG = 'WELCH';
CONTRASTS = {'DURING_vs_PRE','POST_vs_PRE'};

LEVEL_NAMES = { ...
    'Level1_WholeScalp34_43', ...
    'Level2_30ch34_43', ...
    'Level3_150BroadBand'};

LEVEL_LABELS = { ...
    'Level 1: whole-scalp 34-43 Hz', ...
    'Level 2: 30-channel 34-43 Hz', ...
    'Level 3: 150-feature multiband'};

PAIR_INDEX = [1 2; 1 3; 2 3];

All = struct();

%% ============================ MAIN LOOP =================================
for c = 1:numel(CONTRASTS)

    contrast = CONTRASTS{c};

    fprintf('\n============================================================\n');
    fprintf('%s | %s\n',METHOD_TAG,contrast);
    fprintf('============================================================\n');

    D = cell(1,3);
    fileNames = strings(1,3);

    %% ---------------------- LOAD LEVEL FILES -----------------------------
    for L = 1:3

        pattern = sprintf('%s_RS5CV_PermFoldFixed_Level%d_%s_details.mat', ...
            METHOD_TAG,L,contrast);

        fname = resolve_exact_or_latest(pattern);
        fileNames(L) = string(fname);

        fprintf('Loading Level %d: %s\n',L,fname);

        D{L} = load(fname);

        % Only require what is actually necessary for the final comparison.
        requiredFields = {'Rep','nullMeanAUC','meanAUC','SubjectOrder'};

        for rr = 1:numel(requiredFields)
            assert(isfield(D{L},requiredFields{rr}), ...
                'Missing field "%s" in %s.',requiredFields{rr},fname);
        end

        % Normalize shapes
        D{L}.nullMeanAUC = double(D{L}.nullMeanAUC(:));
        D{L}.meanAUC = double(D{L}.meanAUC);

        % Recalculate representation-level permutation p if absent.
        if ~isfield(D{L},'pPerm') || isempty(D{L}.pPerm)
            nPermTmp = numel(D{L}.nullMeanAUC);
            D{L}.pPerm = ...
                (1 + sum(D{L}.nullMeanAUC >= D{L}.meanAUC)) / (nPermTmp + 1);
        end
    end

    %% ------------------------ VALIDATION ---------------------------------
    subjRef = string(D{1}.SubjectOrder.SubjectUID);
    yRef    = double(D{1}.SubjectOrder.Label(:));

    for L = 2:3

        assert(isequal(subjRef,string(D{L}.SubjectOrder.SubjectUID)), ...
            ['Subject order differs between Level 1 and Level %d for %s.\n' ...
             'Do not perform the paired comparison until the same subjects ' ...
             'and sorted SubjectUID order are used.'],L,contrast);

        assert(isequal(yRef,double(D{L}.SubjectOrder.Label(:))), ...
            'Labels differ between Level 1 and Level %d for %s.',L,contrast);
    end

    nPerm = numel(D{1}.nullMeanAUC);

    for L = 2:3
        assert(numel(D{L}.nullMeanAUC)==nPerm, ...
            'Number of permutations differs across levels for %s.',contrast);
    end

    % Optional fold validation if any version stored folds.
    [hasFold1,fold1] = get_saved_fold_matrix(D{1});

    for L = 2:3
        [hasFoldL,foldL] = get_saved_fold_matrix(D{L});

        if hasFold1 && hasFoldL
            assert(isequal(fold1,foldL), ...
                'Observed fold assignments differ between Levels 1 and %d for %s.', ...
                L,contrast);
        end
    end

    fprintf(['Validation PASSED: identical subjects, labels, and Nperm=%d.\n' ...
             'Observed-fold validation: %s\n'], ...
             nPerm,fold_validation_text(D));

    %% ---------------- REPRESENTATION-LEVEL RESULTS -----------------------
    obsAUC = zeros(3,1);
    sdAUC  = zeros(3,1);
    rawP   = zeros(3,1);
    nullMu = zeros(3,1);
    nullSD = zeros(3,1);

    for L = 1:3
        obsAUC(L) = D{L}.meanAUC;
        sdAUC(L)  = std(double(D{L}.Rep.AUC));
        rawP(L)   = double(D{L}.pPerm);
        nullMu(L) = mean(D{L}.nullMeanAUC);
        nullSD(L) = std(D{L}.nullMeanAUC);
    end

    holmP = holm_adjust(rawP);

    RepresentationTable = table( ...
        repmat(string(METHOD_TAG),3,1), ...
        repmat(string(contrast),3,1), ...
        (1:3)', ...
        string(LEVEL_NAMES(:)), ...
        string(LEVEL_LABELS(:)), ...
        obsAUC,sdAUC,rawP,holmP,nullMu,nullSD, ...
        'VariableNames',{ ...
        'Method','Contrast','Level','FeatureSpace','Representation', ...
        'Mean_CV_AUC','SD_CV_AUC','Permutation_p','Holm_Adjusted_p', ...
        'PermutationNullMeanAUC','PermutationNullSD_AUC'});

    fprintf('\nRepresentation-level results:\n');
    disp(RepresentationTable);

    %% ---------------- DIRECT PAIRED AUC COMPARISONS ---------------------
    nPairs = size(PAIR_INDEX,1);

    Pair = strings(nPairs,1);
    A = strings(nPairs,1);
    B = strings(nPairs,1);

    AUC_A = nan(nPairs,1);
    AUC_B = nan(nPairs,1);
    ObsDifference_AminusB = nan(nPairs,1);

    TwoSidedPermutation_p = nan(nPairs,1);
    OneSided_AgreaterB_p  = nan(nPairs,1);

    NullDifferenceMean = nan(nPairs,1);
    NullDifferenceSD = nan(nPairs,1);
    NullDifferenceCI95_Lower = nan(nPairs,1);
    NullDifferenceCI95_Upper = nan(nPairs,1);

    for pp = 1:nPairs

        a = PAIR_INDEX(pp,1);
        b = PAIR_INDEX(pp,2);

        Pair(pp) = sprintf('L%d_vs_L%d',a,b);
        A(pp) = string(LEVEL_LABELS{a});
        B(pp) = string(LEVEL_LABELS{b});

        AUC_A(pp) = obsAUC(a);
        AUC_B(pp) = obsAUC(b);

        obsDiff = obsAUC(a) - obsAUC(b);
        ObsDifference_AminusB(pp) = obsDiff;

        % Pair null results by permutation index.
        % The final Level 1/2/3 scripts used identical deterministic
        % label-permutation and permutation-fold seeds.
        nullDiff = D{a}.nullMeanAUC - D{b}.nullMeanAUC;

        % PRIMARY: two-sided direct comparison
        TwoSidedPermutation_p(pp) = ...
            (1 + sum(abs(nullDiff) >= abs(obsDiff))) / (nPerm + 1);

        % Supplementary directional A > B comparison
        OneSided_AgreaterB_p(pp) = ...
            (1 + sum(nullDiff >= obsDiff)) / (nPerm + 1);

        NullDifferenceMean(pp) = mean(nullDiff);
        NullDifferenceSD(pp)   = std(nullDiff);

        q = prctile(nullDiff,[2.5 97.5]);
        NullDifferenceCI95_Lower(pp) = q(1);
        NullDifferenceCI95_Upper(pp) = q(2);
    end

    PairwiseHolm_p = holm_adjust(TwoSidedPermutation_p);

    PairwiseTable = table( ...
        repmat(string(METHOD_TAG),nPairs,1), ...
        repmat(string(contrast),nPairs,1), ...
        Pair,A,B,AUC_A,AUC_B,ObsDifference_AminusB, ...
        TwoSidedPermutation_p,PairwiseHolm_p,OneSided_AgreaterB_p, ...
        NullDifferenceMean,NullDifferenceSD, ...
        NullDifferenceCI95_Lower,NullDifferenceCI95_Upper, ...
        'VariableNames',{ ...
        'Method','Contrast','Comparison','Representation_A','Representation_B', ...
        'MeanAUC_A','MeanAUC_B','Observed_AUC_Difference_AminusB', ...
        'PairedPermutation_p_TwoSided','Holm_Adjusted_p_TwoSided', ...
        'PairedPermutation_p_OneSided_AgreaterB', ...
        'NullDifferenceMean','NullDifferenceSD', ...
        'NullDifference95_Lower','NullDifference95_Upper'});

    fprintf('\nDirect paired AUC-difference tests:\n');
    disp(PairwiseTable);

    %% -------------------------- STORE -----------------------------------
    All.(contrast).RepresentationTable = RepresentationTable;
    All.(contrast).PairwiseTable = PairwiseTable;
    All.(contrast).FileNames = fileNames;
end

%% ----------------------- COMBINE AND SAVE -------------------------------
RepresentationSummary = [
    All.DURING_vs_PRE.RepresentationTable
    All.POST_vs_PRE.RepresentationTable];

PairwiseSummary = [
    All.DURING_vs_PRE.PairwiseTable
    All.POST_vs_PRE.PairwiseTable];

writetable(RepresentationSummary, ...
    sprintf('%s_FINAL_ML_Representation_Summary.csv',METHOD_TAG));

writetable(PairwiseSummary, ...
    sprintf('%s_FINAL_ML_Pairwise_AUC_Comparisons.csv',METHOD_TAG));

save(sprintf('%s_FINAL_ML_AllResults.mat',METHOD_TAG), ...
    'RepresentationSummary','PairwiseSummary','All');

%% ------------------------- CLEAN CONSOLE SUMMARY ------------------------
fprintf('\n\n============================================================\n');
fprintf('FINAL ML REPRESENTATION SUMMARY\n');
fprintf('============================================================\n');

for c = 1:numel(CONTRASTS)

    contrast = CONTRASTS{c};
    R = All.(contrast).RepresentationTable;
    P = All.(contrast).PairwiseTable;

    fprintf('\n%s\n',contrast);

    for L = 1:height(R)
        fprintf('  L%d: AUC=%.3f, perm p=%.4f, Holm p=%.4f\n', ...
            R.Level(L),R.Mean_CV_AUC(L), ...
            R.Permutation_p(L),R.Holm_Adjusted_p(L));
    end

    fprintf('  Direct paired AUC comparisons:\n');

    for pp = 1:height(P)
        fprintf('    %s: delta AUC=%+.3f, two-sided p=%.4f, Holm p=%.4f\n', ...
            P.Comparison(pp), ...
            P.Observed_AUC_Difference_AminusB(pp), ...
            P.PairedPermutation_p_TwoSided(pp), ...
            P.Holm_Adjusted_p_TwoSided(pp));
    end
end

fprintf('\nSaved:\n');
fprintf('  %s_FINAL_ML_Representation_Summary.csv\n',METHOD_TAG);
fprintf('  %s_FINAL_ML_Pairwise_AUC_Comparisons.csv\n',METHOD_TAG);
fprintf('  %s_FINAL_ML_AllResults.mat\n',METHOD_TAG);
fprintf('============================================================\n');


%% =========================== LOCAL FUNCTIONS ============================

function filename = resolve_exact_or_latest(baseName)

if isfile(baseName)
    filename = baseName;
    return;
end

[~,stem,ext] = fileparts(baseName);
d = dir([stem '*' ext]);

if isempty(d)
    error(['Could not find required file:\n  %s\n' ...
           'Run this script in the folder containing the final Level 1-3 *_details.mat files.'], ...
           baseName);
end

[~,ix] = max([d.datenum]);
filename = d(ix).name;
end


function pAdj = holm_adjust(p)
% Holm-Bonferroni adjusted p-values.

p = double(p(:));
m = numel(p);

[pSort,idx] = sort(p,'ascend');
adjSort = nan(m,1);

for i = 1:m
    adjSort(i) = (m-i+1)*pSort(i);
end

for i = 2:m
    adjSort(i) = max(adjSort(i),adjSort(i-1));
end

adjSort = min(adjSort,1);

pAdj = nan(m,1);
pAdj(idx) = adjSort;
end


function [hasFold,foldMatrix] = get_saved_fold_matrix(S)
% Accept either field name used across earlier script versions.

hasFold = false;
foldMatrix = [];

if isfield(S,'obsFoldID')
    hasFold = true;
    foldMatrix = S.obsFoldID;
elseif isfield(S,'foldID')
    hasFold = true;
    foldMatrix = S.foldID;
end
end


function txt = fold_validation_text(D)

allHave = true;

for L = 1:3
    [hasFold,~] = get_saved_fold_matrix(D{L});
    allHave = allHave && hasFold;
end

if allHave
    txt = 'PASSED (saved fold matrices were present and matched)';
else
    txt = ['not required; one or more files did not save fold matrices. ' ...
           'Subject/label/permutation-length checks passed.'];
end
end
