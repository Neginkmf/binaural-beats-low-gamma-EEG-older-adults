%% PRIMARY: Whole-scalp 34-43 Hz DURING vs PRE
% Hypothesis-driven primary EEG analysis.
% Input: 06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx
%
% Primary test:
%   Experiment vs Control on WholeScalp_Avg_dB for DURING_vs_PRE
%   using Welch's independent-samples t-test.
%
% IMPORTANT:
% - Includes all participants by default.
% - Only add an exclusion below if it was independently defined by QC,
%   not because it changes statistical significance.
% - The script verifies the exact 30-channel order in the input file and
%   verifies that WholeScalp_Avg_dB is the mean of the 30 channel-wise dB
%   changes.

clear; clc; close all;

EXPECTED_CHANNELS = {'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4','P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4','TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};
EXCLUDE_SUBJECTS = {'E_15', 'C_02', 'E_17'};   % e.g. {'E_15'} ONLY if pre-specified QC exclusion

filename = resolve_file('06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx');
T = readtable(filename,'VariableNamingRule','preserve');

verify_narrow_file(T, EXPECTED_CHANNELS);

T = T(strcmp(string(T.Contrast),'DURING_vs_PRE'),:);
T = apply_exclusions(T, EXCLUDE_SUBJECTS);

labels = double(T.Label);
y = T.WholeScalp_Avg_dB;

assert(all(ismember(unique(labels),[0 1])) && numel(unique(labels))==2, ...
    'Label must contain both 0=Control and 1=Experiment.');

C = y(labels==0);
E = y(labels==1);

% Welch t-test: call as E,C so CI is Experiment-Control.
[~,pWelch,ci,stats] = ttest2(E,C,'Vartype','unequal');

% Standardized effect size.
[g,d] = hedges_g(E,C);

% Rank-based sensitivity analysis.
pRank = ranksum(E,C);

% Descriptives.
nC = numel(C); nE = numel(E);
meanC = mean(C); meanE = mean(E);
sdC = std(C); sdE = std(E);
medC = median(C); medE = median(E);
difference = meanE - meanC;

fprintf('\n============================================================\n');
fprintf('PRIMARY EEG: WHOLE-SCALP 34-43 Hz, DURING - PRE\n');
fprintf('============================================================\n');
fprintf('Control:    n=%d, mean=%+.3f dB, SD=%.3f, median=%+.3f\n', ...
    nC,meanC,sdC,medC);
fprintf('Experiment: n=%d, mean=%+.3f dB, SD=%.3f, median=%+.3f\n', ...
    nE,meanE,sdE,medE);
fprintf('Experiment-Control difference = %+.3f dB\n',difference);
fprintf('Welch t(%0.2f) = %+.3f, p = %.6f\n',stats.df,stats.tstat,pWelch);
fprintf('95%% CI for Experiment-Control = [%+.3f, %+.3f] dB\n',ci(1),ci(2));
fprintf('Cohen d = %+.3f; Hedges g = %+.3f\n',d,g);
fprintf('Mann-Whitney/ranksum sensitivity p = %.6f\n',pRank);
fprintf('============================================================\n\n');

Results = table("DURING_vs_PRE",nC,nE,meanC,sdC,medC,meanE,sdE,medE, ...
    difference,stats.tstat,stats.df,pWelch,ci(1),ci(2),d,g,pRank, ...
    'VariableNames',{'Contrast','N_Control','N_Experiment', ...
    'Mean_Control_dB','SD_Control_dB','Median_Control_dB', ...
    'Mean_Experiment_dB','SD_Experiment_dB','Median_Experiment_dB', ...
    'Difference_EminusC_dB','Welch_t','Welch_df','Welch_p', ...
    'CI95_Lower_dB','CI95_Upper_dB','Cohens_d','Hedges_g','Ranksum_p'});
writetable(Results,'Primary_34_43Hz_DURING_vs_PRE_results.csv');

% Individual-subject figure.
figure('Name','Primary 34-43 Hz DURING-PRE');
hold on;
rng(42);
scatter(1 + 0.08*(rand(size(C))-0.5), C, 42, 'filled');
scatter(2 + 0.08*(rand(size(E))-0.5), E, 42, 'filled');
errorbar(1,meanC,sdC/sqrt(nC),'k','LineWidth',1.5,'CapSize',10);
errorbar(2,meanE,sdE/sqrt(nE),'k','LineWidth',1.5,'CapSize',10);
yline(0,'--');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Control','Experiment'});
ylabel('Whole-scalp 34-43 Hz change (dB)');
title(sprintf('DURING-PRE: E-C = %+.2f dB, Welch p = %.3f',difference,pWelch));
box on; hold off;
saveas(gcf,'Primary_34_43Hz_DURING_vs_PRE.png');


%% -------------------------- LOCAL FUNCTIONS -----------------------------
function filename = resolve_file(baseName)
    if isfile(baseName)
        filename = baseName;
        return;
    end
    [~,stem,ext] = fileparts(baseName);
    d = dir([stem '*' ext]);
    if isempty(d)
        error('Could not find %s in the current folder.',baseName);
    end
    [~,ix] = max([d.datenum]);
    filename = d(ix).name;
    fprintf('Using file: %s\n',filename);
end

function verify_narrow_file(T,expectedChannels)
    meta = {'SubjectUID','SubjectID','Group','Label','Contrast'};
    for i=1:numel(meta)
        assert(ismember(meta{i},T.Properties.VariableNames), ...
            'Missing required column: %s',meta{i});
    end

    expectedVars = strcat(expectedChannels,'_Gamma34_43');
    actualVars = T.Properties.VariableNames(6:35);
    assert(isequal(actualVars,expectedVars), ...
        ['CHANNEL ORDER CHECK FAILED. Expected: ' strjoin(expectedChannels,', ')]);
    assert(strcmp(T.Properties.VariableNames{36},'WholeScalp_Avg_dB'), ...
        'Column 36 must be WholeScalp_Avg_dB.');

    X = table2array(T(:,expectedVars));
    err = max(abs(mean(X,2)-T.WholeScalp_Avg_dB));
    assert(err < 1e-8, ...
        'WholeScalp_Avg_dB does not equal the mean of the 30 channel-wise dB changes.');

    fprintf('Channel-order verification PASSED (30 channels).\n');
    fprintf('Whole-scalp mean verification PASSED (max error %.3g).\n',err);
end

function T = apply_exclusions(T,subjects)
    if isempty(subjects), return; end
    rm = ismember(string(T.SubjectUID),string(subjects));
    if any(rm)
        fprintf('QC exclusions: %s\n',strjoin(cellstr(string(T.SubjectUID(rm))),', '));
    end
    T = T(~rm,:);
end

function [g,d] = hedges_g(E,C)
    nE=numel(E); nC=numel(C);
    sp = sqrt(((nE-1)*var(E)+(nC-1)*var(C))/(nE+nC-2));
    if sp==0
        d=NaN; g=NaN; return;
    end
    d=(mean(E)-mean(C))/sp;
    J=1-3/(4*(nE+nC)-9);
    g=J*d;
end
