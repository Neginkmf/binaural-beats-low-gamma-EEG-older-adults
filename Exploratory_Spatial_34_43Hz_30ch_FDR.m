%% EXPLORATORY SPATIAL: 30 channels, 34-43 Hz + BH-FDR
% Input: 06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx
%
% Runs BOTH contrasts separately:
%   1) DURING_vs_PRE
%   2) POST_vs_PRE
%
% For each of the 30 channels:
%   - Welch t-test (Experiment vs Control)
%   - Hedges g
%   - Benjamini-Hochberg FDR across the 30 channels WITHIN that contrast
%
% The script explicitly verifies the channel order against the current
% extraction file before doing any statistics.

clear; clc; close all;

CHANNELS = {'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4','P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4','TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};
CONTRASTS = {'DURING_vs_PRE','POST_vs_PRE'};


filename=resolve_file('06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx');
Tall=readtable(filename,'VariableNamingRule','preserve');

gammaVars=strcat(CHANNELS,'_Gamma34_43');
assert(isequal(Tall.Properties.VariableNames(6:35),gammaVars), ...
    ['CHANNEL ORDER CHECK FAILED. File does not match: ' strjoin(CHANNELS,', ')]);
assert(strcmp(Tall.Properties.VariableNames{36},'WholeScalp_Avg_dB'), ...
    'Expected WholeScalp_Avg_dB as the last column.');
fprintf('Channel-order verification PASSED:\n%s\n\n',strjoin(CHANNELS,' -> '));

for cc = 1:numel(CONTRASTS)

    contrastName = CONTRASTS{cc};

    % Select contrast
    T = Tall(strcmp(string(Tall.Contrast),contrastName),:);

    % Contrast-specific exclusions
    if strcmp(contrastName,'DURING_vs_PRE')

        % E15 already absent from FFT extraction
        % C02 = DURING-specific exclusion
        EXCLUDE_SUBJECTS = {'E_15','C_02', 'E_17', 'C_11'};

    elseif strcmp(contrastName,'POST_vs_PRE')

        % C02 stays in POST
        EXCLUDE_SUBJECTS = {'E_15'};

    else
        EXCLUDE_SUBJECTS = {};
    end

    T = apply_exclusions(T,EXCLUDE_SUBJECTS);

    fprintf('\nContrast: %s\n',contrastName);
    fprintf('Excluded subjects: %s\n', ...
        strjoin(EXCLUDE_SUBJECTS,', '));

    labels=double(T.Label);
    X=table2array(T(:,gammaVars));

    nCh=numel(CHANNELS);
    tstat=nan(nCh,1);
    df=nan(nCh,1);
    pRaw=nan(nCh,1);
    ciLow=nan(nCh,1);
    ciHigh=nan(nCh,1);
    g=nan(nCh,1);
    meanC=nan(nCh,1);
    meanE=nan(nCh,1);

    for ch=1:nCh
        C=X(labels==0,ch);
        E=X(labels==1,ch);

        [~,p,ci,st]=ttest2(E,C,'Vartype','unequal');
        pRaw(ch)=p;
        tstat(ch)=st.tstat;
        df(ch)=st.df;
        ciLow(ch)=ci(1);
        ciHigh(ch)=ci(2);
        meanC(ch)=mean(C);
        meanE(ch)=mean(E);
        g(ch)=hedges_g(E,C);
    end

    pFDR=bh_fdr(pRaw);
    diffMean=meanE-meanC;

    R=table(string(CHANNELS(:)),meanC,meanE,diffMean,tstat,df,pRaw,pFDR, ...
        ciLow,ciHigh,g, ...
        'VariableNames',{'Channel','Mean_Control_dB','Mean_Experiment_dB', ...
        'Difference_EminusC_dB','Welch_t','Welch_df','p_raw','p_FDR_BH', ...
        'CI95_Lower_dB','CI95_Upper_dB','Hedges_g'});
    R=sortrows(R,'p_raw');

    fprintf('\n=== 34-43 Hz spatial analysis: %s ===\n',contrastName);
    disp(R(1:min(10,height(R)),:));
    fprintf('Raw p < .05: %d / 30\n',sum(pRaw<0.05));
    fprintf('BH-FDR q < .05: %d / 30\n',sum(pFDR<0.05));

    outName=sprintf('Spatial_34_43Hz_%s_30channels_FDR.csv',contrastName);
    writetable(R,outName);

    % Effect-size-by-channel figure in the VERIFIED file order.
    figure('Name',['34-43 Hz spatial effects - ' contrastName]);
    bar(g);
    yline(0,'--');
    xticks(1:nCh);
    xticklabels(CHANNELS);
    xtickangle(45);
    ylabel('Hedges g (Experiment - Control)');
    title(sprintf('34-43 Hz spatial effect sizes: %s',strrep(contrastName,'_',' ')));
    box on;
    saveas(gcf,sprintf('Spatial_34_43Hz_%s_HedgesG.png',contrastName));

    % -log10 p plot: raw threshold and FDR-significant channels.
    figure('Name',['34-43 Hz spatial p values - ' contrastName]);
    stem(1:nCh,-log10(pRaw),'filled');
    hold on;
    yline(-log10(0.05),'--','raw p=.05');
    sig=find(pFDR<0.05);
    if ~isempty(sig)
        scatter(sig,-log10(pRaw(sig)),70,'filled');
    end
    xticks(1:nCh);
    xticklabels(CHANNELS);
    xtickangle(45);
    ylabel('-log10(raw p)');
    title(sprintf('Channel-wise Welch tests; BH-FDR across 30: %s',strrep(contrastName,'_',' ')));
    box on; hold off;
    saveas(gcf,sprintf('Spatial_34_43Hz_%s_pvalues.png',contrastName));

    % Optional EEGLAB topoplot if locks3.loc and readlocs/topoplot exist.
    % The statistics DO NOT depend on this plot.
    if exist('topoplot','file')==2 && exist('readlocs','file')==2 && isfile('locks3.loc')
        try
            chanlocs=readlocs('locks3.loc');
            if numel(chanlocs)==30
                % R is sorted; use the original channel-order g vector.
                figure('Name',['34-43 Hz scalp map - ' contrastName]);
                topoplot(g,chanlocs,'electrodes','on');
                colorbar;
                title(sprintf('Hedges g, %s',strrep(contrastName,'_',' ')));
                saveas(gcf,sprintf('Spatial_34_43Hz_%s_topoplot.png',contrastName));
            end
        catch ME
            warning('Topoplot skipped: %s',ME.message);
        end
    end
end


%% -------------------------- LOCAL FUNCTIONS -----------------------------
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

function g=hedges_g(E,C)
    nE=numel(E); nC=numel(C);
    sp=sqrt(((nE-1)*var(E)+(nC-1)*var(C))/(nE+nC-2));
    if sp==0, g=NaN; return; end
    d=(mean(E)-mean(C))/sp;
    g=(1-3/(4*(nE+nC)-9))*d;
end

function q=bh_fdr(p)
    p=p(:);
    [ps,ord]=sort(p);
    m=numel(ps);
    qs=ps.*m./(1:m)';
    qs=flipud(cummin(flipud(qs)));
    qs=min(qs,1);
    q=nan(size(p));
    q(ord)=qs;
end
