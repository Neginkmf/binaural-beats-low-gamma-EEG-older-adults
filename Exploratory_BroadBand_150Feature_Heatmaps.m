%% EXPLORATORY BROAD-BAND: 150-feature Welch tests + heatmaps
% Input: 05_AllBands_Changes_dB_vs_PRE.xlsx
%
% 30 channels x 5 bands = 150 exploratory features.
% Runs DURING_vs_PRE and POST_vs_PRE separately.
%
% For each contrast:
%   - Welch t-test for every feature
%   - Hedges g
%   - BH-FDR across ALL 150 features
%   - raw-p heatmap (channels x frequency bands)
%   - effect-size heatmap
%
% NOTE: raw p-value stars are exploratory. Inferential conclusions should
% use the FDR-adjusted q values.

clear; clc; close all;

CHANNELS = {'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4','P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4','TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};
BANDS = {'Delta','Theta','Alpha','Beta','Gamma'};
CONTRASTS = {'DURING_vs_PRE','POST_vs_PRE'};


filename=resolve_file('05_AllBands_Changes_dB_vs_PRE.xlsx');
Tall=readtable(filename,'VariableNamingRule','preserve');

% Build and VERIFY the exact expected 150-column order.
expectedVars=cell(1,numel(CHANNELS)*numel(BANDS));
k=1;
for ch=1:numel(CHANNELS)
    for b=1:numel(BANDS)
        expectedVars{k}=[CHANNELS{ch} '_' BANDS{b}];
        k=k+1;
    end
end

actualVars=Tall.Properties.VariableNames(6:end);
assert(numel(actualVars)==150,'Expected exactly 150 broad-band features.');
assert(isequal(actualVars,expectedVars), ...
    ['150-FEATURE ORDER CHECK FAILED. Expected channel-major order ' ...
     '(Fp1 Delta..Gamma, Fp2 Delta..Gamma, ... Oz Delta..Gamma).']);

fprintf('150-feature order verification PASSED.\n');
fprintf('Channel order:\n%s\n\n',strjoin(CHANNELS,' -> '));

for cc = 1:numel(CONTRASTS)

    contrastName = CONTRASTS{cc};

    % Select this contrast
    T = Tall(strcmp(string(Tall.Contrast), contrastName), :);

    % ---------------------------------------------------------
    % CONTRAST-SPECIFIC EXCLUSIONS
    % ---------------------------------------------------------
    if strcmp(contrastName,'DURING_vs_PRE')

        % E15: bad PRE
        % C02: DURING-specific sensitivity exclusion
        EXCLUDE_SUBJECTS = {'E_15','C_02', 'C_11', 'E_17'};

    elseif strcmp(contrastName,'POST_vs_PRE')

        % C02 remains included because its issue is DURING only
        EXCLUDE_SUBJECTS = {'E_15'};

    else
        EXCLUDE_SUBJECTS = {};
    end

    % Apply exclusions
    T = apply_exclusions(T,EXCLUDE_SUBJECTS);

    fprintf('\nContrast: %s\n',contrastName);
    fprintf('Excluded subjects: %s\n', ...
        strjoin(EXCLUDE_SUBJECTS,', '));

    labels=double(T.Label);
    X=table2array(T(:,expectedVars));

    pRaw=nan(150,1);
    tstat=nan(150,1);
    df=nan(150,1);
    g=nan(150,1);
    meanC=nan(150,1);
    meanE=nan(150,1);

    for f=1:150
        C=X(labels==0,f);
        E=X(labels==1,f);
        [~,p,~,st]=ttest2(E,C,'Vartype','unequal');
        pRaw(f)=p;
        tstat(f)=st.tstat;
        df(f)=st.df;
        meanC(f)=mean(C);
        meanE(f)=mean(E);
        g(f)=hedges_g(E,C);
    end

    pFDR=bh_fdr(pRaw);
    diffMean=meanE-meanC;

    R=table(string(expectedVars(:)),meanC,meanE,diffMean,tstat,df,pRaw,pFDR,g, ...
        'VariableNames',{'Feature','Mean_Control_dB','Mean_Experiment_dB', ...
        'Difference_EminusC_dB','Welch_t','Welch_df','p_raw','p_FDR_BH','Hedges_g'});
    R=sortrows(R,'p_raw');

    fprintf('\n=== Broad-band exploratory analysis: %s ===\n',contrastName);
    disp(R(1:15,:));
    fprintf('Raw p < .05: %d / 150\n',sum(pRaw<0.05));
    fprintf('BH-FDR q < .05: %d / 150\n',sum(pFDR<0.05));

    writetable(R,sprintf('BroadBand_150features_%s_results.csv',contrastName));

    % The feature order is channel-major, band-minor.
    pGrid=reshape(pRaw,[numel(BANDS),numel(CHANNELS)])';
    qGrid=reshape(pFDR,[numel(BANDS),numel(CHANNELS)])';
    gGrid=reshape(g,[numel(BANDS),numel(CHANNELS)])';

    figure('Name',['Broad-band raw p heatmap - ' contrastName]);
    imagesc(pGrid);
    colorbar;
    xlabel('Frequency band'); ylabel('Channel');
    xticks(1:numel(BANDS)); xticklabels(BANDS);
    yticks(1:numel(CHANNELS)); yticklabels(CHANNELS);
    title(sprintf('Raw Welch p-values (exploratory): %s',strrep(contrastName,'_',' ')));
    hold on;
    for ch=1:numel(CHANNELS)
        for b=1:numel(BANDS)
            if qGrid(ch,b)<0.05
                text(b,ch,'F','HorizontalAlignment','center','FontWeight','bold');
            elseif pGrid(ch,b)<0.05
                text(b,ch,'*','HorizontalAlignment','center');
            end
        end
    end
    hold off;
    saveas(gcf,sprintf('BroadBand_%s_rawP_heatmap.png',contrastName));

    figure('Name',['Broad-band Hedges g heatmap - ' contrastName]);
    imagesc(gGrid);
    colorbar;
    xlabel('Frequency band'); ylabel('Channel');
    xticks(1:numel(BANDS)); xticklabels(BANDS);
    yticks(1:numel(CHANNELS)); yticklabels(CHANNELS);
    title(sprintf('Hedges g (Experiment-Control): %s',strrep(contrastName,'_',' ')));
    saveas(gcf,sprintf('BroadBand_%s_HedgesG_heatmap.png',contrastName));
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
