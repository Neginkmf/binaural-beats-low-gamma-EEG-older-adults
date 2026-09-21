%% Extract_EEG_PRE_DURING_POST_FFT_Comparable.m
% =========================================================================
% LEGACY / COMPARISON PIPELINE:
% Whole-record FFT band-power extraction, made directly comparable to the
% current Welch pipeline.
%
% PURPOSE
% -------
% This script is for a METHOD COMPARISON with the Welch extraction.
% It keeps the following identical to the Welch pipeline:
%   - exact same PRE / DURING / POST files
%   - exact same subject/group/montage assignments
%   - exact same final 30-channel canonical order
%   - exact same common-average reference
%   - exact same broad bands and 34-43 Hz target band
%   - exact same dB change definition:
%         10*log10(DURING/PRE)
%         10*log10(POST/PRE)
%
% The ONLY intended major difference is the spectral estimator:
%   WELCH pipeline:
%       4-s Hann-window PSD + window averaging + PSD integration
%
%   THIS pipeline:
%       one FFT over the entire recording, rectangular frequency mask,
%       inverse FFT reconstruction, then mean-square band energy.
%
% IMPORTANT
% ---------
% This reproduces the FIXED version of the older FFT method:
% power is normalized by the REAL recording length L, not an NFFT-dependent
% length. Therefore it does NOT reproduce the old zero-padding normalization
% bug.
%
% This is a comparison/sensitivity pipeline, not a recommendation to replace
% the Welch pipeline simply because one method yields smaller p-values.
%
% OUTPUT FILES
% ------------
% 00_FileManifest_FFT.xlsx
% 01_AllBands_RawPower_PRE_DURING_POST.xlsx
% 02_NarrowGamma_34_43Hz_RawPower_PRE_DURING_POST.xlsx
% 03_AllBands_LogPower_dB_PRE_DURING_POST.xlsx
% 04_NarrowGamma_34_43Hz_LogPower_dB_PRE_DURING_POST.xlsx
% 05_AllBands_Changes_dB_vs_PRE.xlsx
% 06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx
% 07_FFT_QC.xlsx
% EEG_FFT_AllResults.mat
%
% BEST PRACTICE FOR FAIR COMPARISON
% ---------------------------------
% Point manifestFile to the 00_FileManifest.xlsx produced by your FINAL
% corrected Welch extraction. That guarantees both methods process the
% exact same files with the exact same montage labels.
% =========================================================================

clear; clc;

%% =========================== USER SETTINGS ==============================
Fs = 500;
doAverageReference = true;

% EDIT THIS:
% Select the FINAL 00_FileManifest.xlsx from your Welch extraction
[fileName, folderName] = uigetfile( ...
    '*.xlsx', ...
    'Select FINAL Welch 00_FileManifest.xlsx');

if isequal(fileName,0)
    error('No manifest selected.');
end

manifestFile = fullfile(folderName,fileName);

% Create a new FFT-comparison folder beside the Welch results
outputDir = fullfile(folderName,'FFT_Extraction_COMPARABLE');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

fprintf('\nManifest:\n%s\n',manifestFile);
fprintf('\nFFT results will be saved in:\n%s\n\n',outputDir);


% Subjects intentionally excluded from the final PRE-based analyses.
% These rows are removed from the manifest BEFORE checking file existence.
SOURCE_EXCLUDE_SUBJECTS = {'E_15'};

%% ========================= CHANNEL DEFINITIONS ==========================
% VERIFIED common 30-channel order used by the current extraction files.
canonicalLabels = { ...
    'Fp1','Fp2','F3','F4','FC3','FC4','C3','C4','CP3','CP4', ...
    'P3','P4','O1','O2','F7','F8','FT7','FT8','T3','T4', ...
    'TP7','TP8','T5','T6','Fz','FCz','Cz','CPz','Pz','Oz'};

locks3Labels = canonicalLabels;

% locks4 order BEFORE removing FPz:
locks4FullLabels = { ...
    'Fp1','F3','FC3','C3','Fz','FCz','Cz','FT8','T4','TP8', ...
    'T6','F8','FPz','O2','P4','Fp2','FC4','C4','CP4','Pz', ...
    'F4','Oz','CPz','F7','O1','P3','CP3','FT7','T3','TP7','T5'};

locks4NoFpzLabels = locks4FullLabels(~strcmpi(locks4FullLabels,'FPz'));

bandNames  = {'Delta','Theta','Alpha','Beta','Gamma'};
bandLimits = [0.5 4; 4 8; 8 13; 13 30; 30 80];

narrowLimits = [34 43];

%% ========================== LOAD MANIFEST ===============================
if ~isfile(manifestFile)
    error(['Manifest not found:\n%s\n\nPoint manifestFile to the ' ...
           '00_FileManifest.xlsx from the FINAL Welch extraction.'], ...
           manifestFile);
end

manifest = readtable(manifestFile,'VariableNamingRule','preserve');

requiredCols = {'SubjectUID','SubjectID','Group','Label','Time','Montage','FilePath'};
for k = 1:numel(requiredCols)
    assert(ismember(requiredCols{k},manifest.Properties.VariableNames), ...
        'Manifest is missing required column: %s',requiredCols{k});
end

% Convert relevant columns to strings for stable comparisons.
manifest.SubjectUID = string(manifest.SubjectUID);
manifest.SubjectID  = string(manifest.SubjectID);
manifest.Group      = string(manifest.Group);
manifest.Time       = string(manifest.Time);
manifest.Montage    = string(manifest.Montage);
manifest.FilePath   = string(manifest.FilePath);

% Remove intentionally excluded subjects BEFORE duplicate/file checks.
if ~isempty(SOURCE_EXCLUDE_SUBJECTS)
    rm = ismember(manifest.SubjectUID,string(SOURCE_EXCLUDE_SUBJECTS));
    if any(rm)
        fprintf('\nRemoving intentionally excluded subject(s) from FFT source manifest: %s\n', ...
            strjoin(cellstr(string(unique(manifest.SubjectUID(rm),'stable'))),', '));
        manifest = manifest(~rm,:);
    end
end

% Verify no duplicate subject/time sessions.
key = manifest.SubjectUID + "|" + upper(manifest.Time);
assert(numel(unique(key)) == height(manifest), ...
    'Manifest contains duplicate SubjectUID x Time rows.');

% Verify file existence before starting.
missing = ~arrayfun(@(x)isfile(char(x)), manifest.FilePath);
if any(missing)
    fprintf('\nMissing files in manifest:\n');
    disp(manifest(missing,{'SubjectUID','Time','FilePath'}));
    error('Fix missing FilePath entries before running.');
end

writetable(manifest,fullfile(outputDir,'00_FileManifest_FFT.xlsx'));

fprintf('Loaded manifest with %d sessions.\n',height(manifest));
fprintf('This script will use the SAME files/montages as the Welch extraction.\n\n');

%% ============================= EXTRACTION ===============================
nSessions = height(manifest);
nCh = numel(canonicalLabels);
nBands = numel(bandNames);

allBandRaw = nan(nSessions,nCh*nBands);
narrowRaw  = nan(nSessions,nCh);

durationSec = nan(nSessions,1);
nSamples = nan(nSessions,1);
nFFTused = nan(nSessions,1);

for i = 1:nSessions

    fprintf('\n[%d/%d] %s | %s | %s\n', ...
        i,nSessions,manifest.SubjectUID(i), ...
        manifest.Time(i),manifest.Montage(i));

    % -------- Load exactly the same preprocessed ASCII EEG --------
    X = readAsciiEEG(char(manifest.FilePath(i)), ...
                     char(manifest.Montage(i)));

    % -------- Map both montages to the exact same 30 channels -----
    X = standardizeMontage(X,char(manifest.Montage(i)), ...
        canonicalLabels,locks3Labels, ...
        locks4FullLabels,locks4NoFpzLabels);

    assert(size(X,1)==30, ...
        'Final standardized EEG must contain exactly 30 channels.');

    % -------- Same reference as Welch pipeline --------------------
    if doAverageReference
        X = X - mean(X,1,'omitnan');
    end

    if any(~isfinite(X(:)))
        error('NaN/Inf found after standardization: %s', ...
            manifest.FilePath(i));
    end

    nSamples(i) = size(X,2);
    durationSec(i) = nSamples(i)/Fs;
    nFFTused(i) = 2^nextpow2(nSamples(i));

    % -------- Legacy fixed whole-record FFT broad-band power ------
    B = wholeRecordFFT_BroadBands(X,Fs,bandLimits);
    % B = channels x 5

    % Flatten channel-major:
    % Fp1_Delta ... Fp1_Gamma, Fp2_Delta ... etc.
    allBandRaw(i,:) = reshape(B.',1,[]);

    % -------- Legacy fixed whole-record FFT 34-43 Hz --------------
    narrowRaw(i,:) = wholeRecordFFT_SingleBand( ...
        X,Fs,narrowLimits).';

    fprintf('  %d ch x %d samples | %.3f s | NFFT=%d\n', ...
        size(X,1),size(X,2),durationSec(i),nFFTused(i));
end

%% ============================ FEATURE NAMES =============================
broadFeatureNames = cell(1,nCh*nBands);
k = 1;

for ch = 1:nCh
    for b = 1:nBands
        broadFeatureNames{k} = matlab.lang.makeValidName( ...
            [canonicalLabels{ch} '_' bandNames{b}]);
        k = k+1;
    end
end

narrowFeatureNames = cellfun( ...
    @(x)matlab.lang.makeValidName([x '_Gamma34_43']), ...
    canonicalLabels,'UniformOutput',false);

% Explicit order check.
expectedFirstFive = {'Fp1_Delta','Fp1_Theta','Fp1_Alpha','Fp1_Beta','Fp1_Gamma'};
assert(isequal(broadFeatureNames(1:5),expectedFirstFive), ...
    'Unexpected feature flattening order.');

fprintf('\nFeature/channel order verification PASSED.\n');
fprintf('%s\n',strjoin(canonicalLabels,' -> '));

%% ========================== SESSION OUTPUTS =============================
% ----------------------------- RAW POWER --------------------------------
broadRawTable = [manifest(:,1:6), ...
    array2table(allBandRaw,'VariableNames',broadFeatureNames)];

narrowRawTable = [manifest(:,1:6), ...
    array2table(narrowRaw,'VariableNames',narrowFeatureNames)];

narrowRawTable.WholeScalp_MeanRaw = mean(narrowRaw,2,'omitnan');

writetable(broadRawTable,fullfile(outputDir, ...
    '01_AllBands_RawPower_PRE_DURING_POST.xlsx'));

writetable(narrowRawTable,fullfile(outputDir, ...
    '02_NarrowGamma_34_43Hz_RawPower_PRE_DURING_POST.xlsx'));

% --------------------------- ABSOLUTE LOG POWER --------------------------
allBandLogdB = 10*log10(allBandRaw);
narrowLogdB  = 10*log10(narrowRaw);

broadLogTable = [manifest(:,1:6), ...
    array2table(allBandLogdB,'VariableNames',broadFeatureNames)];

narrowLogTable = [manifest(:,1:6), ...
    array2table(narrowLogdB,'VariableNames',narrowFeatureNames)];

% Same definition as Welch:
% mean of CHANNEL-WISE log powers.
narrowLogTable.WholeScalp_Avg_dB = ...
    mean(narrowLogdB,2,'omitnan');

writetable(broadLogTable,fullfile(outputDir, ...
    '03_AllBands_LogPower_dB_PRE_DURING_POST.xlsx'));

writetable(narrowLogTable,fullfile(outputDir, ...
    '04_NarrowGamma_34_43Hz_LogPower_dB_PRE_DURING_POST.xlsx'));

%% =========================== CHANGE SCORES ==============================
% Same definition as Welch:
%
% target log-power - PRE log-power
% = 10*log10(target/PRE)
%
% This is deliberately NOT the old raw POST-PRE subtraction.

broadChange = makeChangeTable( ...
    broadLogTable,broadFeatureNames,{'DURING','POST'},'PRE');

narrowChange = makeChangeTable( ...
    narrowLogTable,narrowFeatureNames,{'DURING','POST'},'PRE');

if ~isempty(narrowChange)
    narrowChange.WholeScalp_Avg_dB = ...
        mean(narrowChange{:,narrowFeatureNames},2,'omitnan');
end

writetable(broadChange,fullfile(outputDir, ...
    '05_AllBands_Changes_dB_vs_PRE.xlsx'));

writetable(narrowChange,fullfile(outputDir, ...
    '06_NarrowGamma_34_43Hz_Changes_dB_vs_PRE.xlsx'));

%% ================================ QC ====================================
% Whole-record FFT has no Welch-window rejection.
% These fields make the methodological difference explicit.

QC = manifest(:,1:7);
QC.Duration_sec = durationSec;
QC.NSamples = nSamples;
QC.NFFT = nFFTused;
QC.Method = repmat("WholeRecordFFT_Fixed",height(QC),1);
QC.WindowRejectionApplied = false(height(QC),1);

writetable(QC,fullfile(outputDir,'07_FFT_QC.xlsx'));

%% ================================ SAVE ==================================
save(fullfile(outputDir,'EEG_FFT_AllResults.mat'), ...
    'manifest','canonicalLabels','bandNames','bandLimits','narrowLimits', ...
    'allBandRaw','allBandLogdB','narrowRaw','narrowLogdB', ...
    'broadChange','narrowChange','QC','Fs', ...
    'doAverageReference','-v7.3');

fprintf('\n============================================================\n');
fprintf('FFT COMPARISON EXTRACTION COMPLETE\n');
fprintf('============================================================\n');
fprintf('Output folder:\n%s\n\n',outputDir);
fprintf('Estimator: fixed whole-record FFT + rectangular band mask\n');
fprintf('Reference: common average over the same 30 channels = %d\n', ...
    doAverageReference);
fprintf('Change metric: 10*log10(target/PRE), same as Welch\n');
fprintf('Broad bands: Delta, Theta, Alpha, Beta, Gamma\n');
fprintf('Target band: 34-43 Hz\n');
fprintf('============================================================\n');


%% ============================ LOCAL FUNCTIONS ===========================

function X = readAsciiEEG(filePath,montage)
% Read the exported preprocessed numeric matrix exactly as text regardless
% of extension (.mat may actually be ASCII in this dataset).

    try
        M = readmatrix(filePath,'FileType','text','Delimiter','\t');
    catch ME
        error('Could not read %s as tab-delimited ASCII: %s', ...
            filePath,ME.message);
    end

    M(all(isnan(M),2),:) = [];
    M(:,all(isnan(M),1)) = [];

    if isempty(M)
        error('No numeric EEG data found in %s',filePath);
    end

    if any(isnan(M(:)))
        error(['Partial NaNs/non-numeric cells in %s. Expected a pure ' ...
               'numeric EEG matrix with no text/time/index column.'], ...
               filePath);
    end

    if strcmpi(montage,'locks3')
        allowed = 30;
    elseif strcmpi(montage,'locks4')
        allowed = [30 31];
    else
        error('Unknown montage: %s',montage);
    end

    rowMatch = ismember(size(M,1),allowed);
    colMatch = ismember(size(M,2),allowed);

    if rowMatch && ~colMatch
        X = double(M);         % channels x samples
    elseif colMatch && ~rowMatch
        X = double(M.');       % transpose
    elseif rowMatch && colMatch
        error('Ambiguous matrix shape %dx%d in %s.', ...
            size(M,1),size(M,2),filePath);
    else
        error(['Unexpected EEG size %dx%d in %s. Expected 30 channels ' ...
               'for locks3 or 30/31 for locks4.'], ...
               size(M,1),size(M,2),filePath);
    end
end


function Xout = standardizeMontage( ...
    X,montage,canonical,locks3,locks4full,locks4noFpz)
% Map the preprocessed EEG to the VERIFIED common 30-channel order.

    if strcmpi(montage,'locks3')

        if size(X,1)~=30
            error('locks3 input must contain exactly 30 channels.');
        end
        sourceLabels = locks3;

    elseif strcmpi(montage,'locks4')

        if size(X,1)==31
            sourceLabels = locks4full;

            % Remove FPz before mapping to common 30.
            fpz = find(strcmpi(sourceLabels,'FPz'));
            assert(numel(fpz)==1,'Could not uniquely identify FPz.');
            X(fpz,:) = [];
            sourceLabels(fpz) = [];

        elseif size(X,1)==30
            sourceLabels = locks4noFpz;

        else
            error('locks4 input must contain 30 or 31 channels.');
        end

    else
        error('Unknown montage: %s',montage);
    end

    [tf,idx] = ismember( ...
        lower(string(canonical)),lower(string(sourceLabels)));

    if ~all(tf)
        error('Could not map all 30 canonical channels.');
    end

    Xout = X(idx,:);
end


function B = wholeRecordFFT_BroadBands(X,Fs,bandLimits)
% -------------------------------------------------------------------------
% Reproduces the FIXED legacy whole-record FFT method.
%
% For each channel:
%   - one FFT over entire recording
%   - zero pad to 2^nextpow2(L)
%   - keep positive-frequency coefficients in desired band
%   - inverse FFT with 'symmetric'
%   - band power = sum(reconstructed_signal.^2) / REAL L
%
% IMPORTANT:
% Division by L is the FIX that avoids the old NFFT-dependent normalization
% jump.
%
% Output:
%   B = channels x bands
% -------------------------------------------------------------------------

    [nCh,L] = size(X);
    nBands = size(bandLimits,1);
    NFFT = 2^nextpow2(L);

    f = Fs/2*linspace(0,1,NFFT/2+1);

    B = nan(nCh,nBands);

    for ch = 1:nCh

        x = X(ch,:);

        % Same normalization convention used in the fixed old function.
        F = fft(x,NFFT)/L;

        for b = 1:nBands

            low  = bandLimits(b,1);
            high = bandLimits(b,2);

            % Match legacy boundary convention:
            % low inclusive, high exclusive.
            mask = f>=low & f<high;

            Fband = zeros(1,NFFT/2+1);
            Fband(mask) = F(mask)*L;

            bandSignal = ifft(Fband,NFFT,'symmetric');

            % FIXED normalization by REAL sample count.
            B(ch,b) = sum(bandSignal.^2)/L;
        end
    end

    B(B<=0 | ~isfinite(B)) = NaN;
end


function power = wholeRecordFFT_SingleBand(X,Fs,band)
% Same fixed legacy method for a single arbitrary band such as 34-43 Hz.

    [nCh,L] = size(X);
    NFFT = 2^nextpow2(L);

    f = Fs/2*linspace(0,1,NFFT/2+1);

    power = nan(nCh,1);

    for ch = 1:nCh

        x = X(ch,:);
        F = fft(x,NFFT)/L;

        mask = f>=band(1) & f<band(2);

        Fband = zeros(1,NFFT/2+1);
        Fband(mask) = F(mask)*L;

        bandSignal = ifft(Fband,NFFT,'symmetric');

        power(ch) = sum(bandSignal.^2)/L;
    end

    power(power<=0 | ~isfinite(power)) = NaN;
end


function Tout = makeChangeTable( ...
    Tlog,featureNames,targetTimes,baselineTime)
% Create target-vs-PRE dB changes using exactly the same logic as the
% current Welch extraction.

    Tout = table();
    uids = unique(Tlog.SubjectUID,'stable');

    for u = 1:numel(uids)

        Tu = Tlog(Tlog.SubjectUID==uids(u),:);

        baseIdx = strcmpi(Tu.Time,baselineTime);

        if sum(baseIdx)~=1
            warning('%s does not have exactly one %s session. Skipping.', ...
                uids(u),baselineTime);
            continue;
        end

        baseVals = Tu{baseIdx,featureNames};

        for t = 1:numel(targetTimes)

            targetIdx = strcmpi(Tu.Time,targetTimes{t});

            if sum(targetIdx)~=1
                warning('%s does not have exactly one %s session. Skipping contrast.', ...
                    uids(u),targetTimes{t});
                continue;
            end

            dBchange = Tu{targetIdx,featureNames} - baseVals;

            first = Tu(find(baseIdx,1),:);

            meta = table( ...
                first.SubjectUID, ...
                first.SubjectID, ...
                first.Group, ...
                first.Label, ...
                string([targetTimes{t} '_vs_' baselineTime]), ...
                'VariableNames', ...
                {'SubjectUID','SubjectID','Group','Label','Contrast'});

            Tout = [Tout; ...
                meta, ...
                array2table(dBchange,'VariableNames',featureNames)]; %#ok<AGROW>
        end
    end
end
