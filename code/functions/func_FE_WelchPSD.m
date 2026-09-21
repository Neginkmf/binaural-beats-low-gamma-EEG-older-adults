function [PSD, freqs, rejectedFrac, nWindowsUsed, nWindowsTotal] = func_FE_WelchPSD(data, Fs, winSec, threshMAD)
%%---------------------------------------------------------------------------------------------------------
%%% Welch-style PSD per channel using fixed-length non-overlapping windows
%%% plus robust segment-level artifact rejection.
%%%
%%% IMPORTANT: if fewer than 4 clean windows survive, this function keeps
%%% the QC information and returns NaN PSD instead of silently reusing all
%%% rejected windows.
%%---------------------------------------------------------------------------------------------------------

if nargin < 3 || isempty(winSec);    winSec = 4; end
if nargin < 4 || isempty(threshMAD); threshMAD = 5; end

data = double(data);
if any(~isfinite(data(:))); error('EEG contains NaN or Inf values.'); end

[nCh, L] = size(data);
winLen = round(winSec * Fs);
if winLen < 2; error('winSec is too short for the supplied sampling rate.'); end

nWin = floor(L / winLen);
nWindowsTotal = nWin;
if nWin < 4
    error('Recording too short for %.1f-second windows: only %d complete windows available.', winSec, nWin);
end

data = data(:,1:nWin*winLen);
segs = reshape(data, nCh, winLen, nWin);

%% ---------------------- SEGMENT-LEVEL ARTIFACT QC -----------------------
ptp = squeeze(max(segs,[],2) - min(segs,[],2));
if nCh == 1; ptp = ptp'; end

badPerChannel = false(nCh,nWin);
for ch = 1:nCh
    medPTP = median(ptp(ch,:));
    madPTP = mad(ptp(ch,:),1);
    robustScale = 1.4826 * madPTP;

    if robustScale == 0
        badPerChannel(ch,:) = false;
    else
        threshold = medPTP + threshMAD * robustScale;
        badPerChannel(ch,:) = ptp(ch,:) > threshold;
    end
end

fracChannelsBad = mean(badPerChannel,1);
goodWin = fracChannelsBad <= 0.20;

nWindowsUsed = sum(goodWin);
rejectedFrac = 1 - nWindowsUsed/nWin;

%% ------------------------- FREQUENCY AXIS -------------------------------
NFFT = winLen;
nFreqBins = floor(NFFT/2) + 1;
freqs = (0:nFreqBins-1) * (Fs/NFFT);

if nWindowsUsed < 4
    warning(['Only %d of %d clean windows survived (%.1f%% rejected). ' ...
             'PSD returned as NaN; inspect this recording manually.'], ...
             nWindowsUsed, nWin, 100*rejectedFrac);
    PSD = nan(nCh,nFreqBins);
    return;
end

%% ------------------------- WELCH AVERAGING ------------------------------
hannWin = 0.5 * (1 - cos(2*pi*(0:winLen-1)'/(winLen-1)));
U = sum(hannWin.^2);
PSD = zeros(nCh,nFreqBins);
goodIdx = find(goodWin);

for ch = 1:nCh
    accum = zeros(1,nFreqBins);

    for w = goodIdx
        x = squeeze(segs(ch,:,w))';
        x = x - mean(x);
        xw = x .* hannWin;
        Xw = fft(xw,NFFT);

        Pxx = (1/(Fs*U)) * abs(Xw).^2;
        Pxx = Pxx(1:nFreqBins);

        if rem(NFFT,2) == 0
            if nFreqBins > 2; Pxx(2:end-1) = 2*Pxx(2:end-1); end
        else
            if nFreqBins > 1; Pxx(2:end) = 2*Pxx(2:end); end
        end

        accum = accum + Pxx';
    end

    PSD(ch,:) = accum / nWindowsUsed;
end
end
