function bandPower = func_FE_IntegrateBandPower(PSD, freqs, bandEdges)
%%---------------------------------------------------------------------------------------------------------
%%% Integrates absolute band power from an already-computed Welch PSD.
%%% Uses trapezoidal integration and includes both lower and upper band
%%% boundaries when those frequency bins exist: low <= f <= high.
%%---------------------------------------------------------------------------------------------------------

if size(bandEdges,2) ~= 2
    error('bandEdges must be an N x 2 matrix of [low high] limits.');
end
if size(PSD,2) ~= numel(freqs)
    error('PSD column count must match the number of frequency bins.');
end

nBands = size(bandEdges,1);
nCh = size(PSD,1);
bandPower = nan(nCh,nBands);

for b = 1:nBands
    low  = bandEdges(b,1);
    high = bandEdges(b,2);

    if low < 0 || high <= low
        error('Invalid band limits [%.3f %.3f].',low,high);
    end

    inBand = freqs >= low & freqs <= high;

    if sum(inBand) < 2
        error(['Band %.2f-%.2f Hz contains fewer than two frequency bins. ' ...
               'Increase Welch window duration or check sampling rate.'], low, high);
    end

    bandPower(:,b) = trapz(freqs(inBand),PSD(:,inBand),2);
end
end
