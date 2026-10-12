function [saida, info] = PE_func(timeSeries, tauList, D, varargin)
% PE_func  Bandt-Pompe permutation entropy, Jensen-Shannon statistical
% complexity and Fisher information, computed only from ordinal vectors
% that lie entirely inside one continuous recording segment.
%
%   saida = PE_func(x, tauList, D)
%   saida = PE_func(x, tauList, D, 'MinZeroRun', 400)
%   saida = PE_func(x, tauList, D, 'ValidMask', mask)
%   [saida, info] = PE_func(...)
%
% IMPORTANT: x must be passed WITH the zero-amplitude buffers still in
% place (do NOT remove the zeros before calling). The buffers are used to
% locate the boundaries between the originally separated segments.
%
% Inputs
%   x          time series (vector)
%   tauList    vector of embedding delays
%   D          embedding dimension (D <= 9)
%
% Options (name-value)
%   'RespectSegments'  true (default) = discard ordinal vectors that cross
%                a segment boundary. false = use ALL windows of the series
%                exactly as given (same result as the original function;
%                'SegmentId', 'ValidMask' and 'MinZeroRun' are ignored).
%   'SegmentId'  integer vector, same length as x; samples with the same
%                label belong to the same continuous segment, 0 = invalid.
%                A new segment starts whenever the label changes. If given,
%                it overrides 'ValidMask' and automatic zero detection.
%   'ValidMask'  logical vector, same length as x; true = recorded sample.
%                If given, it is used instead of automatic zero detection.
%   'MinZeroRun' minimum length of a run of exact zeros to be treated as
%                a buffer (default 10 samples). Isolated zeros shorter
%                than this are kept as data. NaN samples are always
%                treated as invalid.
%
% An ordinal vector (x_j, x_{j+tau}, ..., x_{j+(D-1)tau}) is used only if
% all samples from j to j+(D-1)tau belong to the same continuous segment.
% Vectors that would cross a buffer are discarded, so no artificial
% transitions between discontinuous portions of the recording are created.
%
% Output
%   saida  [length(tauList) x 3] = [H, C, F]
%   info   struct with segment and window counts (see end of main function)

timeSeries = timeSeries(:);

p = inputParser;
addParameter(p, 'RespectSegments', true);
addParameter(p, 'SegmentId', []);
addParameter(p, 'ValidMask', []);
addParameter(p, 'MinZeroRun', 10);
parse(p, varargin{:});

if p.Results.RespectSegments
    [segId, nSeg] = segment_labels(timeSeries, p.Results.SegmentId, ...
        p.Results.ValidMask, p.Results.MinZeroRun);
else
    % All samples treated as one continuous segment: every window is used,
    % exactly as in the original version of this function.
    segId = ones(numel(timeSeries), 1);
    nSeg  = 1;
end
nValidSamples = sum(segId > 0);

saida = zeros(length(tauList), 3);

nWinUsed     = zeros(length(tauList), 1);
nWinExcluded = zeros(length(tauList), 1);

N = factorial(D);        % Number of possible ordinal patterns
S_max = log2(N);         % Maximum entropy

%   >>>>>   SYMBOL/LEHMER   <<<<<
% Generate all possible ordinal patterns for embedding dimension D.
all_symbols = flipud(perms(1:D));

% Powers used to convert each ordinal pattern into a unique numeric symbol.
% Example: [1 2 3] becomes 123.
powers = 10.^(D-1:-1:0);

% Convert all ordinal patterns into numeric symbols.
all_symbols_flat = sum(all_symbols .* repmat(powers, [N, 1]), 2);

% Create a lookup table that maps each numeric symbol to its Bandt-Pompe index.
index_map = zeros(max(all_symbols_flat), 1);

for i = 1:N
    symbol = all_symbols_flat(i);
    index_map(symbol) = i;
end

tauId = 0;

for tau = tauList

    tauId = tauId + 1;

    % Compute the Bandt-Pompe probability distribution
    % (only windows fully contained in one segment).
    [P, nUsed] = PDF_segments(timeSeries, segId, D, tau, 1, index_map, powers, N);

    % Compute the normalized Shannon entropy.
    H = Shannon(P) ./ S_max;

    % Compute the Jensen-Shannon statistical complexity.
    C = Complexity_Jensen_Shannon(P, N, H);

    % Compute the Fisher information measure.
    F = Fisher(P);

    % Store entropy, complexity, and Fisher information for the current tau.
    saida(tauId, :) = [H, C, F];

    % Windows that the old approach (removing the zeros and concatenating)
    % would have formed across segment boundaries.
    nWinUsed(tauId)     = nUsed;
    nWinExcluded(tauId) = max(nValidSamples - (D-1)*tau, 0) - nUsed;

end

info.RespectSegments    = logical(p.Results.RespectSegments);
info.nSegments          = nSeg;            % continuous segments found
info.nValidSamples      = nValidSamples;   % samples of recorded activity
info.tauList            = tauList(:);
info.nWindowsUsed       = nWinUsed;        % ordinal vectors used, per tau
info.nWindowsExcluded   = nWinExcluded;    % boundary-crossing vectors discarded, per tau

end


%% FUNCTIONS

%   >>>>>   SEGMENT LABELS   <<<<<
function [segId, nSeg] = segment_labels(x, segIn, validMask, minZeroRun)
% segId(t) = k > 0 if sample t belongs to the k-th continuous segment,
% segId(t) = 0 if sample t is a buffer / invalid sample.

M = numel(x);

if ~isempty(segIn)

    segIn = segIn(:);

    if numel(segIn) ~= M
        error('PE_func:SegmentId', 'SegmentId must have the same length as the time series.');
    end

    valid = segIn > 0 & ~isnan(x);

    % New segment when the label changes or after an invalid sample.
    segStart = valid & [true; (segIn(2:end) ~= segIn(1:end-1)) | ~valid(1:end-1)];
    segId    = cumsum(segStart);
    segId(~valid) = 0;
    nSeg     = sum(segStart);
    return;

end

if isempty(validMask)

    valid = ~isnan(x);

    % Mark runs of exact zeros with length >= minZeroRun as buffers.
    isZero   = (x == 0);
    d        = diff([false; isZero; false]);
    runStart = find(d == 1);
    runEnd   = find(d == -1) - 1;
    runLen   = runEnd - runStart + 1;

    for r = find(runLen >= minZeroRun).'
        valid(runStart(r):runEnd(r)) = false;
    end

else

    valid = logical(validMask(:));

    if numel(valid) ~= M
        error('PE_func:ValidMask', 'ValidMask must have the same length as the time series.');
    end

    valid = valid & ~isnan(x);

end

% A new segment starts at every valid sample preceded by an invalid one.
segStart = valid & [true; ~valid(1:end-1)];
segId    = cumsum(segStart);
segId(~valid) = 0;
nSeg     = sum(segStart);

end


%   >>>>>   BANDT-POMPE PDF (segment-aware)   <<<<<
function [P, nUsed] = PDF_segments(x, segId, D, tau, direction, index_map, powers, N)
% Uses every ordinal vector (x_j, x_{j+tau}, ..., x_{j+(D-1)tau}), which is
% exactly the set of windows produced by the previous loop over the tau
% phases k = 1:tau, but keeps only vectors whose first and last samples
% belong to the same segment. Since segments are contiguous, this
% guarantees that every sample in between is also valid.

% Reverse the time series when a negative direction is requested.
if direction < 0
    x     = flipud(x);
    segId = flipud(segId);
end

M    = numel(x);
span = (D-1) * tau;
nWin = M - span;

P = zeros(N, 1);
nUsed = 0;

if nWin <= 0
    P = P ./ sum(P);                     % NaN, as in the original version
    return;
end

j = (1:nWin).';

% Keep only windows fully contained in a single segment.
keep = segId(j) > 0 & segId(j) == segId(j + span);
j    = j(keep);
nUsed = numel(j);

if nUsed == 0
    warning('PE_func:NoWindows', 'No ordinal vector fits inside a single segment for tau = %d.', tau);
    P = P ./ sum(P);                     % NaN
    return;
end

% Build all embedding windows at once.
idx = bsxfun(@plus, j, (0:D-1) * tau);
X   = x(idx);
if nUsed == 1
    X = reshape(X, 1, D);
end

% Sort each window to obtain its ordinal pattern.
[~, order] = sort(X, 2, 'ascend');

% Convert each ordinal pattern into a numeric symbol.
symbols = order * transpose(powers);

% Convert numeric symbols into Lehmer / Bandt-Pompe indices.
lehmer_indices = index_map(symbols);

% Count the occurrences of each ordinal pattern.
P = accumarray(lehmer_indices, 1, [N, 1]);

% Normalize the histogram to obtain a probability distribution.
P = P ./ sum(P);

end


%   >>>>>   SHANNON ENTROPY   <<<<<
function [S] = Shannon(P)

S = 0;

for i = 1:length(P)

    % Ignore zero-probability states to avoid log2(0).
    if P(i) == 0
        continue;
    end

    S = S - P(i) .* log2(P(i));

end

end


%   >>>>>   JENSEN-SHANNON COMPLEXITY   <<<<<
function C = Complexity_Jensen_Shannon(P, N, H)

% Uniform probability distribution.
P_e = (1/N) * ones(N, 1);

% Delta distribution used for normalization.
P_m = zeros(N, 1);
P_m(1) = 1;

% Jensen-Shannon divergence between two probability distributions.
Div_Jensen_Shannon = @(P1, P2) ...
    Shannon((P1 + P2)/2) - Shannon(P1)/2 - Shannon(P2)/2;

% Normalization constant.
Q0 = Div_Jensen_Shannon(P_m, P_e)^(-1);

% Normalized disequilibrium.
Q = Q0 * Div_Jensen_Shannon(P, P_e);

% Statistical complexity.
C = Q * H;

end


%   >>>>>   FISHER INFORMATION   <<<<<
function [S] = Fisher(P)

N = length(P);
S = 0;

% Compute the discrete Fisher information.
for i = 1:N-1
    S = S + (sqrt(P(i+1)) - sqrt(P(i)))^2;
end

% Use a different normalization factor for extreme delta distributions.
if (P(1) == 1 && sum(P(2:N)) == 0) || ...
   (P(N) == 1 && sum(P(1:N-1)) == 0)
    F0 = 1;
else
    F0 = 0.5;
end

S = F0 * S;

end