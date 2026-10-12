% Averages the iEEG channels of each brain state into a single time series
% and computes its entropy (H), statistical complexity (C) and Fisher
% information (F), for each embedding delay tau.
%
% The set of channels contributing to the average changes whenever a
% channel enters or leaves a zero-amplitude buffer. A new segment is
% started at every such change, so that no ordinal vector spans samples
% whose average was computed over different sets of channels.
%
% Input
%   MatlabFile_Filtered.mat, containing one matrix per brain state
%   (data_W, data_N2, data_N3, data_R). Each matrix has one column per
%   channel and one row per sample (200 Hz). Recorded segments are
%   separated by zero-amplitude buffers, and the series are zero-padded
%   to a fixed length.
%
% Output
%   One file per method and brain state, named
%   HC_mean_series_<method>_<state>.mat, containing
%     HCF      [length(tauList) x 3], where the second dimension is [H, C, F]
%     info     segment and window counts returned by WPE_func / PE_func
%     tauList  embedding delays
%     D        embedding dimension

tic

All_files = load('MatlabFile_Filtered.mat');

D = 6;                    % embedding dimension
tauList = 1:100;          % embedding delays, in samples (5 to 500 ms at 200 Hz)
fields = ["data_W", "data_N2", "data_N3", "data_R"];   % Wake, N2, N3, REM
method_list = ["WPE", "PE"];   % weighted and conventional permutation entropy

% Delay used only for the console report of excluded ordinal vectors.
% It must be one of the values in tauList (tau = 8 is 40 ms at 200 Hz).
tau_report = 8;
assert(any(tauList == tau_report), 'tau_report must be in tauList');

% A run of at least min_zero_run consecutive zeros is treated as a buffer
% (10 samples = 50 ms at 200 Hz). This is also the default used inside
% WPE_func and PE_func, and the two criteria must match.
min_zero_run = 10;

% The mean series is truncated to its first N_VALID samples in which at
% least one channel has recorded activity (54 s at 200 Hz).
N_VALID = 10800;

for field = fields

    data = All_files.(field);

    % Mark the recorded samples of each channel and replace the buffers
    % by NaN, so that they are ignored when averaging.
    valid = buffer_mask(data, min_zero_run);   % true = recorded activity
    data(~valid) = NaN;

    % Average across channels, using only the channels with recorded
    % activity at each time point.
    mean_series = mean(data, 2, 'omitnan');

    % Start a new segment whenever the set of channels contributing to the
    % average changes. The first sample always starts a segment.
    set_changed = [true; any(valid(2:end, :) ~= valid(1:end-1, :), 2)];
    seg_id = cumsum(set_changed);

    % Time points without any valid channel belong to no segment (label 0).
    seg_id(~any(valid, 2)) = 0;

    % Truncate right after the N_VALID-th sample with recorded activity.
    last_idx = find(cumsum(seg_id > 0) == N_VALID, 1);
    assert(~isempty(last_idx), ...
        'Mean series of %s has fewer than %d valid samples', field, N_VALID);
    mean_series = mean_series(1:last_idx);
    seg_id      = seg_id(1:last_idx);

    for method = method_list

        % Compute [H, C, F] for every tau, using only ordinal vectors that
        % lie entirely inside one segment of the mean series.
        if method == "WPE"
            [HCF, info] = WPE_func(mean_series, tauList, D, 'SegmentId', seg_id);
        else
            [HCF, info] = PE_func(mean_series, tauList, D, 'SegmentId', seg_id);
        end

        filename = sprintf('HC_mean_series_%s_%s.mat', method, field);
        save(filename, 'HCF', 'info', 'tauList', 'D');

        % Report how many ordinal vectors were discarded for crossing a
        % segment boundary at the reference delay.
        fprintf('%s %s: %d segments; vectors excluded at tau = %d: %d\n', ...
            method, field, info.nSegments, tau_report, ...
            info.nWindowsExcluded(tauList == tau_report));

    end
end

toc

%% Local function

function valid = buffer_mask(data, minZeroRun)
% Returns a logical matrix with the same size as data, true where the
% sample belongs to recorded activity. NaN samples and runs of at least
% minZeroRun consecutive zeros (the buffers between recorded segments)
% are marked false, independently for each channel (column).

valid = ~isnan(data);

for c = 1:size(data, 2)

    % Find the start and end index of every run of exact zeros.
    is_zero   = (data(:, c) == 0);
    d         = diff([false; is_zero; false]);
    run_start = find(d == 1);
    run_end   = find(d == -1) - 1;

    % Mark as invalid only the runs that are long enough to be buffers.
    % Shorter runs of zeros are kept as data.
    for r = find(run_end - run_start + 1 >= minZeroRun).'
        valid(run_start(r):run_end(r), c) = false;
    end

end

end