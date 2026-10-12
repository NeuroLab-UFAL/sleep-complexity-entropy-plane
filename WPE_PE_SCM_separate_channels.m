% Computes the entropy (H), statistical complexity (C) and Fisher
% information (F) of every iEEG channel, separately for each brain state
% and for each embedding delay tau.
%
% Input
%   MatlabFile_Filtered.mat, containing one matrix per brain state
%   (data_W, data_N2, data_N3, data_R). Each matrix has one column per
%   channel and one row per sample (200 Hz). Recorded segments are
%   separated by zero-amplitude buffers, and the series are zero-padded
%   to a fixed length. The zeros must NOT be removed beforehand, because
%   WPE_func and PE_func use them to find the segment boundaries.
%
% Output
%   One file per method and brain state, named
%   HC_separate_channels_<method>_<state>.mat, containing
%     results_all_channels  [length(tauList) x 3 x n_channels], where the
%                           second dimension is [H, C, F]
%     n_segments            [n_channels x 1] continuous segments per channel
%     n_valid_samples       [n_channels x 1] recorded samples per channel
%     n_windows_used        [length(tauList) x n_channels] ordinal vectors used
%     n_windows_excluded    [length(tauList) x n_channels] ordinal vectors
%                           discarded for crossing a segment boundary
%     tauList, D            embedding delays and embedding dimension

tic
clc
clear

file = load('MatlabFile_Filtered.mat');

D = 6;                    % embedding dimension
tauList = 1:20;           % embedding delays, in samples (5 to 100 ms at 200 Hz)
n_channels = 1012;        % channels with data in all four brain states
fields = ["data_W", "data_N2", "data_N3", "data_R"];   % Wake, N2, N3, REM
methods = ["WPE", "PE"];  % weighted (WPE) and conventional (PE) permutation entropy;
                          % PE is the control analysis of Section 4.2

% Every series is truncated to its first N_VALID samples of recorded
% activity, so that all channels and states have the same length.
N_VALID = 10800;          % 54 s at 200 Hz

% A run of at least MIN_ZERO_RUN consecutive zeros is treated as a buffer
% (10 samples = 50 ms at 200 Hz). This must be the same criterion used
% inside WPE_func and PE_func, whose default is also 10.
MIN_ZERO_RUN = 10;

for method = methods
    for field = fields

        data = file.(field);
        assert(size(data, 2) == n_channels, ...
            'Unexpected number of channels in %s', field);

        % Preallocate the outputs for this method and brain state.
        results_all_channels = zeros(length(tauList), 3, n_channels);
        n_segments           = zeros(n_channels, 1);
        n_valid_samples      = zeros(n_channels, 1);
        n_windows_used       = zeros(length(tauList), n_channels);
        n_windows_excluded   = zeros(length(tauList), n_channels);

        for i = 1:n_channels

            % Time series WITH the zero buffers still in place.
            time_series = data(:, i);

            % Locate the recorded samples and cut the series right after
            % the N_VALID-th one. The cut is made on the raw series, so the
            % buffers inside the kept portion are preserved.
            valid = buffer_mask(time_series, MIN_ZERO_RUN);
            last_idx = find(cumsum(valid) == N_VALID, 1);
            assert(~isempty(last_idx), ...
                'Channel %d (%s) has fewer than %d valid samples', ...
                i, field, N_VALID);
            time_series = time_series(1:last_idx);

            % Compute [H, C, F] for every tau, using only ordinal vectors
            % that lie entirely inside one continuous segment.
            if method == "WPE"
                [results_all_channels(:, :, i), info] = ...
                    WPE_func(time_series, tauList, D);
            else
                [results_all_channels(:, :, i), info] = ...
                    PE_func(time_series, tauList, D);
            end

            % Keep the segment and window counts for later quality checks.
            n_segments(i)            = info.nSegments;
            n_valid_samples(i)       = info.nValidSamples;
            n_windows_used(:, i)     = info.nWindowsUsed;
            n_windows_excluded(:, i) = info.nWindowsExcluded;

        end

        % Summary of the segmentation for this method and brain state.
        fprintf('%s %s: segments per channel = %d to %d (median %g); valid samples = %d to %d\n', ...
            method, field, min(n_segments), max(n_segments), median(n_segments), ...
            min(n_valid_samples), max(n_valid_samples));

        filename = sprintf('HC_separate_channels_%s_%s.mat', method, field);
        save(filename, "results_all_channels", "n_segments", "n_valid_samples", ...
            "n_windows_used", "n_windows_excluded", "tauList", "D");

    end
end

toc

%% Local function

function valid = buffer_mask(x, minZeroRun)
% Returns a logical vector, true where the sample belongs to recorded
% activity and false for NaN samples and for runs of at least minZeroRun
% consecutive zeros (the buffers between recorded segments).

valid = ~isnan(x);

% Find the start and end index of every run of exact zeros.
is_zero   = (x == 0);
d         = diff([false; is_zero; false]);
run_start = find(d == 1);
run_end   = find(d == -1) - 1;

% Mark as invalid only the runs that are long enough to be buffers.
% Shorter runs of zeros are kept as data.
for r = find(run_end - run_start + 1 >= minZeroRun).'
    valid(run_start(r):run_end(r)) = false;
end

end