Data_reordered = permute (tongue, [3,2,1]);
SISTool
eeg_data = data.data; % [22 × 750 × n_trials], e.g., 22×750×72
seg_pts_cell = data.segpnts; % n_trials×1 cell array
%% --- LOAD ELECTRODE LABELS AND SPHERICAL COORDINATES
% Define channel labels
labels = {'Fz', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'C5', 'C3', 'C1', 'Cz', ...
          'C2', 'C4', 'C6', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'P1', 'Pz', 'P2', 'POz'};
disp('Defined electrode labels for 22 channels.');

% Define spherical coordinates (theta, phi) in degrees
theta = [0.0, 45.0, 15.0, 0.0, 345.0, 315.0, 90.0, 45.0, 15.0, 0.0, 345.0, 315.0, 270.0, ...
         45.0, 15.0, 0.0, 345.0, 315.0, 15.0, 0.0, 345.0, 0.0];
phi = [90.0, 67.5, 67.5, 67.5, 67.5, 67.5, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, ...
       22.5, 22.5, 22.5, 22.5, 22.5, 0.0, 0.0, 0.0, -22.5];
disp('Defined spherical coordinates (theta, phi) for 22 channels.');
>> % Verify the lengths match
if length(labels) ~= length(theta) || length(labels) ~= length(phi)
    error('Number of labels, theta, and phi must match (should be 22).');
end
>> % Number of electrodes
n_channels = length(labels); % Should be 22
% Verify the toolbox is available
if isempty(which('GetGH'))
    error('CSDtoolbox not found. Please ensure GetGH.m is in the path.');
end

% Create montage struct with field name "lab" as expected by GetGH
M = struct('lab', {labels}, 'theta', theta, 'phi', phi);

% Compute G and H matrices for CSD using GetGH
try
    [G, H] = GetGH(M); % Compute G and H matrices using the montage struct
    lambda = 1e-5; % Default smoothing constant
    head = 10; % Default head radius in cm
    disp('Successfully computed G and H matrices for CSD.');
catch e
    warning('Failed to compute G and H matrices: %s. CSD will be skipped.', e.message);
    G = []; H = [];
end
% Number of trials
[~, ~, n_trials] = size(eeg_data);
connectivity_graphs = cell(0, 1); % Flat list of 22×22 matrices
trial_indices = cell(0, 1); % Metadata: trial index for each matrix
segment_indices = cell(0, 1); % Metadata: segment index for each matrix
actual_segment_counts = zeros(n_trials, 1);

% Sampling rate for conversion (250 Hz)
sampling_rate = 250; % Hz
ms_to_samples = sampling_rate / 1000; % 0.25 samples/ms
>> % Initialize additional cell arrays for segment timing
segment_start_samples = cell(0, 1); % Start sample index for each graph
segment_end_samples = cell(0, 1);   % End sample index for each graph
segment_start_ms = cell(0, 1);      % Start time in milliseconds
segment_end_ms = cell(0, 1);        % End time in milliseconds
t_total = tic;
seg_times_ms = [];
for trial_idx = 1:n_trials
    trial_data = eeg_data(:, :, trial_idx); % 22×750
    
    % Convert seg_pts_cell from milliseconds to sample indices
    seg_trial = round(seg_pts_cell{trial_idx} * ms_to_samples) + 1;
    seg_trial(seg_trial < 1) = 1; % Clip lower bound
    seg_trial(seg_trial > 750) = 750; % Clip to trial length (750 samples)
    seg_trial = unique(seg_trial, 'stable'); % Remove duplicates
    n_segments_trial = numel(seg_trial) - 1;
    
    % Store the actual number of segments
    actual_segment_counts(trial_idx) = n_segments_trial;
    
    % Debug: Print segment points for key trials
    if ismember(trial_idx, [1, 4, 9])
        fprintf('Trial %d seg_pts (samples): %s\n', trial_idx, mat2str(seg_trial));
        fprintf('Trial %d segments: %d\n', trial_idx, n_segments_trial);
    end
    
    if n_segments_trial < 1
        fprintf('Trial %d: No segments\n', trial_idx);
        continue;
    end
    
    % Compute connectivity graphs for each segment
    for s = 1:n_segments_trial
        t_seg = tic;	
        idx_start = seg_trial(s);
        idx_end = seg_trial(s+1) - 1;
        if idx_end < idx_start
            idx_end = seg_trial(s+1);
        end
        segment = trial_data(:, idx_start:idx_end);
        if size(segment, 2) < 2
            fprintf('Trial %d, Segment %d: Skipped due to insufficient samples (%d)\n', ...
                    trial_idx, s, size(segment, 2));
            continue;
        end
        
        % Apply CSD transformation if G and H are available
        if ~isempty(G) && ~isempty(H)
            try
                [~, X] = CSD(segment, G, H, lambda, head); % X is the CSD-transformed data
            catch e
                fprintf('CSD failed for Trial %d, Segment %d: %s\n', trial_idx, s, e.message);
                X = segment; % Fallback: Use raw data
            end
        else
            X = segment; % Fallback: Use raw data
        end
        
        % Zero-mean each channel
        X = X - mean(X, 2);
        
        % Low-rank SVD reconstruction
        [U, S, V] = svd(X, 'econ');
        singular_values = diag(S);
        explained_variance = cumsum(singular_values.^2) / sum(singular_values.^2);
        k = find(explained_variance >= 0.95, 1, 'first');
        if isempty(k) || k < 1
            k = 1;
        end
        X_reconstructed = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
        
        % Compute Pearson correlation matrix
        C = corrcoef(X_reconstructed');
        if any(isnan(C(:)))
            fprintf('Trial %d, Segment %d: Skipped due to NaN in correlation\n', trial_idx, s);
            continue;
        end
        
        % Store the connectivity graph and timing information
        connectivity_graphs{end+1, 1} = C;
        trial_indices{end+1, 1} = trial_idx;
        segment_indices{end+1, 1} = s;
        segment_start_samples{end+1, 1} = idx_start;
        segment_end_samples{end+1, 1} = idx_end;
        segment_start_ms{end+1, 1} = idx_start * 1000 / sampling_rate;
        segment_end_ms{end+1, 1} = idx_end * 1000 / sampling_rate;
        seg_times_ms(end+1,1) = toc(t_seg) * 1000;
    end
end
total_ms = toc(t_total) * 1000;

fprintf('\n=== TIMING (Connectivity graph generation) ===\n');
fprintf('TOTAL time: %.2f ms (%.2f s)\n', total_ms, total_ms/1000);
if ~isempty(seg_times_ms)
    fprintf('Per-seg (ms): mean=%.3f | median=%.3f | p95=%.3f | max=%.3f | n=%d\n', ...
        mean(seg_times_ms), median(seg_times_ms), prctile(seg_times_ms,95), max(seg_times_ms), numel(seg_times_ms));
end
fprintf('=============================================\n\n');
% Save all data, including segment timing
if ~isempty(connectivity_graphs)
    save('connectivity_graphs_A02E_right.mat', 'connectivity_graphs', 'trial_indices', ...
         'segment_indices', 'labels', 'segment_start_samples', 'segment_end_samples', ...
         'segment_start_ms', 'segment_end_ms');
    fprintf('Saved %d connectivity graphs with timing information to connectivity_graphs_A01E_tongue.mat\n', ...
            length(connectivity_graphs));
else
    fprintf('No connectivity graphs available to save.\n');
end

