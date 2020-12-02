function bandpow_struct = calculate_eeg_bandpow(input_EEG, varargin)
%   The function 'calculate_eeg_bandpow' calculates the Band Powers for each row
% being a channel using the Spectrogram function. Arguments follow the same
% as that of the MATLAB function spectrogram.
% 
% Reference: <a href="https://www.mathworks.com/help/signal/ref/spectrogram.html">spectrogram</a>.

fprintf('Calculating bandpowers for %d Channels and %d Samples.\n',...
    size(input_EEG, 1), size(input_EEG, 2));

% Calculate Spectrogram
[s_out, freq_out] = my_multi_channel_spectrogram(input_EEG, varargin{:});

% Frequency Band Powers
f_var = {...
    'delta',        'theta',        'alpha',...
    'low_beta',     'high_beta',	'beta',         'gamma'};

% Frequency Band Values
freq_vals = [...
    00.5, 03.5;     04.0, 07.5;     08.0, 12.5;...
    12.5, 21.0;     21.0, 30.0;     12.5, 30.0;     30.0, 50.0	];

% Frequency Band Value Indices
f_val = nan(size(freq_vals));
for i = 1:numel(freq_vals)
    f_val(i) = find(freq_out >= freq_vals(i), 1);
end, clear i freq_*;

% Calculate and Store Band Powers
bandpow_struct = struct();
for i = 1:size(f_val, 1)
    if ~isnan(mean(f_val(i, :)))
        bandpow_struct.(f_var{i}) =...
            permute(rms(s_out(f_val(i,1):f_val(i,2), :, :)), [3,2,1]);
    end
end, clear i;
end

function [s_out, freq_out, t_out, ps_out, fc_out, tc_out]...
    = my_multi_channel_spectrogram(X, varargin)
%   <strong>Multi_channel_spectrogram</strong> is a wrapper function over MATLAB's function
% spectrogram, scripted to work with multi-channel data. The function works
% in the same manner as the original MATLAB function, except the output
% '<strong>s_out</strong>' is a 3D complex cell-matrix of size f frequencies x p samples x
% n channels. All the arguments used for the function 'spectrogram' remain
% unchanged for the wrapper function as well.
% 
% Check: <a href="https://www.mathworks.com/help/signal/ref/spectrogram.html">spectrogram</a>.
if nargin < 1 || nargin > 8, error('Inconsistent set of inputs.'); end

% Get the Basic Values other than 's_out'
[s, freq_out, t_out, ps_out, fc_out, tc_out] = spectrogram(X(1,:), varargin{:});

% Initialize the complex value (nan) matrix
s_out = nan(size(s,1), size(s,2), size(X, 1));

% Fill in the nan-values in 's_out'
for j = 1:17, s_out(:,:,j) = spectrogram(X(j,:), varargin{:}); end

% Clear unrequired variables
clear j s;
end