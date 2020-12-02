function [output_metrics, single_values] =...
    calculate_eeg_metrics( bandpow_input, channel_locs )
%   This function calculates the Lateral Alpha-Asymmetry, the Theta-Beta
% Ratio (TBR), and the Beta-Alpha Ratio (BAR) across electrodes and
% the lobes; Frontal, Central, Perietal, Occipital, and Temporal (if
% applicable). The function is a wrapper on the functions
% calcAlphaAsymmetry(...), calcTBR(...), and calcBAR(...).
% 
%   Usage:
%       output_metrics = calculate_eeg_metrics( bandpow_input, channel_locs )
%       [output_metrics, single_values] = calculate_eeg_metrics( ... )
% 
%   Inputs:
%       
%       bandpow_input     : Structure of double-values as calculated by the
%                           function calcBandPow(...), with each power as a
%                           P x N double-array for P electrodes and N
%                           windows in the time-series data.
%       
%       channel_locs      : P x 1 array of structures as stored in EEGLAB
%                           variable 'EEG' as chanlocs, or derived from a
%                           '.ced' or '.loc'  locations file.
% 
%   Outputs:
%       
%       output_metrics    : structure of metric values based on the band
%                           powers as defined in the input structure for
%                           each lateral electrode-pair in the time-series
%                           data along with the regional metrics as a mean
%                           power over the given region.
%       
%       single_values     : structure similar to output_metrics with only a
%                           single value for each variable.
% 
if nargin < 2, error('Required Arguments Missing.');
elseif isempty(channel_locs), error('Channel Locations cannot be Empty.');
elseif istable(channel_locs), eloc = table2struct(channel_locs);
else, eloc = channel_locs;
end

% Initialize
output_metrics = bandpow_input;
single_values = struct();

% Fill-in the Single Powers
fvars = fieldnames(bandpow_input);
for i = 1:length(fvars)
    single_values.(fvars{i}) = mean(bandpow_input.(fvars{i}), 2);
end

% Alpha Asymmetry
if isfield(bandpow_input, 'alpha')
    [met, sval] = calcAlphaAsymmetry(bandpow_input.alpha, eloc);
    fvars = fieldnames(met);
    for i = 1:length(fvars)
        output_metrics.(sprintf('%s_alp_asym', fvars{i})) = met.(fvars{i});
        single_values.(sprintf('%s_alp_asym', fvars{i})) = sval.(fvars{i});
    end
end

% Theta-Beta Ratio (TBR)
if isfield(bandpow_input, {'theta','beta'})
    [met, sval] = calcTBR(bandpow_input.theta, bandpow_input.beta, eloc);
    fvars = fieldnames(met);
    for i = 1:length(fvars)
        output_metrics.(sprintf('%s_tbr', fvars{i})) = met.(fvars{i});
        single_values.(sprintf('%s_tbr', fvars{i})) = sval.(fvars{i});
    end
end

% Beta-Alpha Ratio (BAR)
if isfield(bandpow_input, {'alpha','beta'})
    [met, sval] = calcBAR(bandpow_input.beta, bandpow_input.alpha, eloc);
    fvars = fieldnames(met);
    for i = 1:length(fvars)
        output_metrics.(sprintf('%s_bar', fvars{i})) = met.(fvars{i});
        single_values.(sprintf('%s_bar', fvars{i})) = sval.(fvars{i});
    end
end

clear met sval i fvars;
end