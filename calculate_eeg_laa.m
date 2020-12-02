function [output_metrics, single_values] =...
    calculate_eeg_laa( alpha_powers, channel_locs )
%   This function calculates Lateral Alpha Asymmetries across electrodes
% and lobes; Frontal, Central, Perietal, Occipital, and Temporal (if
% applicable). The function calculates the asymmetry according to the
% formula:
% 
%                    Right Electrode/Region Alpha Power
%       LAA = log  (-------------------------------------)
%                10   Left Electrode/Region Alpha Power
% 
% 
%   Usage:
%       output_metrics = calculate_eeg_laa( alpha_powers, channel_locs )
%       [output_metrics, single_values] = calculate_eeg_laa( ... )
% 
%   Inputs:
%       
%       alpha_powers      : P x N double-array of alpha-powers, collected
%                           from an EEG system with P channels/electrodes
%                           and N time series windowed values.
%       
%       channel_locs      : P x 1 array of structures as stored in EEGLAB
%                           variable 'EEG' as chanlocs, or derived from a
%                           '.ced' or '.loc'  locations file.
% 
%   Outputs:
%       
%       output_metrics    : structure of lateral alpha asymmetry values for
%                           each window for each lateral electrode-pair in
%                           the time-series data along with the regional
%                           lateral asymmetry values as a mean power over
%                           the given region.
%       
%       single_values     : structure similar to output_metrics with only a
%                           single value for each variable.
% 
output_metrics = struct();
if nargin == 1, error('Required arguments missing.');
end, eloc = struct2table(channel_locs);

% Mapping Lobic Values
map = nan(size(eloc.labels, 1), 2); lobe = [];
central_lobe = []; central_map = [];

for i_el = 1:size(map, 1)
    if ~isempty(find(map(:) == i_el, 1))
        continue;
        
    elseif contains(eloc.labels{i_el}, 'z', 'IgnoreCase', true)
        i_val1 = split(eloc.labels{i_el}, {'z', 'Z'});
        central_map = [central_map, i_el];
        central_lobe = [central_lobe; {i_val1{1}}];
        
    else
        i_val1 = split(eloc.labels{i_el},...
            {'1','2','3','4','5','6','7','8','9','0'});
        i_val1 = i_val1{1};
        
        i_val2 = split(eloc.labels{i_el}, i_val1);
        i_val2 = str2double(i_val2{2});
        
        if rem(i_val2, 2), clear i_*; continue; end
        i_val3 = find(strcmpi(eloc.labels, sprintf('%s%d',...
            i_val1, i_val2 + - 1)));
        
        map(i_el,:) = [i_el, i_val3];
        lobe = [lobe; {i_val1}];
        
    end
end, map(isnan(sum(map, 2)), :) = [];

% Cleaning-up the Lobic Region Values
for i_el = 1:length(lobe)
    lobe{i_el} = upper(lobe{i_el});
    if strcmpi(lobe{i_el}, 'fp'), lobe{i_el} = 'Fp'; end
end, clear i_*;

% Calculations
if ~isempty(map)
    
    % Electrode Alpha Asymmetries
    for i = 1:length(lobe)
        output_metrics.(sprintf('%s_%s',...
            eloc.labels{map(i,1)}, eloc.labels{map(i,2)}))...
            = alpha_powers(map(i,1), :) ./ alpha_powers(map(i,2), :);
    end, clear i;
    
    % Lobe and Pairing Lengths
    vals_sz1 = numel(central_lobe);
    vals_sz2 = numel(lobe);
    
    % Frontal Values
    vals_f = [...
        map(find(strcmp(lobe, 'F'), vals_sz2), :);...
        map(find(strcmp(lobe, 'AF'), vals_sz2), :);...
        map(find(strcmp(lobe, 'Fp'), vals_sz2), :)];
    if ~isempty(vals_f), output_metrics.frontal =...
        rms(alpha_powers(vals_f(:, 1), :), 1) ./...
        rms(alpha_powers(vals_f(:, 2), :), 1);
    end
    
    % Central Values
    vals_f = map(find(contains(lobe, 'C'), vals_sz2), :);
    if ~isempty(vals_f), output_metrics.central =...
        rms(alpha_powers(vals_f(:, 1), :), 1) ./...
        rms(alpha_powers(vals_f(:, 2), :), 1);
    end
    
    % Parietal Values
    vals_f = map(find(contains(lobe, 'P'), vals_sz2), :);
    if ~isempty(vals_f), output_metrics.parietal =...
        rms(alpha_powers(vals_f(:, 1), :), 1) ./...
        rms(alpha_powers(vals_f(:, 2), :), 1);
    end
    
    % Temporal Values
    vals_f = map(find(contains(lobe, 'T'), vals_sz2), :);
    if ~isempty(vals_f), output_metrics.temporal =...
        rms(alpha_powers(vals_f(:, 1), :), 1) ./...
        rms(alpha_powers(vals_f(:, 2), :), 1);
    end
    
    % Occipital Values
    vals_f = map(find(contains(lobe, 'O'), vals_sz2), :);
    if ~isempty(vals_f), output_metrics.occipital =...
        rms(alpha_powers(vals_f(:, 1), :), 1) ./...
        rms(alpha_powers(vals_f(:, 2), :), 1);
    end
    
end, clear vals_* i_* *map *lobe;

% Calculate Single-Values
fvars = fieldnames(output_metrics);
single_values = struct();

for i = 1:length(fvars)
    single_values.(fvars{i}) =  log10(mean(output_metrics.(fvars{i}), 2));
    output_metrics.(fvars{i}) =  log10(output_metrics.(fvars{i}));
end, clear i fvars;

end