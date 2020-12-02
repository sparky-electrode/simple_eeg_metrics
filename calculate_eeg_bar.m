function [output_metrics, single_values] =...
    calculate_eeg_bar( beta_powers, alpha_powers, channel_locs )
%   This function calculates Beta-Alpha Ratio (BAR) across electrodes and
% the lobes; Frontal, Central, Perietal, Occipital, and Temporal (if
% applicable). The function calculates the ratio according to the given
% formula:
% 
%                              Beta Band Power
%                      BAR = -------------------.
%                              Alpha Band Power
% 
% 
%   Usage:
%       output_metrics = calculate_eeg_bar( beta_powers, alpha_powers, channel_locs )
%       [output_metrics, single_values] = calculate_eeg_bar( ... )
% 
%   Inputs:
%       
%       beta_powers       : P x N double-array of beta-powers, collected
%                           from an EEG system with P channels/electrodes
%                           and N time series windowed values. The matrix
%                           dimensions need to be the same as that for
%                           alpha_powers.
%       
%       alpha_powers      : P x N double-array of alpha-powers, collected
%                           from an EEG system with P channels/electrodes
%                           and N time series windowed values. The matrix
%                           dimensions need to be the same as that for
%                           beta_powers.
%       
%       channel_locs      : P x 1 array of structures as stored in EEGLAB
%                           variable 'EEG' as chanlocs, or derived from a
%                           '.ced' or '.loc'  locations file.
% 
%   Outputs:
%       
%       output_metrics    : structure of beta-alpha ratio values for each
%                           window for each lateral electrode-pair in the
%                           time-series data along with the regional
%                           ratios as a mean power over the given region.
%       
%       single_values     : structure similar to output_metrics with only a
%                           single value for each variable.
% 
output_metrics = struct();
if nargin == 2, error('Required arguments missing.');
end, eloc = struct2table(channel_locs);

% Beta Alpha Ratio
output_metrics.electrode = beta_powers ./ alpha_powers;

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
    
    % Lobe and Pairing Lengths
    vals_sz1 = numel(central_lobe);
    vals_sz2 = numel(lobe);
    
    % Frontal Values
    vals_f1 = [...
        map(find(strcmp(lobe, 'F'), vals_sz2), :);...
        map(find(strcmp(lobe, 'AF'), vals_sz2), :);...
        map(find(strcmp(lobe, 'Fp'), vals_sz2), :)];
    vals_f2 = [...
        central_map(find(strcmp(central_lobe, 'F'), vals_sz1)),...
        central_map(find(strcmp(central_lobe, 'AF'), vals_sz1)),...
        central_map(find(strcmp(central_lobe, 'Fp'), vals_sz1))];
    vals_f = [vals_f1(:)', vals_f2(:)'];
    
    if ~isempty(vals_f), output_metrics.frontal =...
            rms(beta_powers(vals_f, :), 1) ./...
            rms(alpha_powers(vals_f, :), 1);
    end
    
    % Central Values
    vals_f1 = map(find(contains(lobe, 'C'), vals_sz2), :);
    vals_f2 = central_map(find(strcmp(central_lobe, 'C'), vals_sz1));
    vals_f = [vals_f1(:)', vals_f2(:)'];
    
    if ~isempty(vals_f), output_metrics.central =...
            rms(beta_powers(vals_f, :), 1) ./...
            rms(alpha_powers(vals_f, :), 1);
    end
    
    % Parietal Values
    vals_f1 = map(find(contains(lobe, 'P'), vals_sz2), :);
    vals_f2 = central_map(find(strcmp(central_lobe, 'P'), vals_sz1));
    vals_f = [vals_f1(:)', vals_f2(:)'];
    
    if ~isempty(vals_f), output_metrics.parietal =...
            rms(beta_powers(vals_f, :), 1) ./...
            rms(alpha_powers(vals_f, :), 1);
    end
    
    % Temporal Values
    vals_f1 = map(find(contains(lobe, 'T'), vals_sz2), :);
    vals_f2 = central_map(find(strcmp(central_lobe, 'T'), vals_sz1));
    vals_f = [vals_f1(:)', vals_f2(:)'];
    
    if ~isempty(vals_f), output_metrics.temporal =...
            rms(beta_powers(vals_f, :), 1) ./...
            rms(alpha_powers(vals_f, :), 1);
    end
    
    % Occipital Values
    vals_f1 = map(find(contains(lobe, 'O'), vals_sz2), :);
    vals_f2 = central_map(find(strcmp(central_lobe, 'O'), vals_sz1));
    vals_f = [vals_f1(:)', vals_f2(:)'];
    
    if ~isempty(vals_f), output_metrics.occipital =...
            rms(beta_powers(vals_f, :), 1) ./...
            rms(alpha_powers(vals_f, :), 1);
    end
    
end, clear vals_* *lobe *map;

% Calculate Single-Values
fvars = fieldnames(output_metrics);
for i_el = 1:length(fvars), single_values.(fvars{i_el})...
        = mean(output_metrics.(fvars{i_el}), 2);
    % fprintf('Finished Calculating Values for ''%s''\n', fvars{i_el});
end, clear i_*;

end