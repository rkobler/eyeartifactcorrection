classdef eyeeeg < spatial_filter
%eyeeeg eye artifact correction algorithm according to Plöchl et el. 2012
%
% Plöchl, M., Ossandón, J. P., & König, P. (2012). 
%   Combining EEG and eye tracking: identification, characterization, 
%   and correction of eye movement artifacts in 
%   electroencephalographic data. 
%   Frontiers in Human Neuroscience, 6(October), 1–23.
%   https://doi.org/10.3389/fnhum.2012.00278
%
% Copyright (C) 2019 Reinmar Kobler, Graz University of Technology, Austria
% <reinmar.kobler@tugraz.at>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
    
    %%
    properties

        R = NaN; % one step eye artifact correction matrix
        eeg_chan_idxs = NaN;
        bad_idxs = NaN; % indices of eye artifact releated components
        
        th_ratio = 1.1; % threshold ratio
        
    end
    %%
    methods
        %%
        function fit(obj, data, labels, eeg_chan_idxs, icaw, icawi)
            % projects the data to IC space using icaw, detects eye 
            % artifact related ICs based on thresholding

            obj.eeg_chan_idxs = eeg_chan_idxs;
            
            % find saccade and fixation events
            saccade = labels(:,:) > 0 & labels(:,:) < 6;
            fixation = labels(:,:) == 6;
            
            % transform the IC space
            Z = icaw * data(eeg_chan_idxs,:);
            
            % compute variances during saccade and fixation periodes for
            % all ICs
            var_sacc = var(Z(:,saccade), 0, 2);
            var_fix = var(Z(:,fixation), 0, 2);

            % thresholding
            var_ratio = var_sacc ./ var_fix;
            
            obj.bad_idxs = find(var_ratio > obj.th_ratio);
            
            % compute one step reconstruction matrix
            obj.R = eye(length(eeg_chan_idxs)) -  icawi(:, obj.bad_idxs) * icaw(obj.bad_idxs,:);
            
        end
        
        function data = apply(obj, data)
            % perform the correction step
            data(obj.eeg_chan_idxs,:) = obj.R * data(obj.eeg_chan_idxs, :);
        end
    end 
end



