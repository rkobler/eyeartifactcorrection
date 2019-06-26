% Clears ones that have a shorter duration than a minimum length
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
function [ mask ] = min_duration_threshold( mask, min_length )
%MIN_DURATION_THRESHOLD Clears sequenes of ones in mask that 
%   are shorter than min_length

    [n_comp, ~, n_trials] = size(mask);

    for epoch = 1:n_trials
        
        for comp = 1:n_comp
            
            if length(unique(mask(comp,:,epoch))) > 1
                [~, locs, w] = findpeaks(double(mask(comp,:,epoch)), 'MinPeakWidth',min_length);

                idxs = arrayfun(@(x,y) x:y, locs, locs + w, 'UniformOutput', false);
                idxs = cell2mat(idxs);
                
                mask(comp,:,epoch) = 0;
                mask(comp,idxs,epoch) = 1;

            end
        end

    end
end

