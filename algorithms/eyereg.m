classdef eyereg < spatial_filter
%eyereg eye artifact regression algorithm according to Schlögl et el. 2007
%
% Schlögl, A., Keinrath, C., Zimmermann, D., Scherer, R., Leeb, R., 
%   & Pfurtscheller, G. (2007). 
%   A fully automated correction method of EOG artifacts in EEG recordings. 
%   Clinical Neurophysiology, 118(1), 98–104.
%
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
    
    %% class properties
    properties
        W = NaN; % the weight matrix
        eeg_chan_idxs = NaN; % eeg channel indices within the signal
        eog_chan_idxs = NaN; % eog channel indices within the signal
    end

    %% methods
    methods
        %%
        function fit(obj, data, labels, eeg_chan_idxs, eog_chan_idxs)
            % fits the regression coefficients for eye movement and blink
            % periods
            
            % regress the EEG channel activity (Y) form the EOG channel
            % activity (Y)
            X_art  = data(eog_chan_idxs, labels < 6 );
            Y_art  = data(eeg_chan_idxs, labels < 6 );
            
            Z = cat(1, X_art, Y_art);
            nc_x = size(X_art,1);
            nc_y = size(Y_art,1);
            
            Czz = cov_shrink(Z');
            
            Cxx = Czz(1:nc_x, 1:nc_x);
            Cyx = Czz((1:nc_y)+nc_x, 1:nc_x);
            
            % compute minimum mean squared error weight matrix W
            obj.W = Cyx / Cxx;
            obj.eeg_chan_idxs = eeg_chan_idxs;
            obj.eog_chan_idxs = eog_chan_idxs;
        end
        
        %%
        function data = apply(obj, data)
            % applies eog artifact correction to new data samples
            data(obj.eeg_chan_idxs,:) = data(obj.eeg_chan_idxs, :) ...
                - obj.W * data(obj.eog_chan_idxs, :);
        end
    end
end

