% Detects saccades as periods during which the EOG signal is outside 
% the interval spanned by +/- threshold
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
function [ thresh_signal ] = eyepos2events( EEG, channel_idx, label, threshold, min_length, t_extend )
%eyepos2events Detects saccades in the EEG dataset using the channel 
%   defined by channel_idx for trials with the label 'label' that are
%   outside the interval spanned by +/- threshold

    if nargin < 5
        min_length = 0;
    else
        min_length = min_length * EEG.srate;
    end

    if nargin < 6
        t_extend = 0;
    else
        t_extend = t_extend * EEG.srate;
    end

    thresh_signal = cat(1, EEG.data(channel_idx,:,:) > threshold, ...
        EEG.data(channel_idx,:,:) < -threshold);

    thresh_signal(:,:,EEG.etc.trial_labels ~= label) = false;

    if min_length > 0

        thresh_signal = min_duration_threshold(thresh_signal, min_length);
    end

    if t_extend > 0

        b = ones(round(t_extend),1)'/t_extend;

        thresh_signal(1,:,:) = filtfilt(b, 1, squeeze(double(thresh_signal(1,:,:))));
        thresh_signal(2,:,:) = filtfilt(b, 1, squeeze(double(thresh_signal(2,:,:))));

        thresh_signal = thresh_signal > 0;
    end

end

