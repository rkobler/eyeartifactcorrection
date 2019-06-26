% Detects blinkds as periods during which the EOG signal is outside 
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
function [ blink_signal ] = eogblink2events( EEG, channel_idx, threshold, t_extend, label )
%eogblink2events Detects blinks in the EEG dataset using the channel 
%   defined by channel_idx for trials with the label 'label' that are
%   outside the interval spanned by +/- threshold

% all_event_mat = [];

blink_signal = squeeze(abs(EEG.data(channel_idx,:,:)) > threshold);

blink_signal(:,EEG.etc.trial_labels ~= label) = 0;

if t_extend > 0
    
    b = ones(round(t_extend*EEG.srate),1)';
    
    blink_signal = filtfilt(b, 1, double(blink_signal));
    
    blink_signal = blink_signal > 0;
end

blink_signal = reshape(blink_signal, 1, EEG.pnts, EEG.trials);

end

