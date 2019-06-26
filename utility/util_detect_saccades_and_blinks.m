% Detects saccades and blinks in EOG derivative signals (muV)
% for an EEG dataset containing signals according to the paradigm 
% described by Kobler et al. 2017
%
% Kobler, R. J., Sburlea, A. I., & Müller-Putz, G. R. (2017). 
%   A comparison of ocular artifact removal methods for block design based 
%   electroencephalography experiments. 
%   In Proceedings of the 7th Graz Brain-Computer Interface Conference 
%   pp. 236–241.
%   https://doi.org/10.3217/978-3-85125-533-1-44
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
function EEG = util_detect_saccades_and_blinks(EEG)
%% detect blinks and saccades in trials

% thresholds were optimized for active electrodes with
% with Brainamp (Brain Products GmbH) and g.USBamp (g.tec Gmbh) amplifiers
blink_fp_th = 100;
blink_vert_fp_th = 150;
blink_tp_th = 75;

saccade_fp_th = 100;
saccade_tp_th = 10;


eye_blink_sig = false(1,EEG.pnts,EEG.trials);

eye_blink_sig = eye_blink_sig | eogblink2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), blink_fp_th, 0.75, 1);
eye_blink_sig = eye_blink_sig | eogblink2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), blink_fp_th, 0.75, 2);
eye_blink_sig = eye_blink_sig | eogblink2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), blink_vert_fp_th, 0.75, 3);
eye_blink_sig = eye_blink_sig | eogblink2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), blink_tp_th, 0.025, 4);

eye_lr_sig = false(2,EEG.pnts,EEG.trials);

eye_lr_sig = eye_lr_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'HEOG_lpf'), 1, saccade_fp_th, 0.2, 0.05);
eye_lr_sig = eye_lr_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'HEOG_lpf'), 2, saccade_tp_th, 0.2);
eye_lr_sig = eye_lr_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'HEOG_lpf'), 3, saccade_fp_th, 0.2);

eye_l_sig = eye_lr_sig(1,:,:);
eye_r_sig = eye_lr_sig(2,:,:);

eye_ud_sig = false(2,EEG.pnts,EEG.trials);

eye_ud_sig = eye_ud_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), 1, saccade_fp_th, 0.2, 0.05);
eye_ud_sig = eye_ud_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), 2, saccade_fp_th, 0.2);
eye_ud_sig = eye_ud_sig | eyepos2events(EEG, eeg_chaninds(EEG, 'VEOG_lpf'), 3, saccade_tp_th, 0.2);

eye_u_sig = eye_ud_sig(1,:,:);
eye_d_sig = eye_ud_sig(2,:,:);

eye_fix_sig = repmat(permute(EEG.etc.trial_labels == 1, [1 3 2]), 1, EEG.pnts,1);

art_sig = zeros(size(eye_fix_sig));
label_sig = zeros(size(eye_fix_sig));

%% create artifact signal 

for i_trial = 1:EEG.trials
   
    t_label = EEG.etc.trial_labels(i_trial);
    
    label_sig(1,:, i_trial) = t_label;
    
    switch(t_label)
        case 1
            art_sig(1, :, i_trial) = eye_blink_sig(1, :, i_trial) ...
                                    | eye_l_sig(1, :, i_trial) ...
                                    | eye_r_sig(1, :, i_trial) ...
                                    | eye_d_sig(1, :, i_trial) ...
                                    | eye_u_sig(1, :, i_trial);
            eye_fix_sig(1, :, i_trial) = eye_fix_sig(1, :, i_trial) & ~art_sig(1, :, i_trial);
            eye_fix_sig(1, :, i_trial) = min_duration_threshold(eye_fix_sig(1, :, i_trial), 0.25*EEG.srate);
            art_sig(1, :, i_trial) = ~eye_fix_sig(1, :, i_trial);
        case 2
            art_sig(1, :, i_trial) = eye_blink_sig(1, :, i_trial) ...
                                    | eye_d_sig(1, :, i_trial) ...
                                    | eye_u_sig(1, :, i_trial);
            eye_l_sig(1, :, i_trial) = eye_l_sig(1, :, i_trial) & ~art_sig(1, :, i_trial);
            eye_l_sig(1, :, i_trial) = min_duration_threshold(eye_l_sig(1, :, i_trial), 0.1*EEG.srate);
            eye_r_sig(1, :, i_trial) = eye_r_sig(1, :, i_trial) & ~art_sig(1, :, i_trial);
            eye_r_sig(1, :, i_trial) = min_duration_threshold(eye_r_sig(1, :, i_trial), 0.1*EEG.srate);
            
        case 3
            art_sig(1, :, i_trial) = eye_blink_sig(1, :, i_trial) ...
                                    | eye_l_sig(1, :, i_trial) ...
                                    | eye_r_sig(1, :, i_trial);            
            eye_d_sig(1, :, i_trial) = eye_d_sig(1, :, i_trial) & ~art_sig(1, :, i_trial);
            eye_d_sig(1, :, i_trial) = min_duration_threshold(eye_d_sig(1, :, i_trial), 0.1*EEG.srate);
            eye_u_sig(1, :, i_trial) = eye_u_sig(1, :, i_trial) & ~art_sig(1, :, i_trial);
            eye_u_sig(1, :, i_trial) = min_duration_threshold(eye_u_sig(1, :, i_trial), 0.1*EEG.srate);

        case 4
            eye_blink_sig(1, :, i_trial) = min_duration_threshold(eye_blink_sig(1, :, i_trial), 0.1*EEG.srate);
    end
end

%% add masks as channels

% set masks to zero for non-relevant trials
eye_fix_sig(:, :, EEG.etc.trial_labels ~= 1) = false;
eye_r_sig(:, :, EEG.etc.trial_labels ~= 2) = false;
eye_l_sig(:, :, EEG.etc.trial_labels ~= 2) = false;
eye_u_sig(:, :, EEG.etc.trial_labels ~= 3) = false;
eye_d_sig(:, :, EEG.etc.trial_labels ~= 3) = false;
eye_blink_sig(:, :, EEG.etc.trial_labels ~= 4) = false;

% add the channels to the dataset
EEG.data = cat(1, EEG.data, eye_l_sig, eye_r_sig, eye_d_sig, eye_u_sig, eye_blink_sig, eye_fix_sig, art_sig, label_sig);
EEG.chanlocs(end+1).labels = 'eye-l';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-r';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-d';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-u';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-blink';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-fix';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'eye-art';
EEG.chanlocs(end).type = 'ARTIFACT';
EEG.chanlocs(end+1).labels = 'label';
EEG.chanlocs(end).type = 'STATE';
EEG.nbchan = EEG.nbchan + 8;
