% Eye artifact correction data preprocessing demonstration script
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

st = dbstack;
if length(st) < 2

    clear all
    close all
    clc

    addpath('./utility');
    addpath('./external');
    addpath('./algorithms');    

    file_name = 'demo_trainset';
end


%% load and pre-process the train data

% the dataset was recorded at a rate of 1kHz
% offline, a 50 Hz notch filter was applied before the data was resampled
% at 200 Hz
EEG = pop_loadset([file_name '.set'], '.');

%% apply a 0.4 Hz zero-phase high-pass filter to the eeg and eog signals

d = fdesign.highpass('N,F3db', 2, 0.4, EEG.srate);
hHp = design(d, 'butter', 'SOSScaleNorm', 'linf');    
SOS = hHp.sosMatrix;
G = hHp.ScaleValues;

if exist('eeg_chantype', 'file') % eeglab version 14 and before
    eeg_chan_idxs = eeg_chantype(EEG, 'EEG');
else % eeglab version 2019
    eeg_chan_idxs = eeg_decodechan(EEG.chanlocs, 'EEG', 'type');
end

EEG.data(eeg_chan_idxs,:) = filtfilt(SOS, G, EEG.data(eeg_chan_idxs,:)')';

%% compute the EOG derivatives

leog_chan_idxs = eeg_chaninds(EEG, {'EOGL1', 'EOGL2', 'EOGL3'});
reog_chan_idxs = eeg_chaninds(EEG, {'EOGR1', 'EOGR2', 'EOGR3'});

eog_chan_idxs = EEG.nbchan + (1:3);

% compute the horizontal eog derivative
EEG.data(eog_chan_idxs(1),:) = EEG.data(reog_chan_idxs(1),:) - ...
    EEG.data(leog_chan_idxs(1),:);
% compute the vertical eog derivative as the average of the left and right
% sides
EEG.data(eog_chan_idxs(2),:) = mean(EEG.data([reog_chan_idxs(2), leog_chan_idxs(2)],:) - ...
    EEG.data([reog_chan_idxs(3), leog_chan_idxs(3)],:),1);
% compute the radial eog derivative as the average activity of the eog
% channels
EEG.data(eog_chan_idxs(3),:) = mean(EEG.data([reog_chan_idxs, leog_chan_idxs],:),1);

[EEG.chanlocs(eog_chan_idxs).labels] = deal('HEOG', 'VEOG', 'REOG');
[EEG.chanlocs(eog_chan_idxs).type] = deal('EOGLPF');
EEG.nbchan = EEG.nbchan + 3;

%% apply a low-pass filter to the eog derivative signals for saccade
% and blink detection

d = fdesign.lowpass('N,F3dB', 2, 5, EEG.srate);
h = design(d, 'butter');

SOS = h.sosMatrix;
G = h.ScaleValues;

eoglpf_chan_idxs = EEG.nbchan + (1:length(eog_chan_idxs));

EEG.data(eoglpf_chan_idxs,:) = filtfilt(SOS, G, EEG.data(eog_chan_idxs,:)')';

[EEG.chanlocs(eoglpf_chan_idxs).labels] = deal('HEOG_lpf', 'VEOG_lpf', 'REOG_lpf');
[EEG.chanlocs(eoglpf_chan_idxs).type] = deal('EOGLPF');
EEG.nbchan = EEG.nbchan + 3;

%% epoch the dataset into trials
events_to_epoch = {'1', '2', '3', '4'}; % 1=rest,2=horz,3=vert,4=blink
epoch_slice_times = [1 9]; % omit the first and last second

% adjust the trial events so that they are maintained in the epoched data
[~, epoch_evt_idxs] = pop_selectevent(EEG, 'type',events_to_epoch);

trial_labels = str2double({EEG.event(epoch_evt_idxs).type});

[EEG, accepted_trials] = pop_epoch( EEG, events_to_epoch, epoch_slice_times);
EEG.etc.trial_labels = trial_labels(accepted_trials);

EEG = eeg_checkset(EEG, 'eventconsistency');


% visual inspection and artifact rejection
if exist([file_name '_badepochs.mat'], 'file')
    load([file_name '_badepochs.mat'], 'rej_mask');
    disp(['loading bad epochs from file!. bad epochs: ' num2str(find(rej_mask))]);
else

    tmpreject = trial2eegplot(rej_mask, zeros(EEG.nbchan, EEG.trials), EEG.pnts, [1, 1, 0.783]);

    eegplot( EEG.data, 'srate', EEG.srate, 'title', 'manual artifact rejection', ...
              'limits', [EEG.xmin EEG.xmax]*1000, 'command', ' ', ...
              'ploteventdur', 'off', 'winlength', 1, 'spacing', 75, ...
              'events', EEG.event, 'eloc_file', EEG.chanlocs, ...
              'submean', 'off','winrej', tmpreject);    

    uiwait(gcf)
    
    if ~isempty(TMPREJ)
        rej_mask = logical(eegplot2trial(TMPREJ, EEG.pnts, EEG.trials));
    else
        rej_mask = false(1, EEG.trials);
    end

    save([file_name '_badepochs.mat'], 'rej_mask');
end

rej_mask = logical(rej_mask);
disp(['rejecting epochs: ' num2str(find(rej_mask))]);
if sum(rej_mask) > 0
    
    bad_evts = EEG.etc.trial_labels(rej_mask);
    edges = unique(EEG.etc.trial_labels);
    cnts = histc(bad_evts, edges );
    cnts_all = histc(EEG.etc.trial_labels, edges );
    fprintf('evt %d: %d of %d\n', cat(1, edges(:)', cnts(:)', cnts_all(:)'));
    
    EEG = pop_rejepoch( EEG, rej_mask,0);
    EEG.etc.trial_labels(rej_mask) = [];
end

%% detect the saccades and blinks
EEG = util_detect_saccades_and_blinks(EEG);

mask_blink_idx = eeg_chaninds(EEG, 'eye-blink');
mask_up_idx = eeg_chaninds(EEG, 'eye-u');
mask_down_idx = eeg_chaninds(EEG, 'eye-d');
mask_left_idx = eeg_chaninds(EEG, 'eye-l');
mask_right_idx = eeg_chaninds(EEG, 'eye-r');
mask_rest_idx = eeg_chaninds(EEG, 'eye-fix');

% encode the differnet artifact classes
EEG.etc.eye.labels = {'right', 'left', 'up', 'down', 'blink', 'rest'};
EEG.etc.eye.codes = 1:length(EEG.etc.eye.labels);

EEG.nbchan = EEG.nbchan + 1;
label_chan_idx = EEG.nbchan;

EEG.data(label_chan_idx,:,:) = EEG.data(mask_right_idx,:,:) + 2 * EEG.data(mask_left_idx,:,:) ...
             + 3 * EEG.data(mask_up_idx,:,:) + 4 * EEG.data(mask_down_idx,:,:) ...
             + 5 * EEG.data(mask_blink_idx,:,:) + 6 * EEG.data(mask_rest_idx,:,:);
EEG.chanlocs(label_chan_idx).labels = 'artifactclasses';
EEG.chanlocs(label_chan_idx).type = 'LABEL';


