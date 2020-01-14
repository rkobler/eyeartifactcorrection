% Eye artifact correction demonstration script
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

clear all
close all
clc

addpath('./utility');
addpath('./external');
addpath('./algorithms');  

% check for correct eeglab options

% startup eeglab
eeglab nogui

% IMPORTANT: use double precision in your eeglab settings! Otherwise some
% optimization algorithms might not converge.
% the next lines perform an automatic check
eeglab_options
if option_single
    error(['Wrong eeglab options: please switch off single precision ' ...
        'in the memory options! Execute the command ''pop_editoptions'' '...
        'to edit your configuration.']);
end
clearvars

%% load and preprocess the demonstration train and test set
file_name = 'demo_trainset';
demo_preprocessing
EEGTRN = EEG;

file_name = 'demo_testset';
demo_preprocessing

EEGTST = EEG;
clearvars -except EEGTRN EEGTST

%% load a dataset from the public repository on OSF
% 
% % download the dataset from https://osf.io/2qgrd/ and copy the files 
% % 'study02_p02_prep.set' and 'study02_p02_prep.fdt' to this folder
% 
% % load one recording
% dir = '.';
% file_name = 'study02_p02_prep.set';
% EEG = pop_loadset(file_name, dir);
% 
% trn_trial_idxs = find(EEG.etc.trial_blocks == 1);
% tst_trial_idxs = find(EEG.etc.trial_blocks == 2);
% 
% % partition the data into train and test sets based on on each trial's
% % block ID
% EEGTRN = pop_select(EEG, 'trial', find(EEG.etc.trial_blocks == 1));
% EEGTRN.etc.trial_labels = EEGTRN.etc.trial_labels(trn_trial_idxs);
% EEGTRN.etc.trial_ids = EEGTRN.etc.trial_ids(trn_trial_idxs);
% EEGTRN.etc.trial_blocks = EEGTRN.etc.trial_blocks(trn_trial_idxs);
% 
% EEGTST = pop_select(EEG, 'trial', find(EEG.etc.trial_blocks == 2));
% EEGTST.etc.trial_labels = EEGTST.etc.trial_labels(tst_trial_idxs);
% EEGTST.etc.trial_ids = EEGTST.etc.trial_ids(tst_trial_idxs);
% EEGTST.etc.trial_blocks = EEGTST.etc.trial_blocks(tst_trial_idxs);
% clearvars -except EEGTRN EEGTST

%% extract the features and fit the model to the training data

label_chan_idx = eeg_chaninds(EEGTRN, 'artifactclasses');
if exist('eeg_chantype', 'file') % eeglab version 14 and before
    eeg_chan_idxs = eeg_chantype(EEGTRN, 'EEG');
else % eeglab version 2019
    eeg_chan_idxs = eeg_decodechan(EEGTRN.chanlocs, 'EEG', 'type');
end

X_trn = EEGTRN.data;
y_trn = EEGTRN.data(label_chan_idx, :, :);

% fit the eye regression approach according to Schloegl et al. 2007
% algo = eyereg();
% eogder_chan_idxs = eeg_chaninds(EEGTRN, {'HEOG', 'VEOG'});
% algo.fit(X_trn, y_trn, eeg_chan_idxs, eogder_chan_idxs);

% fit the eye eeg algorithm according to Ploechl et al. 2012
% algo = eyeeeg();
% EEGTRN = pop_runica(EEGTRN, 'extended',1,'interupt','on', 'chanind', eeg_chan_idxs);
% algo.fit(X_trn, y_trn, eeg_chan_idxs, EEGTRN.icaweights * EEGTRN.icasphere, EEGTRN.icawinv);

% fit the eye subspace subtraction algorithm according to Parra et al. 2005
% algo = eyesubtract();
% algo.fit(X_trn, y_trn, eeg_chan_idxs);

% fit the generalized eye subspace subtraction algorithm according to 
% Kobler et al. 2017
% algo = geyesub();
% algo.fit(X_trn, y_trn, eeg_chan_idxs);

% fit the sparse generalized eye subspace subtraction algorithm according to 
% Kobler et al. 2019
algo = sgeyesub();
algo.fit(X_trn, y_trn, eeg_chan_idxs);


%% apply the model to the test dataset

data2 = algo.apply(EEGTST.data);

% plot the result at single sample level
pop_eegplot( EEGTST, 1, 1, 1, [], 'dispchans', 64, 'winlength', 1, ...
    'spacing', 50, 'data2', data2);

%% some evaluation plots

raweog_chan_idxs = eeg_chaninds(EEGTST, {'EOGL1', 'EOGL2', 'EOGL3', 'EOGR1', 'EOGR2', 'EOGR3'});
eeg_chan_idxs = setdiff(eeg_chan_idxs, raweog_chan_idxs);

eoglpf_chan_idxs = eeg_chaninds(EEGTST, {'HEOG_lpf', 'VEOG_lpf'});

eval_table = nan(length(eeg_chan_idxs), 4, 2);
metrics = {'REST (rmse)', 'HORZ (r)', 'VERT (r)', 'BLINK (r)'};
conditions = {'uncorrected', 'corrected'};

% compute the root mean squared error for rest trials
rest_samples = EEGTST.data(:,:, EEGTST.etc.trial_labels == 1);
rest_samples_c = data2(:,:, EEGTST.etc.trial_labels == 1);
eval_table(:,1,1) = sqrt(mean((rest_samples(eeg_chan_idxs,:) - rest_samples(eeg_chan_idxs,:)).^2,2));
eval_table(:,1,2) = sqrt(mean((rest_samples_c(eeg_chan_idxs,:) - rest_samples(eeg_chan_idxs,:)).^2,2));

% compute the correlations with eog derivatives for the other trials
horz_samples = EEGTST.data(:,:, EEGTST.etc.trial_labels == 2);
horz_samples_c = data2(:,:, EEGTST.etc.trial_labels == 2);

eval_table(:,2,1) = corr(horz_samples(eeg_chan_idxs,:)',horz_samples(eoglpf_chan_idxs(1),:)');
eval_table(:,2,2) = corr(horz_samples_c(eeg_chan_idxs,:)',horz_samples(eoglpf_chan_idxs(1),:)');

vert_samples = EEGTST.data(:,:, EEGTST.etc.trial_labels == 3);
vert_samples_c = data2(:,:, EEGTST.etc.trial_labels == 3);

eval_table(:,3,1) = corr(vert_samples(eeg_chan_idxs,:)',vert_samples(eoglpf_chan_idxs(2),:)');
eval_table(:,3,2) = corr(vert_samples_c(eeg_chan_idxs,:)',vert_samples(eoglpf_chan_idxs(2),:)');

blink_samples = EEGTST.data(:,:, EEGTST.etc.trial_labels == 4);
blink_samples_c = data2(:,:, EEGTST.etc.trial_labels == 4);

eval_table(:,4,1) = corr(blink_samples(eeg_chan_idxs,:)',blink_samples(eoglpf_chan_idxs(2),:)');
eval_table(:,4,2) = corr(blink_samples_c(eeg_chan_idxs,:)',blink_samples(eoglpf_chan_idxs(2),:)');

% plotting
figure
for i_plot = 1:size(eval_table,2)
    ax1 = subplot(2,size(eval_table,2),i_plot);
    topoplot(eval_table(:,i_plot,1), EEGTST.chanlocs(eeg_chan_idxs),'conv', 'on');
    title(metrics{i_plot})
    
    ax2 = subplot(2,size(eval_table,2),i_plot+size(eval_table,2));
    topoplot(eval_table(:,i_plot,2), EEGTST.chanlocs(eeg_chan_idxs),'conv', 'on');
    if i_plot == 1
        subplot(ax1)
        h_txt = text(-1,0.,conditions{1},'HorizontalAlignment','center', 'FontSize', 12);
        set(h_txt, 'rotation', 90)
        subplot(ax2)
        
        h_txt = text(-1,0.,conditions{2},'HorizontalAlignment','center', 'FontSize', 12);
        set(h_txt, 'rotation', 90)        
        
        set(ax1, 'clim', [0, 10])
        set(ax2, 'clim', [0, 10])
    else
        set(ax1, 'clim', [-1, 1])
        set(ax2, 'clim', [-1, 1])
    end
end