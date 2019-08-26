% Eye movement and blink correction reference paradigm as presented in
% Kobler et al. 2017.
%
% The implementation is using PsychToolbox <version 3.0.14>
% (Brainard, 1997; Pelli, 1997; Kleiner et al, 2007). 
% You find the latest version on the project website
% http://psychtoolbox.org/
%
% Kobler, R. J., Sburlea, A. I., & Müller-Putz, G. R. (2017). 
%   A comparison of ocular artifact removal methods for block design based 
%   electroencephalography experiments. 
%   In Proceedings of the 7th Graz Brain-Computer Interface Conference 
%   pp. 236–241.
%   https://doi.org/10.3217/978-3-85125-533-1-44
%
% Brainard, D. H. (1997) .
%   The Psychophysics Toolbox, Spatial Vision 10:433-436.
%
% Pelli, D. G. (1997) 
%   The VideoToolbox software for visual psychophysics: Transforming 
%   numbers into movies, Spatial Vision 10:437-442.
%
% Kleiner M., Brainard D., Pelli D. (2007).
%   "What’s new in Psychtoolbox-3?" Perception 36 ECVP Abstract Supplement.
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

clearvars
close all
clc
sca % if this fails, you probably miss PsychToolbox on your matlab path

%% set global parameters

% refreshrate of the paradigm
loop_fs = 60; % Hz

photodiode = false;
photo_dimension = 120;
photo_dur = 0.15;

%% set the codes and trials per condition
start_trial_code = 1000;
end_trial_code = -1000;

rest_code = 1; n_rest = 9;
horz_code = 2; n_horz = 6;
vert_code = 3; n_vert = 6;
blink_code = 4; n_blink = 6;

%% define the timings of each trial

rest_dur = 1;
track_dur = 10;
t_rise = 2;
min_break = 2;
max_break = 3;


%% instantiate the library
if exist('lsl_loadlib', 'file')
    disp('Loading the library...');
    liblsl = lsl_loadlib();
    
    % create streams
    disp('Creating a marker stream...');
    info_marker = lsl_streaminfo(liblsl,'eyeblock-marker','marker',...
        1,0,'cf_float32','myuniquesourceid1355678');

    disp('Opening a marker outlet...');
    outlet_marker = lsl_outlet(info_marker,0, 1000);

    disp('Creating a xy position stream...');
    info_stimuli = lsl_streaminfo(liblsl,'eyeblock-stimuli','stimuli',...
        3,loop_fs,'cf_float32','myuniquesourceid1599874');

    % append some meta-data

    channels = info_stimuli.desc().append_child('channels');

    channels.append_child('channel') ...
        .append_child_value('label', 'target.X') ...
        .append_child_value('type', 'Scale');
    channels.append_child('channel') ...
        .append_child_value('label', 'target.Y') ...
        .append_child_value('type', 'Scale');    
    channels.append_child('channel') ...
        .append_child_value('label', 'target.S') ...
        .append_child_value('type', 'Scale');

    disp('Opening a xy position outlet...');
    outlet_stimuli = lsl_outlet(info_stimuli,0, 1000);    

    pause(0.5)
    
    logging = true;
else
    warning('WARNING: LSL not configured. Logging of markers and stimuli disabled!\n');
    fprintf('Press any key to confirm and proceed anyway...\n')
    logging = false;
    pause
end

%% generate the sequence of tirals

n_trials = n_rest + n_horz + n_vert + n_blink;

sequence = [rest_code*ones(1,n_rest) horz_code*ones(1,n_horz) ...
            vert_code*ones(1,n_vert) blink_code*ones(1,n_blink)];

n_rep_max = 2;
        
while true
    indexes = randperm(n_trials); % compute random permutation

    % check for long (> n_rep_max) sequences of the same condition
    cand_sequence = sprintf('%d', abs(diff(sequence(indexes))));
    rep_idxs = strfind(diff(cand_sequence), repmat('0', 1, n_rep_max+1));

    % accept the sequence if there are not too long single condition
    % sequences
    if isempty(rep_idxs)
        sequence = sequence(indexes);
        break
    end
end

%% open the psychotoolbox screen
AssertOpenGL;
% PsychDebugWindowConfiguration();

screens = Screen('Screens');
screenNumber = max(screens); % figure will be displayed in the second screen
[w, rect] = Screen('OpenWindow', screenNumber, 0); % get window number w and coordinates of the screen's rectangle
Screen('TextSize', w, 30); % setting text size for future text

white = WhiteIndex(w);
Priority(MaxPriority(w));

%% generate signals

% frequency moving dot in vertical and horizontal eye movements
f_signal = 0.5; %Hz

% generate signals for horz and vert trials
t_trans = (0:loop_fs-1)/loop_fs;
win_trans = (1 - cos(2 *pi * f_signal * t_trans)) * 0.5;

window = cat(2, win_trans, ones(1,(track_dur - 1 - 1)*loop_fs), 1 - win_trans);

t_task = (0:track_dur*loop_fs-1)/loop_fs;

eyesignal =  sin(2 * pi * f_signal * t_task) .* window;

% generate signal for blink trials
% blink_trigger = ((1.0 - bartlett(floor(0.25 * loop_fs))) + 0.5) / 1.5;
t_trans = linspace(0,0.25,loop_fs*0.25);
blink_trigger = cos(2*pi * 4 * t_trans) / 4 + 0.75;

blink_pause = ones(1, floor(0.75 * loop_fs));
blink_window = cat(2,blink_trigger, blink_pause);

blinksignal = repmat(blink_window, 1, ceil(track_dur*loop_fs / length(blink_window)));

break_durations = rand(n_trials)*(max_break - min_break) + min_break;

%% set dot size, colors etc.

dotsize = 40;
dotsize_inner = dotsize/5;
color_active = [2 122 198]';
color_inner = [0 0 0]';

% range of the trajectories and their central point
max_w = rect(3);
max_h = rect(4);

[x_center, y_center] = RectCenter(rect);

x_range = floor(min(rect(3:4))/2) - dotsize; y_range = x_range;

oval_pos = [x_center-dotsize/2 y_center-dotsize/2 ...
    x_center+dotsize/2 y_center+dotsize/2];

oval_pos_inner = [x_center-dotsize_inner/2 y_center-dotsize_inner/2 ...
    x_center+dotsize_inner/2 y_center+dotsize_inner/2];


%% pause to start the labrecorder
if logging
    fprintf('\n------------------------------\n')
    fprintf('NOTICE: streams open. Start recording and press any key to proceed...')
    fprintf('\n------------------------------\n')
    pause
end

%% start of paradigm

% Do initial flip
vbl = Screen('Flip', w);

% pause for 5 seconds before the paradigm starts
pause(5)

for k_trial = 1:n_trials
    
    % experimenter logging
    fprintf('trial %02d/%02d: cond: %d\n', k_trial, n_trials, sequence(k_trial));
    
    if logging % send marker at the beginning of the trial
        outlet_marker.push_sample(start_trial_code); 
    end
    
    %% 'rest' part of the trial
    t_start = GetSecs();
    t_rel = 0;
    while t_rel < rest_dur
        
        t_iter = GetSecs();
        Screen('FillOval', w, color_active, oval_pos);
        Screen('FillOval', w, color_inner, oval_pos_inner);

        if photodiode && t_rel < photo_dur % add the photodiode
            Screen('FillRect', w, white,[max_w-photo_dimension 0 max_w photo_dimension]);
        end

        Screen('Flip',w); % update image
        
        if logging
            outlet_stimuli.push_sample([x_center; y_center; dotsize])
        end
        
        t_op = GetSecs() - t_iter;
        t_sleep = max(1/loop_fs - t_op, 0);
        t_rel = WaitSecs(t_sleep) - t_start;
    end

    %% 'tracking' part of the trial

    if logging
        outlet_marker.push_sample(sequence(k_trial));
    end
    
    t_start = GetSecs();
    t_rel = 0; loop_iter = 0;
    while t_rel < track_dur
        
        t_iter = GetSecs();
        loop_iter = loop_iter + 1;

        switch sequence(k_trial)

            case rest_code
                
                Screen('FillOval', w , color_active, oval_pos);                   
                Screen('FillOval', w, color_inner, oval_pos_inner);
                if logging
                    outlet_stimuli.push_sample([x_center; y_center; dotsize])
                end
                Screen('Flip',w); % update image

            case {horz_code,vert_code}

                if sequence(k_trial) == horz_code
                    dx = eyesignal(loop_iter) * x_range;
                    Screen('FillOval', w , color_active, ...
                        [oval_pos(1)+dx oval_pos(2) oval_pos(3)+dx oval_pos(4)]);
                    Screen('FillOval', w, color_inner, ...
                        [oval_pos_inner(1)+dx oval_pos_inner(2) oval_pos_inner(3)+dx oval_pos_inner(4)]);
                    if logging
                        outlet_stimuli.push_sample([x_center + dx; y_center; dotsize])
                    end
                else
                    dy = eyesignal(loop_iter) * y_range;
                    Screen('FillOval', w , color_active, ...
                        [oval_pos(1) oval_pos(2)+dy oval_pos(3) oval_pos(4)+dy]);
                    Screen('FillOval', w, color_inner, ...
                        [oval_pos_inner(1) oval_pos_inner(2)+dy oval_pos_inner(3) oval_pos_inner(4)+dy]);

                    if logging
                        outlet_stimuli.push_sample([x_center; y_center + dy; dotsize])
                    end
                end
                Screen('Flip',w); % update image

            case blink_code

                delta = blinksignal(loop_iter)*dotsize/2;
                delta_inner = blinksignal(loop_iter)*dotsize/5/2;
                
                Screen('FillOval', w , color_active, ...
                    [x_center-delta y_center-delta x_center+delta y_center+delta]);
                Screen('FillOval', w, color_inner, oval_pos_inner);

                if logging
                    outlet_stimuli.push_sample([x_center; y_center; delta*2])
                end
                Screen('Flip',w); % update image
        end
        
        t_op = GetSecs() - t_iter;
        t_sleep = max(1/loop_fs - t_op, 0);
        t_rel = WaitSecs(t_sleep) - t_start;
    end
    
    if logging % push marker sample at the end of the trial
        outlet_marker.push_sample(end_trial_code); 
    end
    
    %% 'break' part of the trial
    
    % clear the screen
    Screen('Flip',w);
    Screen('Flip',w);
    
    % continue logging the outlet stimuli
    t_start = GetSecs();
    t_rel = 0;
    while t_rel < break_durations(k_trial)
        t_iter = GetSecs();
        
        if logging
            outlet_stimuli.push_sample([x_center; y_center; dotsize])
        end
        
        t_op = GetSecs() - t_iter;
        t_sleep = max(1/loop_fs - t_op, 0);
        t_rel = WaitSecs(t_sleep) - t_start;
    end
    
end
%% tidy up

% clear the screen
sca

% clear the lsl variables
if logging
    pause(0.5)
    clear lib

    pause(0.5);

    clear info_marker info_stimuli
    clear outlet_marker
    clear outlet_stimuli
end