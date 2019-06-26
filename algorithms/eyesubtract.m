classdef eyesubtract < spatial_filter
% Reference eye subspace subtraction algorithm according to Parra et al. 2005
%
% The code is based on the reference EEGLAB plugin
% 'eyesubtract v1.0'
% The plugin is available online
%   https://sccn.ucsd.edu/wiki/EEGLAB_Extensions 
% Authors: Xiang Zhou (zhxapple@hotmail.com, 2005)
%          with Adam Gerson (adg71@columbia.edu, 2005),
%          and Lucas Parra (parra@ccny.cuny.edu, 2005),
%          and Paul Sajda (ps629@columbia,edu 2005)
%
% Parra, L. C., Spence, C. D., Gerson, A. D., & Sajda, P. (2005). 
%   Recipes for the linear analysis of EEG. NeuroImage, 28(2), 326–341. 
%   https://doi.org/10.1016/j.neuroimage.2005.05.032
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
        
        A = NaN; % mixing matrix
        V = NaN; % unmixing matrix
        R = NaN; % one step eye artifact correction matrix
        eeg_chan_idxs = NaN;
        
    end
    %%
    methods
        %%
        function fit(obj, data, labels, eeg_chan_idxs)
            % find subspaces for horizontal, vertical and blink components
            % and remove them from the original data

            data_blink = data(eeg_chan_idxs, labels == 5);
            data_left  = data(eeg_chan_idxs, labels == 2);
            data_right = data(eeg_chan_idxs, labels == 1);
            data_up    = data(eeg_chan_idxs, labels == 3);
            data_down  = data(eeg_chan_idxs, labels == 4);

            % Forward model

            a_H  =  obj.difference(data_left, data_right);
            a_V  =  obj.difference(data_up, data_down );
            a_B  =  obj.Maximumpower(data_blink);

            obj.eeg_chan_idxs = eeg_chan_idxs;
            obj.A=[a_H a_V a_B];
            obj.V = (obj.A'*obj.A)\obj.A';
        
            % precompute one step correction matrix
            obj.R = eye(length(eeg_chan_idxs)) - obj.A * obj.V;
        end
        
        function data = apply(obj, data)
            % performs the correction step
            data(obj.eeg_chan_idxs,:) = obj.R * data(obj.eeg_chan_idxs, :);
        end
    end
    %% original functions from eyesutract eeglab plugin v1.0
    
    % Copyright (C) 2005 Xiang Zhou, Adam Gerson, Lucas Parra and Paul Sajda
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program; if not, write to the Free Software
    % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    methods (Access = protected)
        
    function Difvec=difference(obj, EEG_1,EEG_2)

        %Difvec   = mean(EEG_1-EEG_2,2);
        Difvec = mean(EEG_1,2)-mean(EEG_2,2);

        Difvec   = Difvec./norm(Difvec);

    end

    function Max_eigvec=Maximumpower(obj, EEGdata)

        [Channel,SRate,Epoch]=size(EEGdata);

        for e=1:Epoch
            for i=1:Channel
                EEGdata(i,:,e)=EEGdata(i,:,e)-mean(EEGdata(i,:,e));
            end;
        end;


        % Maximum Power method

        a=EEGdata(:,:);

        pow_a=a*a';

        [vb,tmp] = eig(pow_a);
        Max_eigvec=vb(:,end)./norm(vb(:,end));


        end
    
    
    end
    
end



