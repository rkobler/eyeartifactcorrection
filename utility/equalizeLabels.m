% Equalizes the samples of differenty labels within the dataset
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
function [X_out, l_out] = equalizeLabels(X, l, type, do_shuffle)

    if nargin < 3
        type = 'max';
    end
    
    if nargin < 4
        do_shuffle = false;
    end    
    
    if nargin < 5
        except = [];
    end
    
    labels = unique(l);
    
    n_labels = length(labels);
    
    n_l = zeros(n_labels, 1);
    
    for idx = 1:n_labels
        n_l(idx) = sum(l == labels(idx));
    end
    
    if strcmp(type, 'max')
        n_max = max(n_l);
    else
        n_max = min(n_l);
    end
    
    
    l_out = [];
    X_out = [];
    for idx = 1:n_labels
        l_idxs = find(l == labels(idx));
        if do_shuffle
            s_idxs = randperm(n_l(idx));
            l_idxs = l_idxs(s_idxs);
        end
        idxs = 1+mod(0:n_max-1, n_l(idx));
        l_idxs = l_idxs(idxs);
        
        X_out = cat(2, X_out, X(:, l_idxs));
        l_out = cat(2, l_out, labels(idx) * ones(1,n_max));
    end

end