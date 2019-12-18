classdef geyesub < spatial_filter
% Gerneralized eye subspace subtraction algorithm according to 
% Parra et al. 2005 and Kobler et al. 2017
%
% Parra, L. C., Spence, C. D., Gerson, A. D., & Sajda, P. (2005). 
%   Recipes for the linear analysis of EEG. NeuroImage, 28(2), 326–341. 
%   https://doi.org/10.1016/j.neuroimage.2005.05.032
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

    %%
    properties
        
        W = NaN; % unmixing matrix
        A = NaN; % mixing matrix
        V = NaN; % unmixing matrix (generalized least squares)
        R = NaN; % one step correction matrix
        eeg_chan_idxs = NaN;
        
        % penalized logistic regression parameters
        plr_lambda = 1e2; % regularization parameter
        plr_tol = 1e-6;
        plr_maxiter = 1e3;
        
        umix_lambda = 1e-3; % unmixing matrix regularization parameter

    end
    %%
    methods
        %%
        function fit(obj, data, labels, eeg_chan_idxs)
            % find subspaces for horizontal, vertical and blink components
            % with penalized logistic regression, combine the subspaces
            % and find their mixing matrices, and finally find a new
            % unmixing matrix that consideres the covariance during resting
            % periods.
            
            obj.eeg_chan_idxs = eeg_chan_idxs;
            
            X = data(eeg_chan_idxs,:);
            y = labels(:,:);
            
            % equalize the number of observations for each artifact type
            X_rst = X(:, y == 6);
            y_rst = y(:, y == 6);

            X_art = X(:, y > 0 & y < 6);
            y_art = y(:, y > 0 & y < 6);

            [X_art, y_art] = equalizeLabels(X_art, y_art, 'max', false);

            X = cat(2, X_art, X_rst);
            y = cat(2, y_art, y_rst);
            
            % split the data
            X_h = X(:, y == 1 | y == 2);
            y_h = y(y == 1 | y == 2);
            X_v = X(:, y == 3 | y == 4);
            y_v = y(y == 3 | y == 4);
            X_bu = X(:, y == 3 | y == 5);
            y_bu = y(y == 3 | y == 5);
            X_bd = X(:, y == 4 | y == 5);
            y_bd = y(y == 4 | y == 5);       
            X_rst  = X(:,y == 6 );
            
            % fit the individual unmixing vectors
            wh = obj.pen_logreg(X_h, y_h);
            wv = obj.pen_logreg(X_v, y_v);
            wbu = obj.pen_logreg(X_bu, y_bu);
            wbd = obj.pen_logreg(X_bd, y_bd);
            
            obj.W = [wh'; wv'; wbu'; wbd'];
            
            % find the individual mixing vectors
            a_H = obj.wtoa_shrink(wh, X_h);
            a_V = obj.wtoa_shrink(wv, X_v);
            a_B = obj.wtoa_shrink(wbu, X_bu);
            a_B2 = obj.wtoa_shrink(wbd, X_bd);
            
            obj.A=[a_H, a_V, a_B, a_B2];
            
            % recompute the unmixing matrix, considering the noise
            % covarinace matrix
            Lambda = obj.umix_lambda * eye(4);

            R_n = cov_shrink(X_rst');   
            R_n_inv = R_n^-1;

            obj.V = (obj.A'*R_n_inv*obj.A + Lambda)\obj.A'*R_n_inv;            
            
            % compute one step reconstruction matrix
            obj.R = eye(length(eeg_chan_idxs)) - obj.A*obj.V;
        end
        %%
        function data = apply(obj, data)
            % perform the correction step
            data(obj.eeg_chan_idxs,:) = obj.R * data(obj.eeg_chan_idxs, :);
        end
    
        %%
        function [w_plr] = pen_logreg(obj, X, labels)
            % penalized logistic regression according to the algorithm
            % outlined in Parra et al. 2005
            
            [ndim, nobs] = size(X);
            classes = unique(labels);
            % only allow applicable for a two-class problem
            assert(length(classes) == 2);
            % check for balanced set
            assert(sum(labels == classes(1)) == sum(labels == classes(2)));
            
            tmp = labels;
            labels(tmp==classes(1)) = 1;
            labels(tmp==classes(2)) = 0;
            clear tmp
            
            Lambda = obj.plr_lambda * eye(ndim);
            
            f = @(x) 1./(1 + exp(-x));
            
            w_plr = zeros(ndim, 1);
            dw_norm = 1;
            iter = 1;
            fprintf('PLR working')
            while dw_norm > obj.plr_tol
                
                if iter > obj.plr_maxiter || dw_norm > 1e6
                    warning('partial least squares did not converge!');
                    w_plr = nan(ndim,1);
                    break
                end
                
                if mod(iter, 10) == 0
                    fprintf('dw=%e', dw_norm);
                end
                p = f(w_plr' * X);
                P = sparse(1:nobs, 1:nobs, p .* (1-p));
                g = X * (labels - p)' - Lambda * w_plr;
                H = X * P * X' + Lambda;
                
                dw = H\g;
                w_plr = w_plr + dw;
                
                dw_norm = norm(dw);
                
                fprintf('.');
                iter = iter +1;
            end
            fprintf('DONE.\n')
            
        end
        
        %%
        function [A, Y] = wtoa_shrink(obj, W, X)
            % find the mixing matrix associated to the weight matrix and
            % singals in X
            
            % compute the estimated source signal Y
            Y = W' * X;
            
            % compute covariance matrices of X and Y
            nc_x = size(X,1);
            nc_y = size(Y,1);
            
            Z = cat(1, X, Y);
            [Czz, ~] = cov_shrink(Z');
            
            Cxx = Czz(1:nc_x, 1:nc_x);
            Cyy = Czz((1:nc_y)+nc_x, (1:nc_y)+nc_x);
            
            % compute the activation patterns
            A = Cxx * W / Cyy;
            
        end
        
    end
    
end



