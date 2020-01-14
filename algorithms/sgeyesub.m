classdef sgeyesub < spatial_filter
% Sparse gerneralized eye subspace subtraction algorithm according to 
% Kobler et al. submitted
%
% Kobler, R. J., Sburlea, A. I., Lopes-Dias, C., Schwarz, A., Hirata M.
%   & Müller-Putz, G. R. 
%   Corneo-retinal-dipole and blink related eye artifacts can be corrected 
%   offline and online in electroencephalographic 
%   and magnetoencephalographic signals. 
%   NeuroImage, submitted
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
        C = NaN; % one step correction matrix
        eeg_chan_idxs = NaN;
        
        % elastic net penalized logistic regression parameters
        plr_lambda_l2 = 1e0;
        plr_lambda_l1 = 1e-2;
        plr_tol = 1e-3;
        plr_maxiter = 1e4;
        
        verbose = false;
        
    end
    %%
    methods
        %%
        function fit(obj, data, labels, eeg_chan_idxs)
            % remove subspaces in a two step procedure
            % first, remove the subspaces related to horizontal and
            % vertical eye movements
            % second, remove the subspace related to eye blinks
            
            obj.eeg_chan_idxs = eeg_chan_idxs;
            
            X = data(eeg_chan_idxs,:);
            y = labels(:,:);
            
            nbchan = length(eeg_chan_idxs);
            
            % estimat noise and eye movement artifact covariance matrix
            X_rst = X(:, y == 6);
            y_rst = y(:, y == 6);
            
            Rn = cov_shrink(X_rst');
            
            %% equalize the eye artifact data
            [X_eye, y_eye] = equalizeLabels(X(:, y >=1 & y <=4 ), y(:, y >=1 & y <=4 ), 'max', false);
            Y_eye = full(ind2vec(y_eye));
            
            R_eye = cov_shrink(X_eye');
            
            % assuming more resting observations that blinks
            % take the minimum of twice the blink dataset and the rest
            % trials
            X_blk = X(:, y == 5);
            y_blk = y(:, y == 5);
            
            [X_rb, y_rb] = equalizeLabels(cat(2, X_rst, X_blk, X_blk), cat(2, y_rst, y_blk, y_blk), 'min', false);
            Y_rb = full(ind2vec(y_rb));
            
            %% compute unmixing weights for the eye movements
            X_lr = X_eye(:, y_eye == 1 | y_eye == 2);
            Y_lr = Y_eye(2, y_eye == 1 | y_eye == 2);
            w_lr = obj.plr_elasticnet(X_lr, Y_lr, Rn);
            
            X_ud = X_eye(:, y_eye == 3 | y_eye == 4);
            Y_ud = Y_eye(4, y_eye == 3 | y_eye == 4);
            w_ud = obj.plr_elasticnet(X_ud, Y_ud, Rn);
            
            W_eye = cat(2, w_ud, w_lr);
            
            R_eye_y = W_eye'*R_eye * W_eye;
            A_eye = R_eye * W_eye / R_eye_y;
            
            %% apply eye movement correction
            
            % compute the correction matrix
            C_eye = (eye(nbchan) - A_eye*W_eye');
            
            X_rb_c = C_eye*X_rb;
            R_rb_c = cov_shrink(X_rb_c');

            X_rst_c = C_eye*X_rst;
            Rn_c = cov_shrink(X_rst_c');
            
            %% compute unmixing matrix for blink correction

            w_blk = obj.plr_elasticnet(X_rb_c, Y_rb(5,:), Rn_c);

            r_blk_y = w_blk' * R_rb_c * w_blk;
            a_blk = R_rb_c * w_blk / r_blk_y;

            % compute the unmixing matrix

            obj.W = cat(2, W_eye, w_blk);
            obj.A = cat(2, A_eye, a_blk);            
            
            % compute one step reconstruction matrix
            obj.C = (eye(nbchan) - a_blk * w_blk')*C_eye;
        end
        %%
        function data = apply(obj, data)
            % perform the correction step
            data(obj.eeg_chan_idxs,:) = obj.C * data(obj.eeg_chan_idxs, :);
        end
        
        %%
        function [w_plr, b_plr] = plr_elasticnet(obj, X, labels, Rn)
            % elastic net regularized penalized logistic regression
            %
            % optimization problem:
            % w_plr,b_plr = argmin_w,b sum_t log (1 + exp(-y_t*z_t))
            %                          + lambda_l2/2*||w||^2_Rn 
            %                          + lambda_l1*||w||_1
            %
            % with z_t = w'x_t+b
            %
            
            [ndim, nsamples] = size(X);
            
            % define simgoid function
            f = @(x) 1./(1 + exp(-x));
            
            % define l1 loss proxmap
            proxmap = @(x, t) max(0, abs(x) - t*obj.plr_lambda_l1).*sign(x);
            
            % define l2 loss
            Lambda = obj.plr_lambda_l2 * Rn;
            
            % initialize weight vars
            w_plr = zeros(ndim, 1);
            b_plr = 0;
            w_plr_old = w_plr;
            b_plr_old = 0;
            dw_norm = 1;
            w_norm = 1;
            
            fprintf('PLR working')
            
            % preallocate variables
            P = speye(nsamples);
            diag_idxs = find(P);
            
            % iterate until convergence
            iter = 1;
            while dw_norm/w_norm > obj.plr_tol
                
                if iter > obj.plr_maxiter
                    warning('Warning: aborting due to maxiter constraint');
                    w_plr = nan(ndim,1);
                    break
                end
                
                % compute the momentum for this iteration
                beta = (iter-1)/(iter+2);
                % and apply it
                w_plr_ = w_plr + beta * (w_plr - w_plr_old);
                b_plr_ = b_plr + beta * (b_plr - b_plr_old);
                w_plr_old = w_plr;
                b_plr_old = b_plr;
                
                % compute the gradient
                p = f(w_plr_' * X + b_plr);
                g = 1/nsamples * X * (labels - p)' - Lambda * w_plr_;
                g_b = 1/nsamples * sum(labels-p);
                
                % estimate the Lipschitz constant every 10th iteration
                if mod(iter, 10) == 1
                    P(diag_idxs) = p.*(1-p);
                    H = 1/nsamples * X * P * X' + Lambda;
                    H_b = 1/nsamples * sum(diag(P));
                    
                    % compute the Lipschitz constant
                    s = svd(H);
                    L = max(s);
                end
                
                % make a descent step
                w_plr = proxmap(w_plr_ + 1/L*g, 1/L);
                b_plr = b_plr_ + 1/H_b*g_b;
                
                % compute the norms
                dw = w_plr - w_plr_old;
                dw_norm = norm(dw);
                db = b_plr - b_plr_old;
                db_norm = norm(db);
                
                w_norm = norm(w_plr);
                
                if mod(iter, 10) == 1
                    fprintf('.');
                end
                
                if mod(iter, 100) == 0
                    c = (labels - 0.5)*2;
                    e_l2 = w_plr'*Lambda*w_plr;
                    e_l1 = obj.plr_lambda_l1*sum(abs(w_plr));
                    e_lr = 1/nsamples * sum(-log(f(c.*(w_plr' * X))));
                    
                    if obj.verbose
                        fprintf('||dw||/||w||=%e, L=%e, ||db||=%e\n', dw_norm/w_norm, L, db_norm);
                        fprintf('e_lr=%e, e_l2=%e, e_l1=%e\n', e_lr, e_l2, e_l1);
                    else
                        fprintf('||dw||/||w||=%f', dw_norm/w_norm);
                    end
                end
                iter = iter +1;
            end
            fprintf('DONE.\n')
            
            
        end
    end
    
end

