function [coord, process] = myoptimize(P, no_dims, condition)
% this function is inspired from tsne_p. We use the gradient descent
% method to optimize our target function specified in our paper.
% condition can be loose or tight, we suggest using loose condition 
% for dataset with over 10000 cells

    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 3;
    end
    if ~exist('condition', 'var')
        condition = 'tight';
    end
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        coord = no_dims;
        no_dims = size(coord, 2);
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    n = size(P, 1);                                     % number of instances
    momentum = 0.5;                                     % initial momentum
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 250;                              % iteration at which momentum is changed                             % iteration at which lying about P-values is stopped
    max_iter = 1000;                                    % maximum number of iterations
    epsilon = 1000;                                     % initial learning rate
    min_gain = .01;                                     % minimum gain for delta-bar-delta
    process = zeros([n, 0]);                            % to save process
    
    % Make sure P-vals are set properly
    P(1:n + 1:end) = 0;                                 % set diagonal to zero
    P = 0.5 * (P + P');                                 % symmetrize P-values
    P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
    const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
    
    % Initialize the solution
    temp=clock;
    temp=sum(temp(4:6))*sum(temp(2:3));
    temp=round(temp/10);
    rand('seed',temp);
    if ~initial_solution
        if no_dims == 3
            x = (rand(n, 1) -0.5) * 50;
            y = (rand(n, 1) -0.5) * 50;
            z = (rand(n, 1) -0.5) * 50;
        	coord = [x, y, z];
        elseif no_dims == 2
            x = (rand(n, 1) -0.5) * 50;
            y = (rand(n, 1) -0.5) * 50;
            coord = [x, y];
        else
            coord = [];
            for i = 1 : no_dims
                coord = [coord, (rand(n, 1) -0.5) * 50];
            end
        end
    end
    incs  = zeros(size(coord));
    gains = ones(size(coord));
    
    % Run the iterations
    for iter=1:max_iter
        
        % Compute joint probability that point i and j are neighbors
        sum_current = sum(coord .^ 2, 2);
        d = bsxfun(@plus, sum_current, bsxfun(@plus, sum_current', -2 * (coord * coord')));
        num = 1 ./ (1 + d);                                      % Student-t distribution
        num(1:n+1:end) = 0;                                                 % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
        
        % Compute the gradients (faster implementation)
        P_Q = P - Q;
        P_Q(P_Q > 0 & d <= 0.01) = -0.01;
        L = P_Q .* num;
        grads = 4 * (diag(sum(L, 1)) - L) * coord;
            
        % Update the solution
        gains = (gains + .2) .* (sign(grads) ~= sign(incs)) ...         % note that the grads are actually -grads
              + (gains * .8) .* (sign(grads) == sign(incs));
        gains(gains < min_gain) = min_gain;
        incs = momentum * incs - epsilon * (gains .* grads);
        coord = coord + incs;
        coord = bsxfun(@minus, coord, mean(coord, 1));
         
        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        
        % Print out progress
        if ~rem(iter, 10)
            cost = const - sum(P(:) .* log(Q(:)));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
        end
        % rescale to 50
        range = max(abs(coord(:)));
        if strcmp(condition, 'tight')
            if range > 50 && ~rem(iter, 10)
                coord = coord .* 50/range;
            end
        else
            if range > 50 && ~rem(iter, 1000)
                coord = coord .* 50/range;
            end
        end
        process = [process, coord]; %#ok<AGROW>
    end
end
    
