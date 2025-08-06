function g = sigmoid(z,b)
%SIGMOID Compute sigmoid functoon
%   J = SIGMOID(z) computes the sigmoid of z.

% You need to return the following variables correctly 
g = zeros(size(z));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the sigmoid of each value of z (z can be a matrix,
%               vector or scalar).



%g = 1./(1+exp(-z));
g = 1./( 1 + exp( (-z+0.5)/b) );
%g = 1./( 1 + exp( (-z+0.5) ) );
% =============================================================

end
