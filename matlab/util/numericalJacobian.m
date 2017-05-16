function J = numericalJacobian(f, x, vectorized, delta)
% NUMERICAL_JACOBIAN   compute the Jacobian of a function using forward
%    differences.
%

  %% Check input
  if (nargin < 3)
    vectorized = false;
  end % if
  if (nargin < 4)
    delta = 1e-6;
  end % if
  
  %% Compute Jacobian based on finite differences
  F0 = f(x);
  D_in = numel(x);
  D_out = numel(F0);
  if (vectorized)
    % XXX: there seems to be a transposition bug here -- results conflict
    % with non-vectorized code
    state_D = bsxfun(@plus, x, delta*eye(D_in));
    FD = f(state_D);
    J = bsxfun(@minus, FD, F0) ./ delta;
  else
    J = NaN(D_out, D_in);
    %{
    for d = 1:D_in
      state_d = x;
      state_d(d) = x(d) + delta;
      Fplus = f(state_d);
      state_d(d) = x(d) - delta;
      Fminus = f(state_d);
      J (:, d) = (Fplus(:) - Fminus(:)) ./ (2*delta);
    end % for
    %}
    % {
    for d = 1:D_in
      state_d = x;
      state_d(d) = state_d(d) + delta;
      Fd = f(state_d);
      J(:, d) = (Fd(:) - F0(:)) ./ delta;
    end % for
    %}
  end % if
end % function
