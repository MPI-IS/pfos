function R = semichol (A, uplo)
% SEMICHOL   Cholesky-like decomposition of a positive semidefinite matrix.
%    SEMICHOL(A) returns a matrix R such that A = R'*R. If A is positive
%    definite then R is upper diagonal.
%
%    SEMICHOL(A, 'lower') returns R such that A = R*R'. The default behavior
%    is achived with SEMICHOL(A, 'upper').
%
%    SEMICHOL first attempt to perform the decomposition using CHOL. If
%    this fail A is decompsed into A = P*L*D*L'*P' using LDL. All negative 
%    elements of the diagonal matrix D is set to zero prior to returning
%    R = (P*L*sqrt(D))'. This is computatationally more demanding than the
%    Choleskey decomposition, but still significantly faster than an eigen
%    value decomposition.
%
%    See also CHOL, LDL, EIG

% Author: Soren Hauberg <soren.hauberg@tue.mpg.de>.
% License: public domain (i.e. use this code as you please).
% Feedback: if you have comments, questions or find bugs in the code,
% please contact me.

  %% Check inputs and set default parameter
  if (nargin < 1)
    error ('semichol: not enough input arguments');
  elseif (nargin < 2)
    uplo = 'upper';
  end % if
  
  %R = chol(A+1e-6*eye(size(A)), uplo); return
  
  %% Try first with a Cholesky decomposition
  [R, p] = chol (A, uplo);
  if (p ~= 0)
    %% Cholesky failed so now we try the more expensive LDL^T decomposition
    [L, D, P] = ldl (A, 'vector');
    D = sqrt (max (diag (D), 0)); 
    L = bsxfun (@times, L, D.');
    R (P, P) = L;
    if (strcmpi (uplo, 'upper'))
      R = R.';
    end % if
  end % if
end % function
