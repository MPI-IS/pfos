function Phi11 = sdePhi11Matrix (F, h)
% SDEPHI11MATRIX - Computes the Phi_11 matrix on the big matrix exponential

% sdePhi11Matrix.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-01-26
% Version: 0.1

q = max(determineQs(F));
  
I = eye(size(F));
Phi11 = I + h/q * F;

for k=1:q-1
  Phi11 = I + h/(q-k) * F * Phi11;
end

end