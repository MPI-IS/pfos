function Phi12 = sdePhi12Matrix(F, L, ssq, h)
% SDEPHI12MATRIX - Computes the Phi21 matrix on the big matrix exponential

% sdePhi12Matrix.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-01-26
% Version: 0.1

q = max(determineQs(F));

% LQL = L * diag(ssq) * L.';
M12 = h * L * diag(ssq) * L.';
M11 = h * F;

Phi12 = M12;

for k=1:2*q
  M12 = h/(k+1) * (M11 * L * diag(ssq) * L.' - M12 * F.');
  M11 = h/(k+1) * M11 * F;
  
  Phi12 = Phi12 + M12;
end

end 