function [wQ] = matRad_weightedQuantile(values, percentiles, weight, isSorted, extraPol)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty analysis report generaator function
% 
% call
%   latexReport(ct, cst, pln, nominalScenario, structureStat, param)
%
% input
%   values:             random variable vector
%   percentiles:        percentiles to be calculated
%   weight:             (optional) weight vector (same length as values)
%   isSorted:           (optional) bool: are the values sorted alreay?

% output
%   (binary)            a pdf report will be generated and saved
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % check input arguments
  if ~exist('sortedB', 'var') || isempty(isSorted)
    isSorted = false;
  end
  if ~exist('weight', 'var') || isempty(weight)
    weight = ones(size(values));
  end
  if ~exist('extraPolB', 'var') || isempty(extraPol)
    extraPol = false;
  end
  if ~(percentiles(:) >= 0 & percentiles(:) <= 1)
    error('Quantiles must not be outside [0, 1]');
  end

  % sort values
  if ~isSorted
    [values, sortIx] = sort(values);
    weight = weight(sortIx);
  end
  
  wQtemp = cumsum(weight) - 0.5 * weight;
  wQtemp = wQtemp ./ sum(weight);

  wQ = NaN * ones(size(values,1), 2);

  [x,ia] = unique(wQtemp);

  V = values(ia);
  F = griddedInterpolant(x, V);
  wQ = F(percentiles);

  
  wQ = wQ';

end % eof
