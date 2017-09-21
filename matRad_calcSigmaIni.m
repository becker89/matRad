function [sigmaIni] = matRad_calcSigmaIni(baseData,rays)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates simultaneously the initial sigma of the beam for
% one or more energies
%
% call
%   sigmaIni = matRad_calcSigmaIni(machine.data,stf(i).ray,stf(i).ray(j).SSD);
%
% input
%   baseData:           'machine.data' file
%   rays:            	'stf.ray' file
%   SSD:                source-surface difference
%
% output
%   sigmaIni:        	initial sigma of the ray at certain energy (or
%                       energies). The data is given in 1xP dimensions,
%                       where 'P' represents the number of different
%                       energies
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is not part of the offical matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% helper function for energy selection
round2   = @(a,b)round(a*10^b)/10^b;

sigmaIni = zeros(1,numel(rays.focusIx));

for i = 1:numel(rays.focusIx)
   energyIx = max(round2(rays.energy(i),4)) == round2([baseData.energy],4);

   sigmaIni(i) = matRad_interp1(baseData(energyIx).initFocus.dist(rays.focusIx(i),:)',...
                  baseData(energyIx).initFocus.sigma(rays.focusIx(i),:)',rays.SSD);

end

end

