function [structureStat, doseStat, param] = matRad_samplingAnalysis(ct,cst,subIx,mRealizations,w, nominalScenario, multScen, param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty sampling analysis function
% 
% call
%   [structureStat, doseStat] = samplingAnalysis(ct,cst,subIx,mRealizations,w)
%
% input
%   ct:             ct cube
%   cst:            matRad cst struct
%   subIx           set of indices of the cube which are used for dose calculation
%   mRealizations   resulting dose of the individual scenarios
%   w               vector containing probabilities of the scenarios
%
% output
%   structureStat   structure-wise statistics (mean, max, percentiles, ...)
%   doseStat        dose-wise statistics (mean, max, percentiles, ...)
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

%% check integrity of statistics
param.sufficientStatistics = checkIntegrity(multScen, param);
%% calculate mean and std cube
doseStat.meanCube              = zeros(ct.cubeDim);
doseStat.stdCube               = zeros(ct.cubeDim);

doseStat.meanCubeW             = zeros(ct.cubeDim);
doseStat.stdCubeW              = zeros(ct.cubeDim);

doseStat.meanCube(subIx)       = mean(mRealizations,2);   
doseStat.stdCube(subIx)        = std(mRealizations,1,2);  
doseStat.meanCubeW(subIx)      = (sum(mRealizations * diag(w),2) );
doseStat.stdCubeW(subIx)       = std(mRealizations,w,2);

% gamma cube
if isfield(nominalScenario,'RBExD')
    doseCube = nominalScenario.RBExD;
    doseStat.gammaAnalysis.cube1Name = 'nominalScenario.RBExD';
else
    doseCube = nominalScenario.physicalDose;
    doseStat.gammaAnalysis.cube1Name = 'nominalScenario.physicalDose';
end
doseStat.gammaAnalysis.cube1 = doseCube;
doseStat.gammaAnalysis.cube2 = doseStat.meanCubeW;
doseStat.gammaAnalysis.cube2Name = 'doseStat.meanCubeW';
criteria = [2 2];
doseStat.gammaAnalysis.doseAgreement = criteria(1);
doseStat.gammaAnalysis.distAgreement = criteria(2);
doseStat.gammaAnalysis.gammaCube = matRad_gammaIndex(doseCube,doseStat.meanCubeW,[ct.resolution.x ct.resolution.y ct.resolution.z],criteria);

%% percentiles
percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
percentileNames = cell(numel(percentiles),1);
% create fieldnames
for i = 1:numel(percentiles)
    percentileNames{i} = ['P',num2str(percentiles(i)*100)];
end
% create table rownames
metric = vertcat({'mean';'min';'max';'std'},percentileNames{:});

% create statstics where structure based results (QI and DVH) are available
for i = 1:size(cst,1)
    structureStat(i).name = cst{i,2};
    structureStat(i).dvh = cst{i,9};
    structureStat(i).qi = cst{i,10};
    structureStat(i).percentiles = percentiles;
    structureStat(i).metric = metric;
    
    structureStat(i).dvhStat = calcDVHStat(cst{i,9},percentiles,w);
    structureStat(i).qiStat = calcQiStat(cst{i,10},percentiles,w);
end


% dvh statistics

    function dvhStat = calcDVHStat(dvh,percentiles,w)
        doseGrid = dvh{1}(1,:);
        dvhMat = NaN * ones(numel(dvh),size(dvh{1},2));
        for j = 1:numel(dvh)
            dvhMat(j,:) = dvh{j}(2,:);            
        end
        % for statistical reasons, treat NaN as 0
        dvhMat(isnan(dvhMat)) = 0;
        
        dvhStat.mean(1,:) = doseGrid;
        dvhStat.mean(2,:) = wMean(dvhMat,w);
        % [~,argmin] = min(dvhStat.mean(2,:));
        % dvhStat.mean(2,argmin + 1:end) = NaN;
        
        dvhStat.min(1,:) = doseGrid;
        dvhStat.min(2,:) = min(dvhMat);
        % [~,argmin] = min(dvhStat.min(2,:));
        % dvhStat.min(2,argmin + 1:end) = NaN;
        
        dvhStat.max(1,:) = doseGrid;
        dvhStat.max(2,:) = max(dvhMat);
        % [~,argmin] = min(dvhStat.max(2,:));
        % dvhStat.max(2,argmin + 1:end) = NaN;
        
        dvhStat.std(1,:) = doseGrid;
        dvhStat.std(2,:) = std(dvhMat,w);

        dvhStat.percDVH = NaN * ones(numel(percentiles),size(dvh{1},2));
        
        for j = 1:size(dvhMat,2)
            wQ =  matRad_weightedQuantile(dvhMat(:,j), percentiles, w', false, 'none');
            dvhStat.percDVH(:,j) = wQ;
        end

    end % eof calcDVHStat

    % qi statistics
    function qiStat = calcQiStat(qi,percentiles,w)
        fields = fieldnames(qi{1});
        % convert to struct (maybe change structure to struct later)
        qiStruct(1) = qi{1};
        for j = 2:numel(qi)
            qiStruct(j) = qi{j};
        end
        
        % create helper matlab structure which will be converted to table
        qiStatH = struct();
        for j = 1:numel(fields)
            if numel([qiStruct(:).(fields{j})]) == numel(w)
                qiStatH(1).(fields{j}) = wMean([qiStruct(:).(fields{j})],w);
                qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
                qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
                qiStatH(4).(fields{j}) = std([qiStruct(:).(fields{j})],w);
                wQ = matRad_weightedQuantile([qiStruct(:).(fields{j})], percentiles, w', false, 'none');
                for k = 1:numel(wQ)
                    sIx = k + 4;
                    qiStatH(sIx).(fields{j}) = wQ(k);
                end
            else
                for k = 1:(4 + numel(percentiles))
                    qiStatH(k).(fields{j}) = [];
                end
            end
        end
        qiStat = struct2table(qiStatH);
        qiStat.Properties.RowNames = metric;
    end % eof calcQiStat

    function S = wMean(X,w)
        if exist('w','var') || ~isempty(w)
            if isvector(X) && isvector(w)
                S = reshape(w,1,[]) * reshape(X,[],1) / sum(w);
            else
                % row-wise
                S = reshape(w,1,[]) * X ./ sum(w);        
            end

        else
            S = mean(X);
        end
    end

    % check integrity of scenario analysis (i.e. check number of scenarios)
    function statCheck = checkIntegrity(multScen, param)
        if exist('param','var') && isfield(param, 'sufficientStatistics') && ~param.sufficientStatistics
            statCheck = false;
            return
        end
        
        if multScen.numOfScen > 20
            totalNum = true;
        else
            totalNum = false;
        end
        
        statCheck = totalNum; % * .... *
        
    end

end
