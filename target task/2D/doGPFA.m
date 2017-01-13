%-------------------------------------------------------------------------%
%%%% DEFINE THESE FOLLOWING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_conds = [];
%{indices1,indices2};  
% e.g. {[1 3 5],[2 4 6]} means that trials [1 3 5] are from condition 1 and 
% trials [2 4 6] are from condition 2. There can be any number of different
% conditions. If you don't want to separate by condition, set empty 
% (diff_conds = []). 
% *This is for plotting purposes only, the dimensionality
% reduction process uses all (unlabeled) data.*

trial_table = alldays(1).tt; %PUT_TRIAL_TABLE_HERE; 
% Trial table size with rows representing trials

start_column = 5;%PUT_START_COLUMN_HERE; 
% e.g. start_column = 5  means the trajectories will begin aligned to the 
% times listed in column 5 of the trial table

start_offset = 0;%PUT_START_OFFSET_HERE; 
% e.g. start_offset = 200 means the trajectories will start 200ms after the
% times in start_column

end_column = 7;%PUT_END_COLUMN_HERE; 
% Same as start_column

end_offset = 0;%PUT_END_OFFSET_HERE; 
% Same as start_offset

units = alldays(1).PMd_unist;%PUT_UNITS_HERE; % Cell array in which each cell contains spike times
%                         from a single neuron
runIdx = 'COTEST';%'PUT_LABEL_HERE'; % Some identifying tag for the run 
%                            (e.g. 'Target_On to Go_Cue 01_30_2014')
directory = 'C://Users/limblab/Desktop/testing123';%'PUT_DIRECTORY_HERE'; % Directory string for saving data.

% Cross validation takes a long time. For exploratory types of analysis,
% leave this as false.
do_cross_val = false; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
[trial_rast,dat] = trial_raster(units,trial_table,[start_column start_offset],...
                                                  [end_column,end_offset]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===========================================
% 1) Basic extraction of neural trajectories
% ===========================================
method = 'gpfa';

% Select number of latent dimensions
xDim = 8; %find optimal using CV (below section)
kernSD = 30; % find optimal kernal using CV (below section)

% Extract neural trajectories
result = neuralTraj_limblab(runIdx, dat, directory, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD);

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);

% Seq train will NOT be ordered by trial. This makes labeling by trial 
% with indices quite difficult, so we reorder it to be sequential. We also
% keep the original (orig_SeqTrain) in case a shuffled form is desired 
ord_seqTrain = seqTrain;
for tr = 1:length(seqTrain)
    tId = seqTrain(tr).trialId; 
    ord_seqTrain(tId) = seqTrain(tr);
end
orig_seqTrain = seqTrain; seqTrain = ord_seqTrain;

% Plot neural trajectories in 3D space
plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
%plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3,'nPlotMax',1000);

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);

fprintf('\nDone\n');

%% Plot traces 
if ~isempty(diff_conds)
    figure; hold on; 
    cols2plot = distinguishable_colors(length(diff_conds));

    cond_seqTrain = cell(length(diff_conds),1);
    for j = 1:length(diff_conds)
        cond_seqTrain{j} = seqTrain(diff_conds{j});
        plot3D_addon(cond_seqTrain{j}, 'xorth', cols2plot(j,:),...
                    'dimsToPlot', 1:3,'nPlotMax',10000);
    end
end

%%
if do_cross_val
% ========================================================
% 2) Full cross-validation to find:
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================

% Select number of cross-validation folds
numFolds = 4;

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
for xDim = [2 5 8]
  %neuralTraj(runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
  %neuralTraj(runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
  %neuralTraj(runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
end
fprintf('\n');
% NOTES:
% - These function calls are computationally demanding.  Cross-validation 
%   takes a long time because a separate model has to be fit for each 
%   state dimensionality and each cross-validation fold.

% Plot prediction error versus state dimensionality.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
kernSD = 30; % select kernSD for two-stage methods
plotPredErrorVsDim(runIdx, kernSD);
% NOTES:
% - Using this figure, we i) compare the performance (i.e,,
%   predictive ability) of different methods for extracting neural
%   trajectories, and ii) find the optimal latent dimensionality for
%   each method.  The optimal dimensionality is that which gives the
%   lowest prediction error.  For the two-stage methods, the latent
%   dimensionality and smoothing kernel width must be jointly
%   optimized, which requires looking at the next figure.
% - In this particular example, the optimal dimensionality is 5. This
%   implies that, even though the raw data are evolving in a
%   53-dimensional space (i.e., there are 53 units), the system
%   appears to be using only 5 degrees of freedom due to firing rate
%   correlations across the neural population.
% - Analogous to Figure 5A in Yu et al., J Neurophysiol, 2009.

% Plot prediction error versus kernelSD.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
xDim = 8; % select state dimensionality
plotPredErrorVsKernSD(runIdx, xDim);
% NOTES:
% - This figure is used to find the optimal smoothing kernel for the
%   two-stage methods.  The same smoothing kernel is used for all units.
% - In this particular example, the optimal standard deviation of a
%   Gaussian smoothing kernel with FA is 30 ms.
% - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.
end