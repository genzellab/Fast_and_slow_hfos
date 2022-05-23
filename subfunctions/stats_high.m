
function [zmap]=stats_high(freq1,freq2,w)
ntrials=size(freq1.powspctrm,1);
%Requires converting NaNs values into zeros.
no1=freq1.powspctrm;
no2=freq2.powspctrm;

no1(isnan(no1))=0;
no2(isnan(no2))=0;
%%
freq1.powspctrm=no1;
freq2.powspctrm=no2;
%% statistics via permutation testing

% p-value
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));

% number of permutations
n_permutes = 500;

% initialize null hypothesis maps
permmaps = zeros(n_permutes,length(freq1.freq),length(freq1.time));

% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1" 
%   and trials ntrials+1:end are from channel "2"
tf3d = cat(3,reshape(squeeze(freq1.powspctrm(:,w,:,:)),[length(freq1.freq) length(freq1.time)... 
    ntrials ]),reshape(squeeze(freq2.powspctrm(:,w,:,:)),[length(freq1.freq) length(freq1.time)... 
    ntrials ]));

%concatenated in time.
% freq, time, trials

% generate maps under the null hypothesis
for permi = 1:n_permutes
    permi
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
end
%% show non-corrected thresholded maps

diffmap = squeeze(mean(freq2.powspctrm(:,w,:,:),1 )) - squeeze(mean(freq1.powspctrm(:,w,:,:),1 ));

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;
%%

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max


% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end
%%
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

%% now with max-pixel-based thresholding

% find the threshold for lower and upper values
thresh_lo = prctile(max_val(:,1),100-100*pval); % what is the
thresh_hi = prctile(max_val(:,2),100-100*pval); % true p-value?

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

end