function[th]=gaussmix(data,Rat,tr)
%Fits average frequency values into a two mixed gaussian model and
%determines the cutoff frequency to separate gaussians.

if size(data,1)< size(data,2)
 data=data.';   
end

% Fit a gaussian mixture model
% obj = fitgmdist(data,2);
if length(data)>2 %Only when having more than 2 detections.
    obj = fitgmdist(data,2,'SharedCovariance',true);
    % try
    %     obj = fitgmdist(data,2);
    % catch exception
    %     disp('There was an error fitting the Gaussian mixture model')
    % %     error = exception.message
    %     obj = fitgmdist(data,2,'SharedCovariance',true);
    % 
    % end
    idx = cluster(obj,data);
    cluster1 = data(idx == 1,:);
    cluster2 = data(idx == 2,:);

    %Uncomment if running for the first time:
    
    %Find threshold value
%     th=min([max(cluster1) max(cluster2)]);

end

%These values were precomputed with the line above.
if Rat==26
    th=156.7; %Rat 26
end
% 
if Rat==27
    th=158; %Rat 27
end
% 
if Rat==24
    if tr(2)==40
        th=153.5; %Rat 24 (40)
    else
        th=153.6; %Rat 24 (35)
    end
end

%Use 155 threshold for swr disruption extra rats
if Rat ~=26 & Rat ~=27 & Rat ~=24
%      th=min([max(cluster1) max(cluster2)]);
th=155;
end

%  th=156.7; %Rat 26
% th=158; %Rat 27
% th=153.5; %Rat 24 (40)
% th=153.6; %Rat 24 (35)
end