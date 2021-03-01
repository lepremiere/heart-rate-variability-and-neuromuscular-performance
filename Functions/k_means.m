function [means, idx] = k_means(x)   
% This function finds two means within an array.
%
% USAGE
% E.g.:
% [means, idx] = k_means(x)  
%
% INPUT
% x:            Array (n x 1) for which means should be identified.
%
% OUTPUT
% means:        Array (1 x 2) containing the identified means.
% idx:          Array (n x 2) containing logicals to which mean the
%               values in input "x" belong.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------
     
centroids = [mean(x)-1 mean(x)+1];  % Initial centroids
condition = 0;
while condition == 0
    distances = [abs(x - centroids(1)), abs(x - centroids(2))]; % Euclidean distances to the centroids
    [~, min_dist_idx] = min(distances, [], 2);                  % Identifying the closest centroid
    new_means = [mean(x(min_dist_idx == 1)),...                 % Updating centroids
                 mean(x(min_dist_idx == 2))];   
    
    % Terminal condition to break loop
    if(new_means == centroids)
        break;
    else
        centroids = new_means;
    end
end
means = centroids;      % Final means
idx  = min_dist_idx;    % Indeces to which each data point belongs
end