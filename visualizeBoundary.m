function [innerX] =  visualizeBoundary(X, y, model, varargin)
    %   VISUALIZEBOUNDARY gives pol inner approx of non-linear boundary learned by the SVM

    % Plot the training data on top of the boundary
    % plotData(X, y)

    % Make classification predictions over a grid of values
    x1plot = linspace(min(X(:,1)), max(X(:,1)), 500)';
    x2plot = linspace(min(X(:,2)), max(X(:,2)), 500)';
    [X1, X2] = meshgrid(x1plot, x2plot);
    vals = zeros(size(X1));


    %% Need the specific test point flags
    % These test points are then used to approximate the polytope
    % which is the inner approx. of SVM curve and conOld intersect.
    innerX = []; 

    for i = 1:size(X1, 2)
       this_X = [X1(:, i), X2(:, i)];
       vals(:, i) = predict(model, this_X);     
       for k = 1:size(vals,1)
            if vals(k,i) == 1
                innerX =  [innerX ; this_X(k,:)];
            end
       end  
    end

end
