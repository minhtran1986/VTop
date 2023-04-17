function [P,H] = Filter(element, node, fac)
%% Purpose: compute the weight factors for density filter
%%
numElement = size(element,1);
H = sparse(numElement, numElement);
P = sparse(numElement, numElement);
eCenter = sparse(numElement,2);
radius = sparse(numElement,1);
for ii = 1:numElement
    % center of element i
    iconn = element{ii};
    nn = length(iconn);
    iCenter = [sum(node(iconn,1))/nn, sum(node(iconn,2))/nn];
    eCenter(ii,1:2) = iCenter; % coordinates of element center
    % Distance from nodes to center
    r = sqrt( (node(iconn,1) - iCenter(:,1)).^2 + (node(iconn,2) - iCenter(:,2)).^2);
    radius(ii) = mean(r)*fac; % filter radius
end

for ii = 1:numElement
    rmin = radius(ii);
    iCenter = eCenter(ii,:);
    %% Distance from the iCenter to the centers of other elements
    r = sqrt( (eCenter(:,1) - iCenter(:,1)).^2 + (eCenter(:,2) - iCenter(:,2)).^2);
    %%
    test = rmin - r;
    Jdomain = find(test >= 0);
    tmp = (rmin - r(Jdomain))/rmin;
    H(ii, Jdomain') = tmp';
    P(ii,:) = H(ii,:)/sum(tmp);
end
% P = H./ sum(H,2); % this is only applicable for R2021
end