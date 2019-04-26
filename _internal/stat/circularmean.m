function mu = circularmean(points, weights, angleidx)

% mu = circularmean(points, weights, angleidx)

weights = weights(:)'/sum(weights(:));
mu = sum(bsxfun(@times, points, weights),2);
for i=angleidx(:)'
    sinx = sum(sin(points(i,:)).*weights);
    cosx = sum(cos(points(i,:)).*weights);
    mu(i) = atan2(sinx, cosx);
end
