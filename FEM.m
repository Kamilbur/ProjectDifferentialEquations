%
% FEM by Kamil Burkiewicz
% Made for Differential Equations Classes
%
function projekt(n)
tic
    % n -- number of elements for interval to be divide on
    %
    % Problem - find approximation of u(x), that satisfies:
    %             -(a u')' + b u' + c u = f
    % 
    %                     u(0) = 0
    %              u'(right) - u(right) = 0
    %
    % Approximation will be made using Finite Elements Method.
    %
    
        %%%%%%%%%%%%%%%%
        % coefficients %
        %%%%%%%%%%%%%%%%
    a = @(x) 1;
    b = @(x) 0;
    c = @(x) -1;
    f = @(x) 0;
    
        %%%%%%%%%%%%%
        % constants %
        %%%%%%%%%%%%%
    % This constant states for number of Legendre's polynomial used to
    % approximations of integrals (nodes and weights are precomputed at the
    % beginning of the program). Error is small due to specification of
    % base functions.
    pointsGaussLegendre = 3;
    % Boundaries
    left = 0;
    right = 2;
    % Step
    h = (right - left) / n;
    
    % Values necessary for Gauss-Legendre quadrature
    [nodes, weights] = computeGaussLegendreCoefficients(pointsGaussLegendre);
    
    
        %%%%%%%%%%%%%%%%%%
        % filling matrix %
        %%%%%%%%%%%%%%%%%%
    % Arrays to create sparse matrix. This approach lead to speed up of
    % computations and is recomended by MATLAB documentation.
    I = zeros(1, 3 * (n + 1) - 2);
    J = zeros(1, 3 * (n + 1) - 2);
    V = zeros(1, 3 * (n + 1) - 2);
    
    funcToIntegrate = @(i,j,x) a(x) * e_der(i,n,x,left,right) .* e_der(j,n,x,left,right) + ...
                     b(x) * e_der(i,n,x,left,right) * e(j,n,x,left,right) + ...
                     c(x) * e(i,n,x,left,right) * e(j,n,x,left,right);
    for i = 1:n + 1
        for j = i - 1:i + 1
            if (j >= 1 && j <= n + 1)
                % 3 * (i - 1) + j - i + 1 -- this transforms pair (i,j) to
                % the number of the pair in ordering provided by these for
                % loops f.e.   (1,0) -> no idx, (1,1) -> 1, (1,2) -> 2, 
                % (2,1) -> 3, (2,2) -> 4 and so on
                idx = 3 * (i - 1) + j - i + 1;
                I(idx) = i;
                J(idx) = j;
                leftIntegrationBound  = min((i - 1) * h, (j - 1) * h);
                rightIntegrationBound = max((i + 1) * h, (j + 1) * h);
                
               
                V(idx) = GaussianLegendreQuadrature(@(x) funcToIntegrate(i - 1,j - 1,x), ...
                    leftIntegrationBound, rightIntegrationBound, nodes, weights) - ...
                    a(right) * e(i - 1,n,right,left,right) * e(j - 1,n,right,left,right);
            end
        end
    end
    % make sparse matrix from values calculated in loop above 
    A = sparse(I,J,V);
    
        %%%%%%%%%%%%%%%%%%
        % filling vector %
        %%%%%%%%%%%%%%%%%%    
    B = zeros(n + 1, 1);
    funcToIntegrate = @(i,x) f(x) * e(i,n,x,left,right);
    for i = 1:n + 1
        leftIntegrationBound  = (i - 1) * h;
        rightIntegrationBound = (i + 1) * h;
        B(i) = GaussianLegendreQuadrature(@(x) funcToIntegrate(i - 1, x),...
                leftIntegrationBound, rightIntegrationBound, nodes, weights);
    end

        %%%%%%%%%%%%%%%%
        % coefficients %
        %%%%%%%%%%%%%%%%
    
    % finding coefficients is just solving system of linear equations of the form
    % A X = B, where A and B where values in these matices was computed earlier.
    coeff = A \ B;
    
    draw(coeff, left, right, n);
    
toc % tic at the beginning and toc at the end for measurement of running time 
end



function [nodes, weights] = computeGaussLegendreCoefficients(n)
    % n - number of points to be used in later quadrature
    
    p = zeros(n + 1);
    p(1,end) = 1;
    p(2,end - 1) = 1;
    
        % Calculate n-th Legendre polynomial bottom-up method
    for i = 3:n + 1
        p(i, end) = -(i - 2) * p(i - 2, end) / (i - 1);
        for j = 1:n
                % Formula: n * P_n(x) = (2 * n - 1) * x * P_(n-1) - (n - 1) * P_(n-2)
            p(i, n + 1 - j) = ((2 * i - 3) * p(i - 1, n + 1 - j + 1) - (i - 2) * p(i - 2, n + 1 - j)) / (i - 1);
        end
    end

    L_n = p(end,:);  % n-th Legendre's polynomial
    L_n_derivative = polyder(L_n);
    
        % nodes are roots of Legendre's polynomial
    nodes = transpose(roots(L_n));
    
        % computation of weights
    weights = zeros(1, length(nodes));
    
        % w_i = 2 / ((1 - x_i ^2) * L'_n(x_i))
        % x_i is an i-th element of nodes array
    temp = 1 - nodes.^2;
    for i = 1:length(nodes)
        weights(i) = 2 / (temp(i) * polyval(L_n_derivative, nodes(i))^2);
    end
end



function g = GaussianLegendreQuadrature(integrandFun, l, r, nodes, weights)
    % Quadrature is made using number of points specified in constants.
    % nodes and weights have to be passed as a function argument.
    % integrandFun -- function to integrate
    % l -- left boundary of interval
    % r -- right boundary of interval
    % nodes -- list of nodes necessary to Gauss-Legendre quadrature
    % weithgs -- list of weights necessary to Gauss-Legendre quadrature
    
    n = length(nodes);
    sum = 0;
    
        % if interval is different from [-1,1] we have to translate it
        % using affine transform
    a = 1;
    b = 0;
    if (l ~= -1 || r ~= 1)
        a = (r - l) / 2;
        b = (r + l) / 2;
    end
    
    for i = 1:n
        sum = sum + weights(i) * integrandFun(a * nodes(i) + b);
    end
    
    g = a * sum;
end



function e = e(k, n, x, l, r)
    % l -- left boundary of interval
    % r -- righht boundary of interval
    h = (r - l) / n;
    e = max((1 - abs(x - k * h) /h), 0);
end



function e_der = e_der(k, n, x, l, r)
    % derivatives of e functions
    
    h = (r - l) / n;
    val = zeros(1, length(x));
    for i = 1:length(x)
        if (x(i) >= (k - 1) * h && x(i) < k * h)
            val(i) = 1 / h;
        elseif (x(i) >= k * h && x(i) < (k + 1) * h)
            val(i) = -1 / h;
        end
    end
    
    e_der = val;
end



function draw(coeff, l, r, n)
    u = @(x) e_linear_combination(coeff, n, x, l, r);
    points = l:(r-l)/n:r;
    values = zeros(n + 1, 1);
    
    for i = 1:n + 1
        values(i) = u(points(i));
    end

    plot(points,values,'rx');
    hold on;
    plot(points,values);
    hold off;

end



function lc = e_linear_combination(coeff, n, x, l, r)
    % linear combination of base functions e determined by coefficients
    % from coeff list
    sum = 0;
    for i = 1:length(coeff)
        sum = sum + coeff(i) * e(i, n, x, l, r);
    end
    
    lc = sum;
end
