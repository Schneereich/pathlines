function [x,y] = bahn(x0,y0,T,N,Z,U,V,X,Y)
%% BAHN  Return path coordinates of particle
%
% Assumptions: tbd.
%
% Inputs:
%   x0, y0 coordinates of starting point
%   X, Y in R^(n x m) equidistant grid
%   U, V cell-array with N entries as velocity profile
%   T time span
%   N number of time steps
%   Z number of intermediate steps between two time steps
% Outputs:
%   x, y path coordinates of particle

%
% $Date: February 5, 2020
% ________________________________________

% step size
dt = T/N; 
% additional intermediate evaluation points
h = dt/(Z+1);
x = zeros(1,N+Z*(N-1));
y = zeros(1,N+Z*(N-1));

x(1) = x0;
y(1) = y0;

% time steps to validate particle positions for
t = linspace(0,T,N+Z*(N-1));

for k = 1:length(t)-1
    tPrevIdx = floor(t(k)/dt)+1;
    tNextIdx = ceil(t(k)/dt)+1;
    gamma = (t(k)-(tPrevIdx-1)*dt)/dt;
    
    % in between to nodes
    if gamma ~= 0
        [up,vp] = evalVelocity(X,Y,U{tPrevIdx},V{tPrevIdx},x(k),y(k));
        % TODO: calculation error due to assumption / approximation.
        [un,vn] = evalVelocity(X,Y,U{tNextIdx},V{tNextIdx},x(k),y(k));
        uk = up + gamma*(un-up);
        vk = vp + gamma*(vn-vp);
    % directly on node
    else
       [uk,vk] = evalVelocity(X,Y,U{tPrevIdx},V{tPrevIdx},x(k),y(k));
    end
    
    % case 1: particle exceeded grid boundaries
    if isnan(uk)
        % particle exceeded yMax boundary
        if y(k) > Y(end,1)
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            x(k) = x(k-1) + (Y(end,1)-y(k-1))/s;
            y(k) = Y(end,1);
        % particle exceeded xMax boundary
        elseif x(k) > X(1,end)
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            y(k) = s*(X(1,end)+x(k-1)) + y(k-1);
            x(k) = X(1,end);
        % particle exceeded yMin boundary
        else
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            x(k) = x(k-1) + (Y(1,1)-y(k-1))/s;
            y(k) = Y(1,1);          
        end
        % x, y vectors - cut off residues
        x = x(1:k);
        y = y(1:k);
        return
    end
    
    % case 2: particle in inner circle
    if isinf(uk)
        % x, y vectors - cut off residues
        x = x(1:k);
        y = y(1:k);
        return
    end
    
    % case 3: particle position valid
    % explicit Euler
    x(k+1) = x(k) + h * uk;
    y(k+1) = y(k) + h * vk;
end
end