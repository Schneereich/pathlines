function [uk,vk] = evalVelocity(X,Y,U,V,xk,yk)
%% evalVelocity  Return velocities
%
% Assumptions: tbd.
%
% Inputs:
%   xk, yk coordinates of evaluated point
%   X, Y in R^(n x m) equidistant grid
%   U, V cell-array with N entries as velocity profile
% Outputs:
%   uk, vk velocities of particle

%
% $Date: February 5, 2020
% ________________________________________

% maye define this variables globally
m = size(X,2);
n = size(X,1);
dx = X(1,end)/m;
dy = Y(end,1)/n;

%% find element indices
leftIdx = floor(xk/dx)+1;
rightIdx = ceil(xk/dx)+1;
downIdx = floor(yk/dy)+1;
upIdx = ceil(yk/dy)+1;

%% case distinction
% case 1: particle outside of grid
if leftIdx < 1 || rightIdx > m || downIdx < 1 || upIdx > n
    warning("Particle is out of grid boundaries. Setting velocity to NaN");
    uk = NaN;
    vk = NaN;
    return
end

% case 2: particle inside inner circle
if any(isnan([U(leftIdx,downIdx),U(leftIdx,upIdx), ...
              U(rightIdx,downIdx),U(rightIdx,upIdx)]))
    warning(["Particle is inside inner circle." ...
             "Setting velocity to infinity."]);
    uk = inf;
    vk = inf;
    return
end
    
% case 3: particle position valid, interpolate values
alpha = (xk-X(1,leftIdx))/dx;
beta  = (yk-Y(downIdx,1))/dy;


%% interpolation between edge values
u1 = U(leftIdx,downIdx)+alpha*(U(rightIdx,downIdx)-U(leftIdx,downIdx));
u2 = U(leftIdx,  upIdx)+alpha*(U(rightIdx,  upIdx)-U(leftIdx,  upIdx));
v1 = V(leftIdx,downIdx)+alpha*(V(rightIdx,downIdx)-V(leftIdx,downIdx));
v2 = V(leftIdx,  upIdx)+alpha*(V(rightIdx,  upIdx)-V(leftIdx,  upIdx));


uk = u1+beta*(u2-u1);
vk = v1+beta*(v2-v1);

end

