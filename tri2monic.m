function ugrid = tri2monic(coordinates,elements,u,x,y)
%*** interpolate piecewise affine function onto grid (x,y)
sx = x(1,:); sy = y(:,1);
m = size(x,2); n = size(x,1);
hx = sx(2)-sx(1); hy = sy(2)-sy(1);
 
nS = length(u); % number of solutions
ugrid = cell(1, nS);
for j = 1:nS
  ugrid{j} = NaN(size(x,1),size(x,2));
end

%*** Interpolation for each triangle
for j = 1:size(elements,1)
  p = coordinates(elements(j,:),:);
        
  %*** Discretization for x-direction 
  xminDreieck = max(1,floor((min(p(:,1)) - sx(1))/hx)+1);
  xmaxDreieck = min(m,ceil( (max(p(:,1)) - sx(1))/hx)+1);
  dsx = sx(xminDreieck:xmaxDreieck);
     
  %*** Discretization for y-direction
  yminDreieck = max(1,floor((min(p(:,2)) - sy(1))/hy)+1);
  ymaxDreieck = min(n,ceil( (max(p(:,2)) - sy(1))/hy)+1);
  dsy = sy(yminDreieck:ymaxDreieck);
     
  %*** Extract solutions 
  for k = 1:nS
    sol(:,k) = (u{k}(j,:))';
  end
  %*** Create submesh
  [dx,dy] = meshgrid(dsx,dsy);
  %*** Interpolate on submesh
  tgrid = tri2grid(p,sol,dx,dy);
  %*** Incorporate submesh into global mesh 
  for k = 1:nS
    ugrid{k}(yminDreieck:ymaxDreieck,xminDreieck:xmaxDreieck) = ...
        min(tgrid{k},ugrid{k}(yminDreieck:ymaxDreieck,xminDreieck:xmaxDreieck));
  end
end

function ugrid = tri2grid(P,u,x,y)
% interpolate from triangle with vertices P to monic grid (x,y)
% x, y, and all ugrid's{:} have to have the same size
[m,n] = size(x);
for k = 1:size(u,2)
  ugrid{k} = NaN(m,n);
end
x = x(:); y = y(:);
area2 = (P(2,1)-P(1,1))*(P(3,2)-P(1,2))...
       -(P(2,2)-P(1,2))*(P(3,1)-P(1,1));
L = zeros(m*n,3);
L(:,1) = (P(2,1)-x).*(P(3,2)-y) - (P(3,1)-x).*(P(2,2)-y);
L(:,2) = (P(3,1)-x).*(P(1,2)-y) - (P(1,1)-x).*(P(3,2)-y);
L(:,3) = area2 - L(:,1) - L(:,2);
tol = 1e-9 * area2;
if area2 > 0 
  ind = find (L(:,1)>=-tol & L(:,2)>=-tol & L(:,3)>=-tol);
else
  ind = find (L(:,1)<=-tol & L(:,2)<=-tol & L(:,3)<=-tol);
end

if ~isempty(ind)
  for j = 1:size(u,2)
    ugrid{j}(ind) = (L(ind,:)*u(:,j))./area2;
    ugrid{j} = reshape(ugrid{j},m,n);
  end
end