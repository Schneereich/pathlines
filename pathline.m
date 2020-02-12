function [x1,y1] = pathline(u1,v1,u2,v2,x,y,x0,y0,T,N)
  dx = x(2)-x(1); sx = length(x); xmx0 = x0-x(1);
  dy = y(2)-y(1); sy = length(y); ymy0 = y0-y(1);
  dt = T/(N-1);
  x1 = zeros(N,1); y1 = zeros(N,1);
  x1(1) = xmx0; y1(1) = ymy0;
  for k = 1 : N-1
    [up, vp] = getValueMonic(u1,v1,dx,dy,sx,sy,xmx0,ymy0);
    [un, vn] = getValueMonic(u2,v2,dx,dy,sx,sy,xmx0,ymy0);
    
    % TODO: Interpolation
    uk = ...%((N-k)*up + k*un);
    vk = ...%((N-k)*vp + k*vn);
    
    x1(k+1) = x(k) + dt * uk;
    y1(k+1) = y(k) + dt * vk;
  end
  x1 = x1 + x(1);
  y1 = y1 + y(1);
end

function [u,v] = getValueMonic(U,V,dx,dy,sx,sy,xmx0,ymy0)
  [fx,lx] = mapIndex(dx,xmx0);
  [fy,ly] = mapIndex(dy,ymy0);
  if (fx > 0 && fx < sx && fy > 0 && fy < sy)
    u = ly'*U(fy:fy+1,fx:fx+1)*lx;
    v = ly'*V(fy:fy+1,fx:fx+1)*lx;
  else
    u = NaN;
    v = NaN;
  end
end

function [fac,lambda] = mapIndex(dx,x)
  quotient = x/dx;
  fac = ceil(quotient);
  lambda = fac - quotient;
  if fac == 0 && lambda == 0 
    fac = fac + 1;
    lambda = 1;
  end
  lambda = [lambda;1-lambda];
end

