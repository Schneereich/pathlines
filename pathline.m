function [x1,y1] = pathline(u1,v1,u2,v2,x,y,x0,y0,T,N,method)
dx = x(2)-x(1); sx = length(x); xmx0 = x0-x(1);
dy = y(2)-y(1); sy = length(y); ymy0 = y0-y(1);
dt = T/(N-1);
x1 = zeros(N,1); y1 = zeros(N,1);
x1(1) = xmx0; y1(1) = ymy0;
switch method
    % Explizites Eulerverfahren
    case 'euler'
        for k = 1 : N-1
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [u,v] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            x1(k+1) = x1(k) + dt*u;
            y1(k+1) = y1(k) + dt*v;
        end
        
    % Verbessertes Eulerverfahren
    case 'bettereuler'
        for k = 1 : N-1
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [u,v] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            xtemp = x1(k)+dt/2*u;
            ytemp = y1(k)+dt/2*v;
            U = u1 + (k-1/2)/(N-1)*(u2-u1);
            V = v1 + (k-1/2)/(N-1)*(v2-v1);
            [u,v] = getValueMonic(U,V,dx,dy,sx,sy,xtemp,ytemp);
            x1(k+1) = x1(k) + dt*u;
            y1(k+1) = y1(k) + dt*v;
        end
        
    % Euler-Heun (einfaches Praediktor-Korrektor)
    case 'eulerheun'
        for k = 1 : N-1
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [uj,vj] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            xtemp = x1(k)+dt*uj;
            ytemp = y1(k)+dt*vj;
            U = u1 + k/(N-1)*(u2-u1);
            V = v1 + k/(N-1)*(v2-v1);
            [uk,vk] = getValueMonic(U,V,dx,dy,sx,sy,xtemp,ytemp);
            x1(k+1) = x1(k) + dt/2*(uj+uk);
            y1(k+1) = y1(k) + dt/2*(vj+vk);
        end
    % Crank-Nicolson (mehrfaches Praediktor-Korrektor)
    case 'cranknicolson'
        for k = 1 : N-1
            maxIt = 10; tol = 1e-13;
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [uj,vj] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            xj = x1(k)+dt*uj;
            yj = y1(k)+dt*vj;
            U = u1 + k/(N-1)*(u2-u1);
            V = v1 + k/(N-1)*(v2-v1);
            for l = 1:maxIt
                [uk,vk] = getValueMonic(U,V,dx,dy,sx,sy,xj,yj);
                xk = x1(k) + dt/2*(uj+uk);
                yk = y1(k) + dt/2*(vj+vk);
                if norm([xk yk]-[xj yk]) < tol
                    break
                end
                xj = xk; yj = yk;
            end
            x1(k+1) = xk;
            y1(k+1) = yk;
        end
    % Runge-Kutta 4. Ordnung (Newtonsche 3/8 Regel)
    case 'rungekutta4newton'
        for k = 1 : N-1
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [k11,k12] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            U = u1 + (k-1+1/3)/(N-1)*(u2-u1);
            V = v1 + (k-1+1/3)/(N-1)*(v2-v1);
            [k21,k22] = getValueMonic(U,V,dx,dy,sx,sy,x1(k)+dt/3*k11,y1(k)+dt/3*k12);
            U = u1 + (k-1+2/3)/(N-1)*(u2-u1);
            V = v1 + (k-1+2/3)/(N-1)*(v2-v1);
            [k31,k32] = getValueMonic(U,V,dx,dy,sx,sy,...
                x1(k)-dt/3*k11+dt*k21,y1(k)-dt/3*k12+dt*k22);
            U = u1 + k/(N-1)*(u2-u1);
            V = v1 + k/(N-1)*(v2-v1);
            [k41,k42] = getValueMonic(U,V,dx,dy,sx,sy,...
                x1(k)+dt*(k11-k21+k31),y1(k)+dt*(k12-k22+k32));
            x1(k+1) = x1(k) + dt/8*(k11+3*(k21+k31)+k41);
            y1(k+1) = y1(k) + dt/8*(k12+3*(k22+k32)+k42);
        end
    % Runge-Kutta 4. Ordnung (klassisch)
    case 'rungekutta4classic'
        for k = 1 : N-1
            U = u1 + (k-1)/(N-1)*(u2-u1);
            V = v1 + (k-1)/(N-1)*(v2-v1);
            [k11,k12] = getValueMonic(U,V,dx,dy,sx,sy,x1(k),y1(k));
            U = u1 + (k-1/2)/(N-1)*(u2-u1);
            V = v1 + (k-1/2)/(N-1)*(v2-v1);
            [k21,k22] = getValueMonic(U,V,dx,dy,sx,sy,x1(k)+dt/2*k11,y1(k)+dt/2*k12);
            [k31,k32] = getValueMonic(U,V,dx,dy,sx,sy,x1(k)+dt/2*k21,y1(k)+dt/2*k22);
            U = u1 + k/(N-1)*(u2-u1);
            V = v1 + k/(N-1)*(v2-v1);
            [k41,k42] = getValueMonic(U,V,dx,dy,sx,sy,x1(k)+dt*k31,y1(k)+dt*k32);
            x1(k+1) = x1(k) + dt/6*(k11+2*(k21+k31)+k41);
            y1(k+1) = y1(k) + dt/6*(k12+2*(k22+k32)+k42);
        end
    otherwise
        error('Method not implemented yet');
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

