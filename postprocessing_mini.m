function varargout = postprocessing_mini(varargin)
%*** varargout{1} = timerVal (runtime of pathline iteration)

%*** parse additional arguments
p = inputParser;
addParameter(p,'N',[])
addParameter(p,'subSteps',100)
addParameter(p,'first',100)
addParameter(p,'method',"explicitEuler")
parse(p,varargin{:})

%*** load solution
load solution_iNS_mini_2.mat N coordinates dirichlet neumann elements3 T U P

%*** if defined, overwrite default N
if p.Results.N
    N = p.Results.N;
end
fprintf("Calculating with method '%s' and %d substeps ...\n", ...
        p.Results.method, p.Results.subSteps)

t = linspace(0,T,N);
ifront = abs(coordinates(:,1)-0.15)<1e-5 ...
       & abs(coordinates(:,2)-0.2 )<1e-5;
iback  = abs(coordinates(:,1)-0.25)<1e-5 ...
       & abs(coordinates(:,2)-0.2 )<1e-5;

pdiff = zeros(1,N);
for k=1:N
    pdiff(k) = -(P{k}(ifront) - P{k}(iback));
end

figure(1)
plot(t(30:end),pdiff(30:end))
title('Pressure difference')

figure(2),clf
plot(reshape(coordinates(dirichlet,1),[],2)', ...
     reshape(coordinates(dirichlet,2),[],2)','r-','linewidth',2)
hold on
plot(reshape(coordinates(neumann,1),[],2)', ...
     reshape(coordinates(neumann,2),[],2)','g-','linewidth',2)

%*** plot pathline for one single point first
x0 = [0,0,0,0,0,0,0,0,0];
y0 = [0.045, 0.085,0.125,0.165,0.205,0.245,0.285,0.3250,0.365];
%*** initiate mesh to map solution to
s = linspace(min(coordinates(:,1)),max(coordinates(:,1)),300);
t = linspace(min(coordinates(:,2)),max(coordinates(:,2)),160);
[x,y] = meshgrid(s,t);

%*** start with 1 or 100
first = p.Results.first;
last =  N-1;
%*** initiate the first solution
u = U{first};
ux = reshape(u(elements3,1),[],3);
uy = reshape(u(elements3,2),[],3);
uv1 = tri2monic(coordinates,elements3,{ux,uy},x,y);

%*** initiate array for timing
timerVal = zeros(length(x0),length(first:last));

for k = first:last
    k
    u = U{k+1};
    ux = reshape(u(elements3,1),[],3);
    uy = reshape(u(elements3,2),[],3);
    uv2 = tri2monic(coordinates,elements3,{ux,uy},x,y);
    for j=1:length(x0)
        tic
        [vx,vy] = pathline(uv1{1},uv1{2},uv2{1},uv2{2},s,t,x0(j),y0(j), ...
                           T/(N-1),p.Results.subSteps,p.Results.method);
        timerVal(j,k) = toc;
        x0(j) = vx(end);
        y0(j) = vy(end);
        plot(vx,vy,'k.')
    end
    drawnow
    %*** save uv2 for the next iteration
    uv1 = uv2;
end

%*** return the mean calculation time of all pathlines for each timestep
if nargout == 1
    varargout{1} = mean(timerVal);
end

hold off
axis equal