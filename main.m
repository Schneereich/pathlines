clc
close all

x0 = 0;
y0 = 1; %y-Boundary = [-2,2]
N = 500;
T = 10;
Z = 2;
U = rand(5); %muss cell-Array sein, für jeden Zeitpunkt müssen Werte vorhanden sein
V = rand(5); % ||
X = [0, 1, 2, 3, 4; ...
    0, 1, 2, 3, 4; ...
    0, 1, 2, 3, 4; ...
    0, 1, 2, 3, 4; ...
    0, 1, 2, 3, 4];
Y = [-2, -2, -2, -2, -2; ...
    -1, -1, -1, -1, -1; ...
    0, 0, 0, 0, 0; ...
    1, 1, 1, 1, 1; ...
    2, 2, 2, 2, 2];

figure
plot(X(:,:),Y(:,:),'b')
hold on
grid on
[x,y] = bahn(x0,y0,T,N,Z,U,V,X,Y)
plot(x(:),y(:),'r*')
hold on


