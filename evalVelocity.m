function [uk,vk] = evalVelocity(X,Y,U,V,xk,yk)

% evtl. global definieren
m = size(X,2);
n = size(X,1);
dx = X(1,end)/m;
dy = Y(end,1)/n;

%% find element
leftIdx = floor(xk/dx)+1; %Indizes
rightIdx = ceil(xk/dx)+1;
downIdx = floor(yk/dy)+1;
upIdx = ceil(yk/dy)+1;

% Fallunterscheidung
if leftIdx < 1 || rightIdx > m || downIdx < 1 || upIdx > n %auﬂerhalb des Gitters
    warning("Out of boundaries.");
    uk = NAN;
    vk = NAN;
    return
end

if any(isnan([U(leftIdx,downIdx),U(leftIdx,upIdx),U(rightIdx,downIdx),U(rightIdx,upIdx)])) %im Loch
    warning("Velocity is infinity.");
    uk = inf;
    vk = inf;
    return
end
    
alpha = (xk-X(1,leftIdx))/dx;
beta = (yk-Y(downIdx,1))/dy;


%% interpolation
u1 = U(leftIdx,downIdx)+alpha*(U(rightIdx,downIdx)-U(leftIdx,downIdx));
u2 = U(leftIdx,upIdx)+alpha*(U(rightIdx,upIdx)-U(leftIdx,upIdx));
v1 = V(leftIdx,downIdx)+alpha*(V(rightIdx,downIdx)-V(leftIdx,downIdx));
v2 = V(leftIdx,upIdx)+alpha*(V(rightIdx,upIdx)-V(leftIdx,upIdx));


uk = u1+beta*(u2-u1);
vk = v1+beta*(v2-v1);

end

