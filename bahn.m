function [x,y] = bahn(x0,y0,T,N,Z,U,V,X,Y)

dt = T/N; %Zeitschrittweite
h = dt/(Z+1); %zus�tzliche Punkte zwischen den Zeitschritten
x = zeros(1,N+Z(N-1));
y = zeros(1,N+Z(N-1));

x(1) = x0;
y(1) = y0;

% Zeitpunke, an denen die Position des Partikels bestimmt wird
t = linspace(0,T,N+Z(N-1));

for k = 1:length(t)-1

    tp = floor(t(k)/dt)+1;
    tn = ceil(t(k)/dt)+1;
    gamma = (t(k)-tp)/dt;
    
    if gamma ~= 0 % zwischen den Knotenpunkten
        
        [up,vp] = evalVelocity(X,Y,U{tp},V{tp},x(k),y(k));
        [un,vn] = evalVelocity(X,Y,U{tn},V{tn},x(k),y(k)); % Verfahrensfehler durch Approx.
        uk = up + gamma*(un-up);
        vk = vp + gamma*(vn-vp);
        
    else % auf Knotenpunkt
        
       [uk,vk] = evalVelocity(X,Y,U{tp},V{tp},x(k),y(k));

    end
    
    if isnan(uk)
        if y(k) > Y(end,1) %�ber Ymax Grenze dr�ber geschossen
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            x(k) = x(k-1) + (Y(end,1)-y(k-1))/s;
            y(k) = Y(end,1);
        elseif x(k) > X(1,end) %�ber Xmax Grenze dr�ber geschossen
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            y(k) = s*(X(1,end)+x(k-1)) + y(k-1);
            x(k) = X(1,end);
        else %�ber Ymin Grenze dr�ber geschossen
            s = (y(k)-y(k-1))/(x(k)-x(k-1));
            x(k) = x(k-1) + (Y(1,1)-y(k-1))/s;
            y(k) = Y(1,1);          
        end
        % x, y Vektor - Reste abschneiden
        x = x(1:k);
        y = y(1:k);
        return
    end
    
    if isinf(uk)
        % x, y Vektor - Reste abschneiden
        x = x(1:k);
        y = y(1:k);
        return
    end
    
    % Euler
    x(k+1) = x(k) + h * uk;
    y(k+1) = y(k) + h * vk;
    

end
end


