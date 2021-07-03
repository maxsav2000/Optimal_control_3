%% Практикум задание 3
%% тест функции reachset
alpha =	3;
t = 1;
hold on
[X, Y,lineX,lineY] = reachset(alpha,t);
plot(X,Y,'LineWidth',2,'color','r');
plot(lineX,lineY,'b.');
legend('Множество достижимости','Кривая переключения')
xlabel('x_1');
ylabel('x_2');
hold off
%% тест функции reachset
alpha =	10;
t = 0.3;
[X, Y] = reachset(alpha,t);
plot(X,Y,'LineWidth',2,'color','r');
legend('Множество достижимости')
xlabel('x_1');
ylabel('x_2');
%% тест функции reachsetdyn (видео abc.avi)
alpha =	1;
t1 = 1;
t2 = 3;
N =	7;
namefile = 'abc.avi';
mov = reachsetdyn(alpha, t1, t2, N, namefile);

%% Запуск видео на котором изменяется множество достижимости
movie(mov);

%%
[X, Y] = reachset(alpha,t);
plot(X,Y,'LineWidth',2,'color','r');
legend('Множество достижимости')
xlabel('x_1');
ylabel('x_2');
function mov = reachsetdyn(alpha, t1, t2, N, filename)
    mov(1:N) = struct('cdata', [],'colormap', []);
    hold on

    t = linspace(t1, t2, N);
    for i=1:N
        [X, Y] = reachset(alpha,t(i));
        plot(X,Y,'LineWidth',2,'color','r');
        mov(i) = getframe();
    end
    legend('Множество достижимости')
    xlabel('x_1');
    ylabel('x_2');
    v = VideoWriter(filename);
    v.FrameRate = 1;
    open(v);
    writeVideo(v,mov);
    close(v);
end


function [X, Y, lineX,lineY] = reachset(alpha, t)
    % S+
    lineX = [];
    lineY = [];
    X = [];
    Y = [];
    tspan1 = [0 t];
    x0 = [0, 0];
    
    opts1 = odeset('Events',@StopEvents_plus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
    [t1,y1] = ode45(@(t,y) odefcn_plus(t,y,alpha), tspan1, x0, opts1);
    %plot(y1(:,1),y1(:,2),'g');
    
    X = [X y1(:,1)'];
    Y = [Y y1(:,2)'];
    lineX = [lineX y1(:,1)'];
    lineY = [lineY y1(:,2)'];
    n = size(t1);
    n = n(1);
    res_plus_1 = zeros(1,n);
    res_plus_2 = zeros(1,n);
    for i =2:n-1
        x0 = [y1(i,1) y1(i,2) 1 0];
        tspan2 = [t1(i) t];
        opts2 = odeset('Events',@StopEvents_psi_minus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
        [t2,y2] = ode45(@(t,y) odefcn_psi_minus(t,y,alpha), tspan2, x0, opts2);
        %plot(y2(:,1),y2(:,2),'k');
        while (1)
            m = length(t2);
            if (t2(m)+0.01<t)
                lineX = [lineX y2(m,1)'];
                lineY = [lineY y2(m,2)'];
                x0 = [y2(m,1) y2(m,2) -1 0];
                tspan2 = [t2(m) t];
                opts3 = odeset('Events',@StopEvents_psi_plus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
                [t2,y2] = ode45(@(t,y) odefcn_psi_plus(t,y,alpha), tspan2, x0, opts3);
                %plot(y2(:,1),y2(:,2),'k');
                X = [X y2(:,1)'];
                Y = [Y y2(:,2)'];
            else
                break;
            end
            m = length(t2);
            if (t2(m)+0.01<t)
                lineX = [lineX y2(m,1)'];
                lineY = [lineY y2(m,2)'];
                x0 = [y2(m,1) y2(m,2) 1 0];
                tspan2 = [t2(m) t];
                [t2,y2] = ode45(@(t,y) odefcn_psi_minus(t,y,alpha), tspan2, x0, opts2);
                %plot(y2(:,1),y2(:,2),'k');
                X = [X y2(:,1)'];
                Y = [Y y2(:,2)'];
            else
                break;
            end
        end
        res_plus_1(i) = y2(length(t2),1);
        res_plus_2(i) = y2(length(t2),2);

    end
    X = [X, res_plus_1];
    Y = [Y, res_plus_2];
    %plot(res_plus_1,res_plus_2,'r.');
    
    
    
    
    % S-
    tspan1 = [0 t];
    x0 = [0, 0];
    opts1 = odeset('Events',@StopEvents_minus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
    [t1, y1] = ode45(@(t,y) odefcn_minus(t,y,alpha), tspan1, x0, opts1);
    %plot(y1(:,1),y1(:,2),'g');
    
    X = [X y1(:,1)'];
    Y = [Y y1(:,2)'];
    lineX = [lineX y1(:,1)'];
    lineY = [lineY y1(:,2)'];
    
    n = size(t1);
    n = n(1);
    res_plus_1 = zeros(1,n);
    res_plus_2 = zeros(1,n);
    for i =2:n-1
        x0 = [y1(i,1) y1(i,2) -1 0];
        tspan2 = [t1(i) t];
        opts2 = odeset('Events',@StopEvents_psi_plus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
        [t2,y2] = ode45(@(t,y) odefcn_psi_plus(t,y,alpha), tspan2, x0, opts2);
        %plot(y2(:,1),y2(:,2),'k');
        while (1)
            m = length(t2);
            if (t2(m)+0.01<t)
                lineX = [lineX y2(m,1)'];
                lineY = [lineY y2(m,2)'];
                x0 = [y2(m,1) y2(m,2) 1 0];
                tspan2 = [t2(m) t];
                opts3 = odeset('Events',@StopEvents_psi_minus,'RelTol',1e-4,'AbsTol',1e-4,'MaxStep',1e-2);
                [t2,y2] = ode45(@(t,y) odefcn_psi_minus(t,y,alpha), tspan2, x0, opts3);
                %plot(y2(:,1),y2(:,2),'k');
                X = [X y2(:,1)'];
                Y = [Y y2(:,2)'];
            else
                break;
            end
            m = length(t2);
            if (t2(m)+0.01<t)
                lineX = [lineX y2(m,1)'];
                lineY = [lineY y2(m,2)'];
                x0 = [y2(m,1) y2(m,2) -1 0];
                tspan2 = [t2(m) t];
                [t2,y2] = ode45(@(t,y) odefcn_psi_plus(t,y,alpha), tspan2, x0, opts2);
                %plot(y2(:,1),y2(:,2),'k');
                X = [X y2(:,1)'];
                Y = [Y y2(:,2)'];
            else
                break;
            end
        end
        res_plus_1(i) = y2(length(t2),1);
        res_plus_2(i) = y2(length(t2),2);

    end
    X = [X res_plus_1];
    Y = [Y res_plus_2];
    %plot(res_plus_1,res_plus_2,'b.');
    if (t<1.5)
    	k = boundary(X',Y',0.05);
    else
        k = boundary(X',Y',0.9);
    end
    %plot(X(k),Y(k));
    X = X(k);
    Y = Y(k);
    
end




function dydt = odefcn_plus(t,y,alpha)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = alpha-y(2)-(sin(2*y(1))-y(1).^3)*cos(2*cos(y(2)));
end
function dydt = odefcn_psi_minus(t,y,alpha)
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = -alpha-y(2)-(sin(2*y(1))-y(1).^3)*cos(2*cos(y(2)));
    dydt(3) = 2*y(4)* cos(y(1))+3*y(4)*(y(1)^2)*cos(2*y(2));
    dydt(4) = -y(3)+y(4)-2*y(4)*(y(1)^3)*sin(2*y(2));
end
function dydt = odefcn_psi_plus(t,y,alpha)
    dydt = zeros(4,1);
    dydt(1) = y(2);
    dydt(2) = alpha-y(2)-(sin(2*y(1))-y(1).^3)*cos(2*cos(y(2)));
    dydt(3) = 2*y(4)* cos(y(1))+3*y(4)*(y(1)^2)*cos(2*y(2));
    dydt(4) = -y(3)+y(4)-2*y(4)*(y(1)^3)*sin(2*y(2));
end
function dydt = odefcn_minus(t,y,alpha)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = -alpha-y(2)-(sin(2*y(1))-y(1).^3)*cos(2*cos(y(2)));
end
function [value,isterminal,direction] = StopEvents_plus(t,y)
value = 0+(y(2)<=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end
function [value,isterminal,direction] = StopEvents_minus(t,y)
value = 0+(y(2)>=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end
function [value,isterminal,direction] = StopEvents_psi_minus(t,y)
value = 0+(y(4)<=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end
function [value,isterminal,direction] = StopEvents_psi_plus(t,y)
value = 0+(y(4)>=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end