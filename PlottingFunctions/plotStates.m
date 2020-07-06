function plotStates(GPOPSoutput)


% 'interpsolution' are data interpolated between collocation points
% 'solution' are the GPOPS output at collocation points
solType = {'interpsolution','solution'};
markerType = {'-','o'};

% initialize figure
f1 = figure('color','w');

% Specify legend fontsize
lgdFS = 14;


P = f1.Position;
% make the figure double the default size
P = [P(1),round(P(2)/2),P(3)*2,P(4)*2];
f1.Position = P;

sp = NaN(6,1);
for i = 1:6
    sp(i) = subplot(3,2,i);
end

% Extract input parameters

aux = GPOPSoutput.result.setup.auxdata;
D = aux.D;
T = aux.T;
g = aux.g;

for i = 1:2
    % extract data
    t = GPOPSoutput.result.(solType{i}).phase.time;
    X = GPOPSoutput.result.(solType{i}).phase.state;

    %%% Plot cartesian positions in time
    axes(sp(1)) %#ok<*LAXES>
    resetcolor % resets colors to default order
    hold on
    box on
    
    x = X(:,1);
    y = X(:,2);
    
    % subtract expected x position based on matching mean velocity
    x_meanvel = D/T*t;
    x_corr = x - x_meanvel;
    
    % plot x position
    yyaxis left
    plot(t,x_corr,markerType{i})
    
    % plot y position
    yyaxis right
    plot(t,y,markerType{i})
    
    if i == 1
    % plot labels
    ylabel('Vertical COM position [m]')
    
    yyaxis left
    ylabel('Horz. COM pos. - $Dt/T$ [m]','interpreter','latex')
    
    end
    
    %%% Plot cartesian velocities in time
    axes(sp(3))
    resetcolor
    hold on
    box on
    
    u = X(:,3); % horizontal velocity
    v = X(:,4); % vertical velocity
    
    % subtract expected x velocity based on matching mean velocity
    u_corr = u - D/T;
    
    % plot horizontal velocity
    plot(t,u_corr,markerType{i})
    
    % plot vertical velocity
    plot(t,v,markerType{i})
    
    if i == 1
        % plot legend
        legend({'Horz. vel. $- D/T$','Vert. vel.'},...
            'interpreter','latex','autoupdate','off','fontsize',lgdFS)
    end
    %%% Plot rotation and rotaional velocities
    axes(sp(5))
    resetcolor
    hold on
    box on
    
    o = X(:,5); % angular position, from horizontal
    w = X(:,6); % angular velocity
    
    % convert to degress
    o_deg = o*180/pi;
    w_deg = w*180/pi;
    
    % plot angular position
    yyaxis left
    plot(t,o_deg,markerType{i})    
    
    % plot angular velocity
    yyaxis right
    plot(t,w_deg,markerType{i})
    
    if i == 1
        % plot y labels
        ylabel('Ang. velocity [deg/s]')
        
        yyaxis left
        ylabel('Ang. position [deg]')
        
        % plot x axis
        xlabel('Time (s)')
    end
    
    %%% Plot forces
    axes(sp(2))
    resetcolor
    hold on
    box on
    
    F = X(:,7:9);
    
    % Divide by gravitational acceleration
    F_norm = F/g;
    plot(t,F_norm,markerType{i})
    if i == 1
        % Set axis lables and legends
        legend({'$F_{tr}$','$F_{ld}$','$F_{ref}$'},...
            'interpreter','latex','AutoUpdate','off','fontsize',lgdFS)
        ylabel('Axial Force [$mg$]','interpreter','latex')
    end
    
    %%% Plot Torques
    axes(sp(4))
    resetcolor
    hold on
    box on
    
    Tau = X(:,10:12);
    
    plot(t,Tau,markerType{i})
    
    if i == 1
        % plot legend
        legend({'$\tau_{tr}$','$\tau_{ld}$','$\tau_{ref}$'},...
            'interpreter','latex','AutoUpdate','off','fontsize',lgdFS)
    end
    
    
    %%% Plot Impulse and impulse-ish states
    axes(sp(6))
    resetcolor
    hold on
    box on
    
    P = X(:,13);
    Q2 = X(:,14);

    % Plot F_ld impulse
    plot(t,P,markerType{i})
    
    % Plot Tau_ld^2 impulse-ish thing
    plot(t,Q2,markerType{i})
    
    if i == 1
        % Plot legend
        legend({'$\int_0^t F_{ld}(\zeta) d\zeta$','$\int_0^t \tau_{ld}(\zeta) d\zeta$'},...
            'interpreter','latex','AutoUpdate','off','fontsize',lgdFS,...
            'location','northwest')
        % Plot x label
        xlabel('Time (s)')
    end
end

