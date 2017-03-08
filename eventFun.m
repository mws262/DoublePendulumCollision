function [value,isterminal,direction] = eventFun(t,z,p)

% Pick conditions based on contact configuration
switch p.contact
    %% No contact yet, look for collisions
    case 1
        % Condition 1 -- Perpendicular distance between link 1 and contact
        % pt goes to 0, and within length of link.
        % Condition 2 -- ... link 2
        [col1a,col1b,col2a,col2b] = collisionCheck(p.l1,p.l2,z(1),z(2),p.xcol,p.ycol);
        value(1,1) = col1a; % Perpendicular distance of point to link.
        value(2,1) = col2a;
        isterminal = [col1b;col2b]; % Only make this a terminal event if the contact point is within the swing arc of that link
        direction = [0;0];
        
    %% Contact with link 1 breaking    
    case 2
        % Normal force goes negative (ie contact point can't suck)
        normalForce = c1force(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,z(1),z(2),z(3),z(4),p.xcol,p.ycol);
        value = p.colside*normalForce; % Force will have a different sign depending on the side of the contact point relative to the link
        isterminal = 1;
        direction = -1;
        
    %% Contact with link 2 breaking
    case 3 
        % Condition 1 -- normal force -> 0
        % Condition 2 --  Link tip slides off the contact pt
        % Condition 3 -- Link joint end slides off the contact pt
        normalForce = c2force(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,p.mu,z(1),z(2),z(3),z(4),p.xcol,p.ycol);
        [c2slideofftip,c2slideoffjoint] = c2slideguard(p.l1,p.l2,z(1),z(2),p.xcol,p.ycol);
        value(1,1) = p.colside*normalForce;
        value(2,1) = c2slideofftip;
        value(3,1) = c2slideoffjoint;
        isterminal = [1;1;1];
        direction = [-1;1;-1];
end

end