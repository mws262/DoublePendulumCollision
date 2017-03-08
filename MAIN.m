%% Double pendulum with a single collision pt
% Hybrid system: continuous integration with discrete collision events
%
% Integration switches between 3 equations of motion depending on the
% contact configuration:
% 1. Free motion
% 2. Collision with link 1 -- pt contact
% 3. Collision with link 2 -- sliding contact with friction
%
% Contact is broken when normal force goes negative OR link 2 slides off
% the contact pt.
% Collisions can be either plastic or elastic.
% Only 1 contact pt at a time
%
% Deriver.m symbolically derives the equations of motion, the contact
% conditions, collision resolutions, etc, and automatically writes them to
% files. MAIN.m, Deriver, RHS.m, and eventFun.m are the only hand-written
% functions.
%
% Matthew Sheen, 2017
clear all; close all;

%% Physical parameters
% Contact parameters
p.xcol = 0.0; % Obstacle location
p.ycol = -1.8;

p.mu = 0.2; % Friction coefficient for the contact
restitution = 0.0; % Bounce. 0 is inelastic, >0 has bounce.
minBounceVel = 1e-2; % How small of a rebound velocity to ignore and just switch to contact.

% Pendulum parameters
p.d1 = 0.5; % Distance from previous joint to COM
p.d2 = 0.5;
p.l1 = 1; % Each link length
p.l2 = 1;
p.m1 = 1; % Link mass
p.m2 = 1;
p.I1 = 1/12*p.m1*p.l1^2; % Link moment of inertia
p.I2 = 1/12*p.m2*p.l2^2;
p.g = 10; % Gravity

%% Initial conditions
% Total simulation time
tend = 10;
% Initial state: [th1, th2, thdot1, thdot2]
inits = [pi/2,0,0,0]';

%% Settings
% Animation
animate = true; % Animation on/off
tfactor = 0.5; % 1 is real time, <1 is slower, >1 is sped up

% Contact-breaking "push"
contactBreakJump = 1e-8; % When contact breaks, kick the position away from the wall by a small amount to prevent it from getting stuck.

% Integration
opts = odeset;
opts.RelTol = 1e-10;
opts.AbsTol = 1e-10;
opts.Events = @(t,z)eventFun(t,z,p); % Handles all collision and breaking contact events

% Create storage variables
ztotal = []; % Keep concatenating with every bounce.
ttotal = [];
phasetotal = 1;
phasetimes = 0;
zcurrent = inits;
tcurrent = 0;
% Initial contact state -- unless you start in contact, you don't need to
% change these
p.contact = 1; % Initially free to move. 2 is link 1 in contact, 3 is link 2 in contact.
p.colside = 1; % 1 is the positive angle side, -1 is the negative angle side.

%% Integrate the hybrid dynamics. Each time contact changes, switch  equations and keep going!
while tcurrent < tend
    opts.Events = @(t,z)eventFun(t,z,p); % Refreshes all the values in p. (Mostly contact mode).
    [tarray,zarray,te,ze,ie] = ode45(@(t,z)odeFun(t,z,p),[tcurrent,tend],zcurrent,opts);
    tcurrent = tarray(end);
    zcurrent = zarray(end,:);
    
    %% Resolve collisions, other events
    % Contact modes:
    % 1 -- no contact
    % 2 -- link 1 touches
    % 3 -- link 2 touches
    if ~isempty(ie)
        %%% No contact %%%
        if p.contact == 1
            %- Event 1: Link 1 hits -%
            if ie(end) == 1
                [newVel1,newVel2] = link1col(p.I2,p.d2,p.l1,p.m2,restitution,zarray(end,1),zarray(end,2),zarray(end,3),zarray(end,4));
                [rj1_o,rj2_o] = positions(p.l1,p.l2,zarray(end-1,1),zarray(end-1,2));
                p.colside = sign(dot([0 0 1],cross([p.xcol;p.ycol;0],rj1_o)));
                
                if ~restitution % Only stay in contact if there's no bounce
                    p.contact = 2;
                elseif newVel1 < minBounceVel % If there's restitution, but the bounce is too small, recalculate with no restitution
                    [newVel1,newVel2] = link1col(p.I2,p.d2,p.l1,p.m2,0,zarray(end,1),zarray(end,2),zarray(end,3),zarray(end,4));
                    p.contact = 2;
                else
                    zcurrent(1) = zcurrent(1) + p.colside*contactBreakJump; % Kick it off the wall a tiny bit if it bounces. Keeps the events from getting screwy.
                end
                zcurrent(3:4) = [newVel1,newVel2];
                %- Event 2: Link 2 hits -%
            elseif ie(end) == 2
                [newVel1,newVel2] = link2col(p.I1,p.I2,p.d1,p.d2,p.l1,p.m1,p.m2,restitution,zarray(end,1),zarray(end,2),zarray(end,3),zarray(end,4),p.xcol,p.ycol);
                
                [rj1_o,rj2_o] = positions(p.l1,p.l2,zarray(end-1,1),zarray(end-1,2));
                p.colside = sign(dot([0 0 1],cross([p.xcol;p.ycol;0] - rj1_o,rj2_o - rj1_o)));;
                if ~restitution % Only stay in contact if there's no bounce.
                    p.contact = 3;
                elseif newVel2 < minBounceVel % If there's restitution, but the bounce is too small, recalculate with no restitution
                    [newVel1,newVel2] = link2col(p.I1,p.I2,p.d1,p.d2,p.l1,p.m1,p.m2,0,zarray(end,1),zarray(end,2),zarray(end,3),zarray(end,4),p.xcol,p.ycol);
                    p.contact = 3;
                else
                    zcurrent(2) = zcurrent(2) + p.colside*contactBreakJump; % Kick it off the wall a tiny bit if it bounces. Keeps the events from getting screwy.
                end
                zcurrent(3:4) = [newVel1,newVel2];
            end
            
            %%% Link 1 in contact %%%
        elseif p.contact == 2 && ie(end) == 1
            %- Event 1: Link 1 normal force -> 0 -%
            p.contact = 1;
            zcurrent(1) = zcurrent(1) + p.colside*contactBreakJump;
            %%% Link 2 in contact %%%
        elseif p.contact == 3
            %- Event 1: Link 2 normal force -> 0 -%
            if ie(end) == 1
                p.contact = 1;
                zcurrent(2) = zcurrent(2) + p.colside*contactBreakJump;
                %- Event2: Link 2 slides off contact pt -%
            else
                p.contact = 1;
            end
        end
    end
    
    % Concatenate the newest continuous phase of the dynamics with the
    % others
    ztotal = [ztotal; zarray(1:end-1,:)];
    ttotal = [ttotal; tarray(1:end-1)];
    phasetotal = [phasetotal;p.contact];
    phasetimes = [phasetimes;tcurrent];
    
end
ztotal = [ztotal; zarray(end,:)]; % Get that last point
ttotal = [ttotal; tarray(end)];

% Make an energy plot
en = energy(p.I1,p.I2,p.d1,p.d2,p.g,p.l1,p.m1,p.m2,ztotal(:,1),ztotal(:,2),ztotal(:,3),ztotal(:,4));
figEn = figure;
plot(ttotal,en);
title('Total energy');
xlabel('Time (s)');
ylabel('Energy (J)');
% plot([phasetimes,phasetimes]',repmat([-1; 1],1,length(phasetimes)));


%% Animation
fig = figure;
rod1 = plot([0,p.l1],[0,p.l1],'b','LineWidth',10);
hold on;
rod2 = plot([0,p.l2],[0,p.l2],'r','LineWidth',10);
plot(p.xcol,p.ycol,'.','MarkerSize',30);
hold off;
axis([-2 2 -2 2]);
daspect([1,1,1]);

currTime = 0;
tic;
while currTime*tfactor < ttotal(end) && ishandle(fig) && animate
    
    currState = interp1(ttotal,ztotal,currTime*tfactor);
    
    rot1 = [cos(currState(1)),-sin(currState(1)); sin(currState(1)), cos(currState(1))];
    rot2 = [cos(currState(1) + currState(2)),-sin(currState(1) + currState(2)); sin(currState(1) + currState(2)), cos(currState(1) + currState(2))];
    
    rod1new = rot1*[0,0;0,-p.l1];
    rod2new = rot2*[0,0;0,-p.l2] + repmat(rod1new(:,2),[1,2]);
    
    rod1.XData = rod1new(1,:);
    rod1.YData = rod1new(2,:);
    
    rod2.XData = rod2new(1,:);
    rod2.YData = rod2new(2,:);
    
    drawnow;
    currTime = toc;
end
