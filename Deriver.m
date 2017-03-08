%% Derive equations of motion, collision events, guards, etc for the double pendulum
% Equations of motion are mostly done via angular momentum balance.
% Collision events are angular impulse-momentum
%
% Matthew Sheen
%
clear all;

% Define state variables and parameters
syms th1 th1dot th2 th2dot th1dotdot th2dotdot m1 m2 I1 I2 l1 l2 d1 d2 g real
% Define collision point.
syms x y real

%% Kinematics
i = [1 0 0]';  
j = [0 1 0]';
k = [0 0 1]';

er1 = sin(th1)*i - cos(th1)*j; % Polar frame travelling with link 1.
eth1 = cos(th1)*i + sin(th1)*j;

er2 = sin(th2 + th1)*i - cos(th2 + th1)*j; % Polar frame travelling with link 2.
eth2 = cos(th2 + th1)*i + sin(th2 + th1)*j;

state = [th1,th2,th1dot,th2dot]';
dstate = [th1dot,th2dot,th1dotdot,th2dotdot]';

rj1_o = l1*er1; % Anchor pt to joint 1
rm1_o = d1*er1; % Anchor pt to COM 1

rj2_j1 = l2*er2; % Joint 1 to end of link 2
rm2_j1 = d2*er2; % Joint 1 to COM 2

rm2_o = rm2_j1 + rj1_o; % Anchor pt to COM 2
rj2_o = rj2_j1 + rj1_o; % Anchor pt to end of link 2

vm1_o = jacobian(rm1_o,state)*dstate; % Inertial velocity of COMs
vm2_o = jacobian(rm2_o,state)*dstate;

am1_o = jacobian(vm1_o,state)*dstate; % Inertial acceleration of COMs
am2_o = jacobian(vm2_o,state)*dstate;

colPt = [x,y,0]'; % Collision point

matlabFunction(rj1_o,rj2_o,'File','positions'); % Write position calculations to file.

%% Nominal case -- swinging freely equations of motion
% Angular momentum about elbow pt
MomL2 = dot(I2*(th2dotdot+th1dotdot)*k + m2*cross(rm2_j1,am2_o) - cross(rm2_j1,-m2*g*j),k);
% Angular momentum about o
MomFull = dot((I1*th1dotdot + I2*(th2dotdot+th1dotdot))*k + m1*cross(rm1_o,am1_o) + m2*cross(rm2_o,am2_o) - cross(rm1_o,-m1*g*j) - cross(rm2_o,-m2*g*j),k);

sol = solve(MomL2,MomFull,[th1dotdot,th2dotdot]); % Solve for angular accelerations

matlabFunction(sol.th1dotdot,sol.th2dotdot,'File','accels');

%% Link 1 in contact -- find EOM and guard
syms Fnorm real; % Define a contact normal force
[th1dotdotc1,th2dotdotc1] = solve(subs(MomL2,[th1dot,th1dotdot],[0,0]), [th1dotdot,th2dotdot]);
matlabFunction(th1dotdotc1,th2dotdotc1,'File','accelC1');

colOnLink1 = dot(colPt,er1)*er1; % If we assume colPt is aligned with er1 (which it should be if there's a collision), then this is the distance down the rod.

% Moment rxF
contactMoment = cross(colOnLink1,Fnorm*eth1);
Fcontact1 = solve(subs(MomFull - dot(contactMoment,k),[th1dot,th1dotdot,th2dotdot],[0,0,sol.th2dotdot]),Fnorm);

matlabFunction(Fcontact1,'File','c1force');

%% Link 2 in sliding contact -- find EOM and guards
syms c2force mu real
c2ptVel = vm2_o + cross([0 0 th2dot + th1dot],(colPt - rm2_o))';
c2constraint = jacobian(dot(c2ptVel,eth2),state)*dstate;

% Angular momentum about elbow pt
c2MomL2 = dot(I2*(th2dotdot+th1dotdot)*k + m2*cross(rm2_j1,am2_o) - cross(rm2_j1,-m2*g*j) - cross(colPt - rj1_o,c2force*eth2),k);
% Angular momentum about o
c2fric = -mu*c2force*dot(c2ptVel,er2)*er2; % Sliding friction force
c2MomFull = dot((I1*th1dotdot + I2*(th2dotdot+th1dotdot))*k + m1*cross(rm1_o,am1_o) + m2*cross(rm2_o,am2_o) - cross(rm1_o,-m1*g*j) - cross(rm2_o,-m2*g*j) - cross(colPt,c2force*eth2) - cross(colPt,c2fric),k);

% Solve for the constrained equations of motion and the constraint force
[th1dotdotc2,th2dotdotc2,c2forcesol] = solve(c2MomL2,c2MomFull,c2constraint,[th1dotdot,th2dotdot,c2force]);

matlabFunction(c2forcesol,'File','c2force'); % Constraint force is a guard for when contact breaks.
matlabFunction(th1dotdotc2,th2dotdotc2,'File','accelC2');

% Find when the bar slides off the contact point
c2slideofftip = dot(colPt - rj1_o,er2) - l2; % Link 2 tip loses contact
c2slideoffjoint = dot(colPt - rj1_o,er2); % Link 2 loses contact at the joint end.

matlabFunction(c2slideofftip,c2slideoffjoint,'File','c2slideguard');

%% Energy - mostly for debugging
KE = 1/2*m1*dot(vm1_o,vm1_o) + 1/2*m2*dot(vm2_o,vm2_o) + 1/2*I1*th1dot^2 + 1/2*I2*(th1dot + th2dot)^2;
PE = m1*g*dot(rm1_o,j) + m2*g*dot(rm2_o,j);

matlabFunction(KE+PE,'File','energy');

%% Collision detection
col1a = dot(cross(colPt,er1),k); % Perpendicular distance to collision pt goes to 0 as collision point aligns with link 1 length.
col1b = norm(colPt) <= l1 & dot(colPt,er1) > 0; % Can't be too far out AND must be on the same side as the collision point.

col2a = dot(cross(colPt-rj1_o,er2),k);
col2b = norm(colPt - rj1_o) <= l2 & dot(colPt - rj1_o,er2) > 0;

matlabFunction(col1a,col1b,col2a,col2b,'file','collisionCheck');

%% Collision resolution (impulse-momentum)
% Map velocities before collision to velocities after

% Collision of first link:
syms th1col th1dot_0 th1dot_1 th2dot_0 th2dot_1 rest real
v1_0 = subs(vm1_o,[th1dot,th2dot],[th1dot_0,th2dot_0]); % Define before and after velocities for both links
v1_1 = subs(vm1_o,[th1dot,th2dot],[th1dot_1,th2dot_1]);
v2_0 = subs(vm2_o,[th1dot,th2dot],[th1dot_0,th2dot_0]);
v2_1 = subs(vm2_o,[th1dot,th2dot],[th1dot_1,th2dot_1]);

% Angular momentum conservation about elbow joint, with restitution as the
% second equation
velAfter = solve(I2*(th2dot_0 + th1dot_0) + dot(cross(rm2_j1,m2*v2_0),k) ...
    == I2*(th2dot_1 + th1dot_1) + dot(cross(rm2_j1,m2*v2_1),k), ...
    th1dot_1 == -rest*th1dot_0, th1dot_1,th2dot_1);

matlabFunction(velAfter.th1dot_1,velAfter.th2dot_1,'File','link1col');

% Collision of second link:
syms J real % Define a variable for collision impulse. I couldn't think of a clever point to make angular momentum be conserved.

% Velocity of link 2 at the contact point before and after
vp0 = v2_0 + cross([0 0 th2dot_0 + th1dot_0],colPt - rm2_o)';
vp1 = v2_1 + cross([0 0 th2dot_1 + th1dot_1],colPt - rm2_o)';

% Collision vector is perpendicular to the link
col2eqn1 = dot(vp1,eth2) == -rest*dot(vp0,eth2);

% Angular momentum about j1 for link 2
col2eqn2 = dot(I2*(th2dot_0 + th1dot_0)*k + cross(rm2_j1,m2*v2_0) - (I2*(th2dot_1 + th1dot_1)*k + cross(rm2_j1,m2*v2_1) + cross(colPt-rj1_o,J*eth2)),k);

% Angular momentum about o for both links
col2eqn3 = dot(I1*th1dot_0*k + I2*(th2dot_0 + th1dot_0)*k + cross(rm1_o,m1*v1_0) + cross(rm2_o,m2*v2_0) - ...
    (I1*th1dot_1*k + I2*(th2dot_1 + th1dot_1)*k + cross(rm1_o,m1*v1_1) + cross(rm2_o,m2*v2_1) + cross(colPt,J*eth2)),k);

col2sol = solve(col2eqn1,col2eqn2,col2eqn3,th1dot_1,th2dot_1,J); % Solve this system

matlabFunction(col2sol.th1dot_1,col2sol.th2dot_1,'File','link2col'); 

