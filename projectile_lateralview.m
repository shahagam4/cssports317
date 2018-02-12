%% Code for volleyball from lateral view
close all;
clear all;
clc;

%% Declaration of global variable
global mass Cd Cl air_density Wx Wy Wz W Area;

%% Value assignment of global variable for Volleyball example
mass=0.27;                      % mass of the projectile
radius=0.105;                   % radius of projectile
Area=pi*(radius^2);             % area of projectile
COR=0.9;                        % coefficient of restitution when strike with table(for T.T. only)
Cd=0.47;                        % Drag coefficient
air_density=1.225;              % Density of air

%% Array of initial conditions for 3 different cases of Volleyball 
init = [-50.26,-6.28,0,1,0,3.5,10,30,1;
        0,0,0,4.5,5,1.5,1,1.5,9;
        0,0,0,4,8.5,3,18,21.13,-20.6
    ];

%% array for color selection. One can also use rand(1,3) for color selection
colorr=[1,0,0;0,1,0;0,0,1;0,1,1;0,1,0;0.95,0.44,0.25;0.23,0.70,0.68;0,1,0;1,1,0];

%% loop for different cases.
for i = [1,2,3]

% set initial value for different cases
Wx=init(i,1);                   % X component of angular velocity rad/sec
Wy=init(i,2);                   % Y component of angular velocity rad/sec
Wz=init(i,3);                   % Z component of angular velocity rad/sec
W=sqrt(Wx^2+Wy^2+Wz^2);         % Value of angular velocity
if (W==0)                       % If W=0 initialise it with very small value.
    W=1e-10;
end
Cl=0.25;                        % Lift coefficient
x=init(i,4);                    % Initial X-coordinate
y=init(i,5);                    % Initial Y-coordinate
z=init(i,6);                    % Initial Z-coordinate
Vx=init(i,7);                   % Initial value of velocity in X direction
Vy=init(i,8);                   % Initial value of velocity in Y direction
Vz=init(i,9);                   % Initial value of velocity in Z direction



% set initial conditions for ode
answer0(1)=x;
answer0(2)=y;
answer0(3)=z;
answer0(4)=Vx;
answer0(5)=Vy;
answer0(6)=Vz;
answer0(7)=0;
% set options for ode
options=odeset('Events',@projectile_sc);
% Set simulation time
maxtime=5;
delta_t=0.03;                   % step size

% ode simulate projectile motion
if i~=2
[t,answer]=ode45(@projectile_fn,0:delta_t:maxtime,answer0,options);
else
delta_t=0.09;                   % different step size for particular case for better result
[t,answer]=ode45(@projectile_fn,0:delta_t:maxtime,answer0,options);
end

% Drop last point for 1st and 3rd case for better result
if i~=2
answ=answer(1:(end-1),:);       % final answers
tt=t(1:(end-1),:);              % total time    
else
answ=answer;                    % final answers
tt=t;                           % total time
end

plot(answ(:,2),answ(:,3),'k');
hold on;

% Calculating energies
potential_energy_factor=5;      % scaling of potential energy for better result

E=(0.5*mass*(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2)+potential_energy_factor*mass*9.8*answ(:,3))+(1/3)*mass*radius^2*W^2; % Kinectic + Potential

TotalEnergy=(E+answ(:,7)+0.0001);                               % Kinetic + Potential + Lost
KineticAndPotential=(E+0.0001);                                 % Kinetic + Potential
Potential=potential_energy_factor*(mass*9.8*answ(:,3)+0.0001);  % Potential

% inner most circle->potential,  middle ring->kinetic , outer most ring->lost energy
size_factor=150;                 % scaling of circle for better result





% Making array for giving shade for Drag force
Drag=0.5*air_density*Area*Cd*(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2);
Ndrag= 0.7*((Drag-min(Drag))/(max(Drag)-min(Drag)))+0.3;
DragArray=[zeros(length(Ndrag),1) Ndrag.*ones(length(Ndrag),1) zeros(length(Ndrag),1)];


if i==1
% Making array for giving shade for Magnus force
Magnus=0.5*air_density*Area*Cl*(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2);
Nmagnus = 0.7*((Magnus-min(Magnus))/(max(Magnus)-min(Magnus)))+0.3;
MagnusArray = [Nmagnus.*ones(length(Nmagnus),1) zeros(length(Nmagnus),1) zeros(length(Nmagnus),1)];

% Ploting scatter 
scatter(answ(:,2),answ(:,3),(answ(:,1)+2.1)*size_factor,MagnusArray,'filled','MarkerEdgeColor','k'); % improve for size
hold on;
scatter(answ(:,2),answ(:,3),((answ(:,1)+2.1).*(KineticAndPotential./TotalEnergy))*size_factor,DragArray,'filled','MarkerEdgeColor','k');
hold on;
scatter(answ(:,2),answ(:,3),((answ(:,1)+2.1).*(Potential./TotalEnergy))*size_factor,colorr((i-1)*3+3,:),'filled','MarkerEdgeColor','k');
hold on;
else
% Ploting scatter
scatter(answ(:,2),answ(:,3),(answ(:,1)+2.1)*size_factor,colorr((i-1)*3+1,:),'filled','MarkerEdgeColor','k'); % improve for size
hold on;
scatter(answ(:,2),answ(:,3),((answ(:,1)+2.1).*(KineticAndPotential./TotalEnergy))*size_factor,DragArray,'filled','MarkerEdgeColor','k');
hold on;
scatter(answ(:,2),answ(:,3),((answ(:,1)+2.1).*(Potential./TotalEnergy))*size_factor,colorr((i-1)*3+3,:),'filled','MarkerEdgeColor','k');
hold on;
end

% Direction of Magnus force using quiver
if i==1
MDX=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wy*answ(:,6)-Wz*answ(:,5))/W)));
MDY=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wz*answ(:,4)-Wx*answ(:,6))/W)));
MDZ=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wx*answ(:,5)-Wy*answ(:,4))/W)));
q=quiver3(answ(:,2),answ(:,3),answ(:,1),MDY,MDZ,MDX);
q.AutoScaleFactor=0.5;
q.Color=[0.4940 0.1840 0.5560];
end
hold on;

end

% lable the plot
xlabel('y(m)-along length of court');
ylabel('z(m)-along height from ground');
grid on;
set(gca,'FontSize',20);
