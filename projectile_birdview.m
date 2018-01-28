%% Code for volleyball from lateral view
close all;
clear all;
clc;

%% Declaration of global variable
global mass Cd Cl air_density Wx Wy Wz W Area;

%% Value assignment of global variable for Volleyball example
mass=0.27;                   % mass of the projectile
radius=0.105;                % radius of projectile
Area=pi*(radius^2);          % area of projectile
Cd=0.47;                     % Drag coefficient
air_density=1.225;           % Density of air

%% set initial value
Wx=-50.26;                   % X component of angular velocity rad/sec
Wy=-6.28;                    % Y component of angular velocity rad/sec
Wz=0;                        % Z component of angular velocity rad/sec
W=sqrt(Wx^2+Wy^2+Wz^2);      % Value of angular velocity
if (W==0)                    % If W=0 initialise it with very small value.
    W=1e-10;
end
Cl=0.25;                     % Lift coefficient
x=1;                         % Initial X-coordinate
y=0;                         % Initial Y-coordinate
z=3.5;                       % Initial Z-coordinate
Vx=10;                       % Initial value of velocity in X direction
Vy=30;                       % Initial value of velocity in Y direction
Vz=1;                        % Initial value of velocity in Z direction



%% set initial conditions for ode
answer0(1)=x;
answer0(2)=y;
answer0(3)=z;
answer0(4)=Vx;
answer0(5)=Vy;
answer0(6)=Vz;
answer0(7)=0;
%% set options for ode
options=odeset('Events',@projectile_sc);

%% Set simulation time
maxtime=5;
delta_t=0.03;                % step size

%% ode simulate projectile motion
[t,answer]=ode45(@projectile_fn,0:delta_t:maxtime,answer0,options);


%% Drop last point for better result
answ=answer(1:(end-1),:);    % final answers
tt=t(1:(end-1),:);           % total time    

%% Calculating energies
potential_energy_factor=5;      % scaling of potential energy for better result

E=(0.5*mass*(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2)+potential_energy_factor*mass*9.8*answ(:,3))+(1/3)*mass*radius^2*W^2; % Kinectic + Potential

TotalEnergy=(E+answ(:,7)+0.0001);                               % Kinetic Energy + Potential Energy + Lost Energy
KineticAndPotential=(E+0.0001);                                 % Kinetic Energy + Potential Energy
Potential=potential_energy_factor*(mass*9.8*answ(:,3)+0.0001);  % Potential Energy

%% Plotting-1
% inner most circle->potential,  middle ring->kinetic , outer most ring->lost energy
size_factor=115;                 % scaling of circle for better result

scatter(answ(:,2),answ(:,1),(answ(:,3)+2.1)*size_factor,'filled','MarkerEdgeColor','k');
hold on;
scatter(answ(:,2),answ(:,1),((answ(:,3)+2.1).*(KineticAndPotential./TotalEnergy))*size_factor,'filled','MarkerEdgeColor','k');
hold on;
scatter(answ(:,2),answ(:,1),((answ(:,3)+2.1).*(Potential./TotalEnergy))*size_factor,'filled','MarkerEdgeColor','k');
hold on;

% Direction of Magnus force using quiver
MDX=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wy*answ(:,6)-Wz*answ(:,5))/W)));
MDY=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wz*answ(:,4)-Wx*answ(:,6))/W)));
MDZ=0.5*air_density*Area*(sqrt(answ(:,4).^2+answ(:,5).^2+answ(:,6).^2).*(Cl*((Wx*answ(:,5)-Wy*answ(:,4))/W)));
q=quiver3(answ(:,2),answ(:,1),answ(:,3),MDY,MDX,MDZ);
q.AutoScaleFactor=3;
hold on;

%% Ploting-2
plot(answ(:,2),answ(:,1));   % ploting simple line plot 
hold on;

%% lable the plot
xlabel('y(m)-along length');
ylabel('x(m)-along width');
grid on;
set(gca,'FontSize',20);
