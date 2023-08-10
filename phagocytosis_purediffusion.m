% NSCM005 CA1
% modelling a pure diffusion model of phagocytosis for a spherical particle
% clear all variables
clear

% Define parameters
DT = 0.0001; % time step, [s]
DR = 0.02; % lattice step, [um]
T_MAX = 30; % max time, [s]
L = 50; % domain length, [um]
D = 1; % diffusion constant, [um^2/s]
p0 = 50; % initial receptor density, [um^-2]
pL = 500; % bound receptor density, [um^-2]
E = 15; % binding energy per receptor-ligand bond
B = 20; % bending modulus
R = 16; % target radius (radius of bead), [um]

% Declare/allocate variables
num_steps = round(T_MAX/DT); % number of time steps [-]
num_latt_pts = round(L/DR); % number of lattice points [-]

p = zeros(num_steps,num_latt_pts); % initialise receptor density array
a = zeros(num_steps, 1); % cup size

p_plus = fzero(@(p_plus)p_plus/pL - log(p_plus/pL)-E+2*B/(pL*R^2)-1,...
    [1e-6 5]); % calculate p_plus from Eq 4

% Initial conditions
for i = 1:num_latt_pts
    p(1,i) = p0; % initialise receptor densities to p0
end
a(1) = 0; % initialise cup size to 0

% Main loop
for i = 2:num_steps
    %cup size in lattice steps
    a_latt_pt = round( a(i-1)/DR + 1 );

    %set p to pL at r<a
    for j=1:a_latt_pt-1
        p(i,j) = pL;
    end

    % set p to p_plus at r=a
    p(i, a_latt_pt) = p_plus;

    %calculate p for r>a (derived using the Euler method on Eq 2a)
    for j = a_latt_pt+1:num_latt_pts-1
        p(i,j) = p(i-1,j) + D*DT/DR * ( 1/(j*DR)*(p(i-1,j+1)-p(i-1,j)) + ...
            1/DR*(p(i-1,j-1)-2*p(i-1,j)+p(i-1,j+1)));
    end

    %no flux boundary condition at x=L
    p(i, num_latt_pts) = p(i-1,num_latt_pts) + D*DT/DR^2 * (p(i-1,num_latt_pts-1)-p(i-1,num_latt_pts));

    %calculate p'_plus (derivative dp/dr evaluated at r=a)
    p_dash_plus = (p(i,a_latt_pt+1)-p(i,a_latt_pt)) /DR;

    %calculate da/dt from Eq 2b
    a_dot = D*p_dash_plus/(pL-p_plus);

    %recalculate cup size using da/dt
    a(i) = a(i-1) + a_dot*DT;
end



%neither fzero nor fsolve could solve for alpha due to the expint function, 
%so I used desmos instead. Nonetheless I keep the code commented below to 
%show the equation that was solved to find alpha.
%alpha = fsolve(@(alpha)alpha^2*exp(alpha^2)*expint(alpha^2)-(p0-p_plus)/(pL-p_plus),[0 5]);
alpha = 0.326338;

% Analytic Solution NOTE: This took over an hour to run on my machine
A = (p0-pL)/expint(alpha^2);
p_an = zeros(num_steps,num_latt_pts);
for i = 1:num_steps
    a_an = 2*alpha*sqrt(D*i*DT); % a = 2*a*sqrt(Dt)
    a_an_latt = round( a_an/DR + 1 ); % find cup size in lattice points
    for j = 1:a_an_latt
        p_an(i,j) = pL; % density at r<a is equal to pL
    end
    for j = a_an_latt:num_latt_pts % density at r>a is equal to p0 - AE1(r^2/4Dt)
        p_an(i,j) = p0-A*expint((j*DR)^2/(4*D*i*DT)); 
    end
end    





%Plot receptor density and cup size (change p to p_an for plots of
%analytical solution)

time_int = [0.5 5 T_MAX];
time_int_steps = round(time_int/DT);

clf;
%subplot(1,2,1);
x = 0:DR:(num_latt_pts-1)*DR;
plot(x,p(time_int_steps(1),:), ...
    x,p(time_int_steps(2),:), ...
    x,p(time_int_steps(3),:),'LineWidth',4);
set(gca,'FontSize',24);
legend(['t=' num2str(time_int(1)) 's'], ...
        ['t=' num2str(time_int(2)) 's'], ...
        ['t=' num2str(time_int(3)) 's']);
xlim([0 16]);
ylim([0 510]);
%xticks([0 2 4 6 8 10 12 14 16])
xlabel('Distance from centre of cup (µm)');
ylabel('Receptor density (µm^{-1})');
title('Receptor density profile during engulfment');

subplot(1,2,2);
t = 0:DT:(num_steps-1)*DT;
plot(t,a,t,2*alpha*sqrt(D*t),'LineWidth',3);
set(gca,'FontSize',24);
xlabel('Time (s)');
ylabel('Cup size (µm)');
legend('Numerical solution','Analytic solution','Location','NorthWest');
title('Phagocytic Cup Growth');