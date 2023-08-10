% NSCM005 CA2 - the Chay-Keizer model of plateau bursting in endocrine
% cells

% Define parameters
DT = 0.9; % time step, [ms]
T_MAX = 40000; % max time, [ms]

g_Ca = 1000; % maximum conductance of Ca channel [pS]
g_KCa = 400; % maximum conductance of K(Ca) channel [pS]
g_K = 2700; % maximum conductance of K channel [pS]
g_KATP = 180; % maximum conductance of K(ATP) channel [pS]
V_Ca = 25; % [mV]
V_K = -75; % [mV]
C_m = 5300; % [fF]
tau_n = 18.7; % [ms]
alpha = 9e-6; % [fA^-1 uM ms^-1]
f = 0.00025; % [-]    % 3 different values
k_PMCA = 0.5; % [ms^-1]
K_d = 0.3; % [uM]
v_n = -12; % [mV]
v_m = -20; % [mV]
s_n = 5; % [mV]
s_m = 12; % [mV]

% Declare/allocate variables
num_steps = round(T_MAX/DT); % number of time steps [-]

V = zeros(num_steps, 1);
n = zeros(num_steps, 1);
c = zeros(num_steps, 1);
I_K = zeros(num_steps, 1);
I_KATP = zeros(num_steps, 1);
I_Ca = zeros(num_steps, 1);
I_KCa = zeros(num_steps, 1);


% Initial conditions
V(1) = -65;
n(1) = 0.0;
c(1) = 0.1;
I_K(1) = 0;
I_KATP(1) = 0;
I_Ca(1) = 0;
I_KCa(1) = 0;
v_n = -12;

% Main Loop
for i = 2:num_steps
    
    % calculate activation functions
    m_inf = ( 1+exp((v_m-V(i-1))/s_m) )^-1;
    n_inf = ( 1+exp((v_n-V(i-1))/s_n) )^-1;
    s_inf = c(i-1)^3/(c(i-1)^3+K_d^3);
    
    % calculate currents
    I_K(i) = g_K*n(i-1)*(V(i-1)-V_K);
    I_KATP(i) = g_KATP*(V(i-1)-V_K);
    I_Ca(i) = g_Ca*m_inf*(V(i-1)-V_Ca);
    I_KCa(i) = g_KCa*s_inf*(V(i-1)-V_K);
    
    % calculate variable changes
    DV = -1/C_m*(I_Ca(i)+I_K(i)+I_KCa(i)+I_KATP(i))*DT;
    Dn = 1/tau_n*(n_inf-n(i-1))*DT;
    Dc = -f*(alpha*I_Ca(i)+k_PMCA*c(i-1))*DT;
    
    % update values
    V(i) = V(i-1)+DV;
    n(i) = n(i-1)+Dn;
    c(i) = c(i-1)+Dc;
end



% Plotting
t = 0:DT:(num_steps-1)*DT;

subplot(1,2,1);
plot(t,V);
set(gca,'FontSize',24);
xlabel('Time (ms)');
ylabel('V (mV)');
ylim([-80 0]);

subplot(1,2,2);
plot(t,c);
set(gca,'FontSize',24);
xlabel('Time (ms)');
ylabel('c (uM)');
ylim([0.06 0.2]);



