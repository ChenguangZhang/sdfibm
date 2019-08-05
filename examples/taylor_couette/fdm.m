ri = 0.3; % inner radius, ri must < 1
ro = 1.0; % outer radius
nr = 31; % radial resolution

U0 = 2*pi; % outer wall velocity 
nu = 0.1;  % kinematic viscosity

rhof = 1.0; % fluid density
rhos = 2.0; % core solid density

r = linspace(ri, 1, nr); % space axis
t = linspace(0, 6, 1201);%  time axix
dt = t(2) - t(1);
dr = r(2) - r(1);

% initial condition of flow field
u = 0*r;

beta = 0.0; % core wall speed
beta_ts = t*0; % time series of core wall speed
beta_ts(1) = 0.0; % initial speed of core wall

% time loop
for n = 2:length(t)
    dudr = [(-3*u(1)+4*u(2)-u(3)) (u(3:nr) - u(1:nr-2)) (3*u(nr)-4*u(nr-1)+u(nr-2))]/(2*dr);
    
    rdudr = r.*dudr;

    rdudrdr = [(-3*rdudr(1)+4*rdudr(2)-rdudr(3)) (rdudr(3:nr) - rdudr(1:nr-2)) (3*rdudr(nr)-4*rdudr(nr-1)+rdudr(nr-2))]/(2*dr);

    rhs = nu * (rdudrdr./r);
    
    % semi-implicit time update
    un = (u./dt + rhs)./(1/dt + nu./(r.^2));
    
    % now it's time to update beta, simplest Euler forward step
    % mind the evaluation of the wall shear stress, the -u(1)/ri was missed
    % previously
    dudr0 = ( (-3*u(1)+4*u(2)-u(3))/(2*dr)) - u(1)/ri;
    beta = beta + 4.0*rhof/rhos*nu/ri*dudr0 * dt;
    beta_ts(n) = beta;
    
    % set boundary condition
    un(1) = beta;
    un(nr)= U0;
    
    % update time
    u = un;
end

%% validation of steady-state (outer speed = 1, inner speed = 0)
% beta = 0 ;
% A = ri*(ri-beta)/(ri*ri-1);
% B = (ri*beta-1)/(ri*ri-1);
% plot(r, u,'ro', r, A./r + B*r, 'k-')

%% validation of  steady-state (outer speed = 1, du/dr = 0 at core)
% A = U0*ri*ri/(1+ri*ri);
% B = U0/(1+ri*ri);
% plot(r, u, 'ro', r, A./r + B*r, 'k-')
%% plot against openfoam
foam = load('./cloud.out');
exact = U0; % t->infty state
%
clf
hold on
plot(t(1:40:end), beta_ts(1:40:end)/ri/(2*pi), 'ko')
plot(foam(:,1), foam(:,16)/(2*pi), 'k-')
plot([0, t(end)], [exact, exact]/(2*pi), 'k--')
hold off
xlim([0 6])
ylim([0 1.1])
xlabel('t')
ylabel('\omega_{i}')
legend('finite difference', 'immersed boundary', '\omega_o = 1')