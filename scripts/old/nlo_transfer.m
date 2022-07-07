%% Define parameters.
m = 1;
c1 = 0.5;
c2 = 0.5;
c = c1 + c2;
d = 1;
xt = 1;
omega = sqrt(c/m);
D = d/(2*sqrt(m*c));

%% Define forcing and IC.
Fmax = 0.9;
t2 = .9;
alpha = Fmax/t2;
x0 = 0.0;

%% Define system and its corresponding transfer function.
A = [0 1; -c/m -d/m];
B = [0 1/m]';
C = [1 0];
sys = ss(A,B,C,0);
[b, a] = ss2tf(A,B,C,0);
G = tf(b, a);

%% Choose which forcing to apply.
t = 0:0.01:100;
n = length(t);
s = tf("s");
choose_forcing = 1;
if choose_forcing == 1
    u = Fmax*ones(n,1);
    U = Fmax/s;
elseif choose_forcing == 2
    u = min(alpha*t, Fmax);
    U = alpha/s^2 * (1 - exp(-t2*s));
elseif choose_forcing == 3
    amp = 0.5;
    omega_sin = omega * sqrt(1-2*D^2);
    u = amp * sin( omega_sin*t );
    U = amp * omega_sin / (s^2 + omega_sin^2);
elseif choose_forcing == 4
    u = Fmax .* ( 1 - exp(-alpha*t) );
    U = alpha / (s * (s + alpha));
else
    u = Fstep*ones(n,1) + min(alpha*t, Fmax);
    U = alpha/s^2 * (1 - exp(-t2*s)) + Fstep/s;
end

%% Compute the resulting transformed function.
% Y0 = (x0 * s + 2*d*omega*x0) * tf(1, a)
Y = G*U % + Y0
Ystat = xt / s

%% Bode plot.
fig2 = figure();
bode(G, U, Y, Ystat)
grid on

%% Lsim-check to verify in time domain.
[y,t,x] = lsim(sys, u, t, [x0, 0]);
fig3 = figure();
subplot(2,1,1);
plot(t, y, t, xt*ones(n,1));
grid on
subplot(2,1,2);
plot(t, u);
grid on

%% Compute Laplace tipping criterion
fig4 = figure();
w = logspace(-3, 3, 1000);
[magY, phaseY, wout] = bode(Y, w);
[magYstat, phaseYstat, wout] = bode(Ystat, w);
ews = (squeeze(magYstat) - squeeze(magY)) ./ squeeze(magYstat);
semilogx(wout, ews);
grid on

tipping_predict = abs(sum(ews(ews < 0)))
tipping_ground_truth = sum( y > xt )
