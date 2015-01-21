close all

gs = 0.5;
muB = 9.27401e-24;
hbar = 1.05457e-34;

dt = 1.e-6;
nsteps = ceil(0.015 / dt * 2);

time    = zeros(1,nsteps);
time(1) = -0.015;

Bx = 1.05e-7;
C = 0.0025;

Bz =-C*time(1);
magB = sqrt(Bx^2 + Bz^2);
nx = Bx / magB;
nz = Bz / magB;

psi      = zeros(2,nsteps);
psi(:,1) = 0.5 * [1 + nx + nz;
             1 + nx - nz] / sqrt(1+nx);
        
Pup = 0.5*[1+nz, nx; nx, 1-nz];
Pdn = 0.5*[1-nz, -nx; -nx, 1+nz];
   
pup    = zeros(1,nsteps);
pup(i) = psi(:,i)' * Pup * psi(:,i);
pdn    = zeros(1,nsteps);
pdn(i) = psi(:,i)' * Pdn * psi(:,i);

for i=2:nsteps
    
    time(i) = time(i-1) + dt;
    
    Bz =-C*time(i);
    magB = sqrt(Bx^2 + Bz^2);
    nx = Bx / magB;
    nz = Bz / magB;

    theta = gs*muB*magB*dt / 2. / hbar;

    U = [cos(theta) - 1i * nz * sin(theta), -1i * nx * sin(theta);
         -1i * nx * sin(theta), cos(theta) + 1i * nz * sin(theta)];
     
    psi(:,i) = U * psi(:,i-1);
    
    Pup = 0.5*[1+nz, nx; nx, 1-nz];
    Pdn = 0.5*[1-nz, -nx; -nx, 1+nz];
    
    pup(i) = real(psi(:,i)' * Pup * psi(:,i));
    pdn(i) = real(psi(:,i)' * Pdn * psi(:,i));
end

pops = conj(psi).*psi;

k = gs*muB*Bx^2 / (hbar*C); 

fprintf('Probability of spin flip is %f\n', exp(-0.5*k*pi) );
fprintf('This simulation predicts %f\n', pops(1,end) );

figure()
plot(time*1e3, pops(1,:), time*1e3, pops(2,:), time*1e3, pops(1,:)+pops(2,:))
ylabel('| \uparrow,\downarrow (t) \rangle')
xlabel('time (ms)')
legend('<\uparrow|\uparrow>','\langle\downarrow\rangle^2','\langle\psi\rangle^2')

figure()
plot(time*1e3, pup, time*1e3, pdn, time*1e3, pup+pdn)
ylabel('| \phi (t) \rangle')
xlabel('time (ms)')
