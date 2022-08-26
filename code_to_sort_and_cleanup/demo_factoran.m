a = sin((1:1000)*0.1) + rand(1000,1)'/10;
b = cos((1:1000)*0.2) + rand(1000,1)'/2;

N = 100;

a = [];
b = [];
noise_amp = 10;
scaling_amp = 10;
noise_a = linspace(noise_amp - noise_amp/2,noise_amp + noise_amp/2,1000);
noise_b = linspace(noise_amp - noise_amp/2,noise_amp + noise_amp/2,1000);
scaling_a = linspace(0,scaling_amp,N);
scaling_b = linspace(0,scaling_amp,N);
phase_a   = linspace(0,pi/4,N);
phase_b   = linspace(3*pi/4,5*pi/4,N);

phase_a = phase_a(randperm(numel(phase_a)));
phase_b = phase_b(randperm(numel(phase_b)));
scaling_a = scaling_a(randperm(numel(scaling_a)));
scaling_b = scaling_b(randperm(numel(scaling_b)));

for idx = 1:100
    t = linspace(0,4*pi,1000);
    a = [a; sin(t + phase_a(idx)) * scaling_a(idx) + noise_a(randperm(numel(noise_a)))];
    b = [b; sin(t + phase_b(idx)) * scaling_b(idx) + noise_b(randperm(numel(noise_b)))];
end

c = a+b;
figure();plot(c')
[LoadingsPM, specVarPM, T, stats, F] = factoran(double(c'), 3);
figure();scatter(LoadingsPM(:,1), scaling_a, 'b', 'filled'); hold on;scatter(LoadingsPM(:,2), scaling_b, 'r', 'filled');
figure();plot(F)