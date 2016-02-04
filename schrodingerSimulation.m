%% Standard wave equation simulation.

%--------------------------------------------------------------------------
% Set me!!!
N = 101;        % Spatial size
numIter = 200;  % Temporal size
pulsePos = 11;  % Initial pulse position
pulseWidth = 5; % Initial pulse width
%--------------------------------------------------------------------------

% Planck units
% c^2 := 1
% dx  := 1
% dt  := 1
x = (1:N)';

% Construct gaussian pulse.
f_pulse = exp(-(x-pulsePos).^2/pulseWidth);
f_pulse = f_pulse/norm(f_pulse);

% Construct 2nd order derivative operator (L = - d^2/dx^2)
e = ones([N,1]);
D = spdiags([-e/2,e/2], [-1,1], N, N);
L = spdiags([-e, 2*e, -e], [-1,0,1], N, N);

% Initial conditions.
f_prev = f_pulse;
f = f_pulse;
fim_prev = zeros([N,1]);
fim = zeros([N,1]);

fprintf(1,'||f||=%d\n', norm(f));
norms = zeros([numIter,1]);
norms(1) = norm(f);

h = figure;
for i=2:numIter
    % Jacobi iteration on d^2/dt^2 f = -Lf.
    f_next = -L*f + 2*f - f_prev;
    f_prev = f;
    f = f_next;
    
    % Just checking...
    norms(i) = norm(f+1i*fim);
    
    figure(h);
    % Why do I have to reset the damn axis?
    plot(x, f, 'b-'); hold on; plot(x, fim, 'r-'); hold off; axis([1,N,-2,2]);
    drawnow;
    pause(0.1);
end

fprintf(1,'||f||=%d\n', norm(f));
figure; plot(norms);
title('Iteration Stability');
ylabel('Instantaneous Wave Norm');
xlabel('Iteration');