function kx = vwa_kxgrid(Nx, dx)
%VWA_KXGRID  kx grid consistent with MATLAB fft ordering (no fftshift).
% Updated to ensure function recognition.
%
% For fft(x), frequency indices correspond to:
%   0, 1, 2, ..., floor(N/2), -ceil(N/2)+1, ..., -1
%
% Returns kx in rad/m.

    if Nx < 2
        kx = 0;
        return;
    end

    dk = 2*pi/(Nx*dx);

    kpos = 0:floor(Nx/2);
    kneg = -ceil(Nx/2)+1:-1;

    k = [kpos, kneg];
    kx = (dk * k).';
end
