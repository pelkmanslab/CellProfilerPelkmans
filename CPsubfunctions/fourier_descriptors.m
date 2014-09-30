function fd = fourier_descriptors(b)
% Fourier decriptors for boundary vector

% 14.3.2008 (C) Pekka Ruusuvuori

N = length(b);
L = 7;
z = complex(b(:,2),-b(:,1));
c = fft(z);
fd = [abs(c(N+1-L:N)); abs(c(3:L+1))] / abs(c(2));
