% Program to demonstrate Talbott effect (or self-imaging) and fractional Talbott effect (or contrast reversal) computationally. 
% An amplitude grating with transmission function
% t(x, y) = 0.5[1 + cos(2 \pi f_0 x)] is considered and its diffraction pattern
% is computed at first and second Talbott distances. 
% The program also demostrates fractional Talbott effect (contrast reversal).

clc;
close all;
pixel = 2 * (10^(-6));          % pixel size in m
lambda = 0.5 * (10^(-6));       % wavelength of em field in m
k = 2*pi / lambda ;             % wavenumber or magnitude of wavevector
f_0 = 7000 ;                    % 'grating frequency'

[x, y] = meshgrid ( -256: 1 : 255);
[fx , fy] = meshgrid ( -0.5 : 1/512 : 0.5-(1/512) ) ;
x = x * pixel ; y = y * pixel ; fx = fx / pixel ; fy = fy / pixel ;

t = 0.5*( 1 + cos( 2*pi*f_0*x ) );  % grating transmission function

L = 1/f_0  ;                                       % period of the grating
n = 1;                                             % n is an integer
z_Talbott = (2 * n * L^2)/lambda ;                 % First Talbott distance
z_frac_Talbott = ((2*n+1) * L^2)/lambda ;          % Fractional Talbott distance
z_second_Talbott = (2 * 2 * L^2)/lambda ;          % Second Talbott distance

grating_FT = fftshift(fft2(ifftshift( t )));
wave_FT_1 = grating_FT;

% Transfer function for angular spectrum method for first Talbott distance
H_ang_Talbott = exp( 1i* z_Talbott *sqrt( k^2 -(  4* pi*pi*(fx.^2 + fy.^2) ) ) );
wave_FT_2_ang_Talbott = wave_FT_1 .* H_ang_Talbott;
wave_2_ang_Talbott = fftshift(ifft2(ifftshift( wave_FT_2_ang_Talbott )));

% Transfer function for angular spectrum method for second Talbott distance
H_ang_second_Talbott = exp( 1i* z_second_Talbott *sqrt( k^2 -(  4* pi*pi*(fx.^2 + fy.^2) ) ) );
wave_FT_2_ang_sec_Talbott = wave_FT_1 .* H_ang_second_Talbott;
wave_2_ang_sec_Talbott = fftshift(ifft2(ifftshift( wave_FT_2_ang_sec_Talbott )));

% Transfer function for angular spectrum method for distance correspondig
% to fractional Talbott effect
H_ang_frac_Talbott = exp( 1i* z_frac_Talbott *sqrt( k^2 -(  4* pi*pi*(fx.^2 + fy.^2) ) ) );
wave_FT_2_ang_frac_Talbott = wave_FT_1 .* H_ang_frac_Talbott;
wave_2_ang_frac_Talbott = fftshift(ifft2(ifftshift( wave_FT_2_ang_frac_Talbott )));

figure,  
subplot(2,2,1), imagesc(t); colormap(gray);
title('Grating transmission function');
subplot(2,2,2), imagesc(abs(wave_2_ang_Talbott)); colormap(gray); 
title('Self imaging: Intensity at first Talbott distance');
subplot(2,2,3), imagesc(abs(wave_2_ang_sec_Talbott)); colormap(gray);
title('Self imaging: Intensity at second Talbott distance');
 subplot(2,2,4), imagesc(abs(wave_2_ang_frac_Talbott)); colormap(gray);
title('Contrast reversal: Fractional Talbott effect');
