% last edited by Lei Tian (lei_tian@alum.mit.edu), 04/29/2015

% reference: 
% Lei Tian and Laura Waller, "Quantitative differential phase contrast
% imaging in an LED array microscope," Opt. Express 23, 11394-11403 (2015)


% % Define Fourier operators
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
% F = @(x) fftshift(fft2(x));
% Ft = @(x) ifft2(ifftshift(x));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wavelength of illumination, assume monochromatic
% R: 624.4nm +- 50nm
% G: 518.0nm +- 50nm
% B: 476.4nm +- 50nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 0.624;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical aperture of the objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA = 0.25;
Np=[1200 1920];
% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA
dx0 = 1/um_m/2;

binning = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnification of the system,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mag = 10;

dpix_c = 5.86*binning; %6.5um pixel size on the sensor plane
% effective image pixel size on the object plane
dpix_m = dpix_c/mag; 
% dpix_m = .498; 

% FoV in the object space
FoV = Np*dpix_m;
% sampling size at Fourier plane is always = 1/FoV
if mod(Np,2) == 1
    du = 1./dpix_m/(Np-1);
else
    du = 1./FoV;
end

% low-pass filter diameter set by the NA = bandwidth of a single measurment
% in index
% N_NA = round(2*um_m/du_m);
% generate cutoff window by NA
m = [1:Np(2)]-round((Np(2)+1)/2); 
n = [1:Np(1)]-round((Np(1)+1)/2);
[mm,nn] = meshgrid(m,n);

um_idx = um_m./du;
ridx = sqrt((mm/um_idx(2)).^2+(nn/um_idx(1)).^2);
% assume a circular pupil function, lpf due to finite NA
w_NA = double(ridx<=1);
% h = fspecial('gaussian',10,5);
% w_NA = imfilter(w_NA,h);

% support of OTF is 2x of ATF(NA)
Ps_otf = double(ridx<=2);

phC = ones(Np);
aberration = ones(Np);
%aberration = ones(N_m);
pupil = w_NA.*phC.*aberration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up image corrdinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncent = [1080,1280]/binning;
nstart = ncent-Np/2;
img_center = (nstart-ncent+Np/2)*dpix_m;
img_start = nstart*dpix_m;
img_end = (nstart+Np)*dpix_m;

%% spatial frequency coordinates
umax = 1/2/dpix_m;
dv = 1/dpix_m/Np(1); du = 1/dpix_m/Np(2);
u0 = [-umax:du:umax-du];
v0 = [-umax:dv:umax-dv];
[u,v] = meshgrid(u0,v0);

%% LED array geometries and derived quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacing between neighboring LEDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_led = 4.5e3; %4mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from the LED to the object
% experientally determined by placing a grating object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_led = 17e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up LED coordinates
% h: horizontal, v: vertical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lit_cenv = 4;
lit_cenh = 4;
vled = [0:7]-lit_cenv;
hled = [0:7]-lit_cenh;

[hhled,vvled] = meshgrid(hled,vled);
rrled = sqrt(hhled.^2+vvled.^2);

% corresponding angles for each LEDs
dd = sqrt((-hhled*ds_led-img_center(1)).^2+(-vvled*ds_led-img_center(2)).^2+z_led.^2);
sin_thetav = (-hhled*ds_led-img_center(1))./dd;
sin_thetah = (-vvled*ds_led-img_center(2))./dd;

tan_thetav = (-hhled*ds_led-img_center(1))./z_led;
tan_thetah = (-vvled*ds_led-img_center(2))./z_led;

illumination_na = sqrt(sin_thetav.^2+sin_thetah.^2);
% corresponding spatial freq for each LEDs
%
vled = sin_thetav/lambda;
uled = sin_thetah/lambda;
% spatial freq index for each plane wave relative to the center
idx_u = round(uled/du);
idx_v = round(vled/dv);




