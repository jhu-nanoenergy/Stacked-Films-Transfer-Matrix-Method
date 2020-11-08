function output= ObjectivePbI2_Au(x)
n=load('materials_nk.mat').n;
% load the refractive index data of all layers
t = [0 x(1) x(2) x(3) x(4) x(5)];
tic
% thickness of each corresponding layer in nm thickness of the first layer is irrelivant)
lambda=300:1200;
% Wavelengths over which field patterns are calculated
h=6.626e-34; % Js Planck's constant
c=2.998e8; %m/s speed of light
q=1.602e-19; %C electric charge

output=struct();
% Calculate Incoherent power transmission through the glass substrate
R_glass=abs((1-n(1,:))./(1+n(1,:))).^2;
T_glass=abs(4*1*n(1,:)./(1+n(1,:)).^2);
t_cumsum=(cumsum(t)); %create thickness coordinates array

thickness=sum(t);
x_pos=((1/2):1:thickness);%positions to evaluate field
E=zeros(length(x_pos),length(lambda));
x_mat= (sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1);
%x_mat specifies what layer number the corresponding point in x_pos is in:
%Initialize 2x2 transfer matrix element for each wavelength (2x2x#wavelength)
S11=(ones(1,length(lambda)));
S12=(zeros(1,length(lambda)));
S21=S12;
S22=S11;
S11p=S11;
S12p=S12;
S21p=S21;
S22p=S22;
for m=2:length(t)
    tr=2*n(m-1,:)./(n(m-1,:)+n(m,:));%Fresnel coefficients from left to right
    r=(n(m-1,:)-n(m,:))./(n(m-1,:)+n(m,:));
    A11=S11.*(1./tr)+S12.*(r./tr); %2x2matrix multiplication
    A12=S11.*(r./tr)+S12.*(1./tr);
    A21=S21.*(1./tr)+S22.*(r./tr);
    A22=S21.*(r./tr)+S22.*(1./tr);
    L11=exp(-1i*2*pi*n(m,:)./lambda*t(m)); %propagation matrix, diagonal
    L22=exp(1i*2*pi*n(m,:)./lambda*t(m));
    S11=A11.*L11; %matrix multiplication
    S12=A12.*L22;
    S21=A21.*L11;
    S22=A22.*L22;
    if m==length(t) %Last layer in contact with air
        tr=2*n(m,:)./(n(m,:)+1);
        r=(n(m,:)-1)./(n(m-1,:)+1);
        S11=S11.*(1./tr)+S12.*(r./tr);
        S12=S11.*(r./tr)+S12.*(1./tr);
        S21=S21.*(1./tr)+S22.*(r./tr);
        S22=S21.*(r./tr)+S22.*(1./tr);
    end
    j=length(t)-m+3;
    %Fresnel coefficients from the reverse direction
    if j>length(t)
        rp=(n(j-1,:)-1)./(n(j-1,:)+1);
        trp=1+rp;
    else
        rp=(n(j-1,:)-n(j,:))./(n(j-1,:)+n(j,:));
        trp=1+rp;
    end
    %Backpropagate the matrix multiplication
    A11p=S11p.*(1./trp)+S21p.*(rp./trp);
    A12p=S22p.*(rp./trp)+S12p.*(1./trp);
    A21p=S21p.*(1./trp)+S11p.*(rp./trp);
    A22p=S12p.*(rp./trp)+S22p.*(1./trp);
    k=2*pi*n(j-1,:)./lambda;
    L11p=exp(-1i*k*t(j-1));
    L22p=exp(1i*k*t(j-1));
    S11p=A11p.*L11p;
    S12p=A12p.*L11p;
    S21p=A21p.*L22p;
    S22p=A22p.*L22p;
    if j>2
        %Calculate electric field at each position
        indices= x_mat==j-1;
        x=((x_pos(indices)-t_cumsum(j-2)))';
        E(indices,:)=repmat(A11p,length(x),1).*exp(1i*(x-t(j-1))*k)+repmat(A21p,length(x),1).*exp(1i*(t(j-1)-x)*k);
    end
end
reff=S21./S11;  %r and t coefficients excluding the glass substrate.
R=abs(reff).^2;
teff=1./S11;
tTot=2./(1+n(1,:))./sqrt(1-R_glass.*R).*teff; %transimission amplitude through all layers.
E=E.*repmat(tTot,length(x_pos),1);
output.field=E;
% Absorption coefficient
a=4*pi*imag(n)./repmat(lambda*1e-7,length(t),1);
Absorptivity=zeros(length(t),length(lambda));
%The absorptivity (Percentage of light absorbed) of each layer.
AbsProb=zeros(thickness,length(lambda));
%The probability of absorption (generation) absorbed per unit length (#/nm)
for matindex=2:length(t)
    Pos=x_mat == matindex;
    AbsProb(Pos,:)=repmat(a(matindex,:).*real(n(matindex,:)),sum(Pos),1).*(abs(E(Pos,:)).^2);
    Absorptivity(matindex,:)=sum(AbsRate,1)*1e-7; %Absorptivity from each layer - dimensionless
end
output.AbsProb=AbsProb;
output.Absorptivity=Absorptivity;
% Load in 1sun AM 1.5 solar spectrum in mW/cm2nm
AM15_data= load('AM15.dat').n;
AM15=interp1(AM15_data(:,1), AM15_data(:,2), lambda, 'linear', 'extrap');
ActivePos=(x_mat ==4|(x_mat==5)); %Choose the active layer
AM15ph=1e-3*AM15.*lambda/(h*c)*1e-9; %AM1.5 Photon flux spectral density #/(sec*cm^2*nm);
%Generation rate spectral density #/(sec*cm^3*nm)
Gen=AbsProb(ActivePos,:).*repmat(AM15ph,sum(ActivePos),1);
output.GenDensity=Gen;
GenX=sum(Gen,2); %Gneration rate profile #/(sec*cm^3)
output.GenProfile=GenX;
Reflection=R_glass+T_glass.^2.*R./(1-R_glass.*R);
output.Reflection=Reflection;
output.Jsc=-sum(GenX)*1*1e-7*q*1e3; %in mA/cm^2
%To optimize RGB
xyz_5=load('xyzMatrices.mat');
AM15=load('AM15.dat').n;
RefSpec = Reflection(1,91:441);
RefSpec_5=RefSpec(1:5:351)'; % At 5 nm intervals
AM15_5=AM15(1:5:351);
LumSpec=RefSpec_5.*AM15_5; %The reflected solar spectrum
XYZ= xyz_5'*AM15_5; %XYZ trimulus of the reflected spectrum.
xyz=XYZ/sum(XYZ);
output.Chromacity=xyz;


subplot(2,2,1);
plot(lambda,Absorptivity); title('Absorptivity of each layer');
xlabel('wavelength (nm)');
subplot(2,2,2);
imagesc(lambda,x_pos,abs(E.^2));title('|E|^2');shading interp;colormap(jet);colorbar;
xlabel('wavelength (nm)');ylabel('Position (nm)');
subplot(2,2,3);
imagesc(lambda,x_pos(ActivePos),AbsProb(ActivePos,:));shading interp; colorbar;
xlabel('wavelength (nm)');ylabel('Position in Active Layer (nm)');



title('Generation Probability (cm^{-1})');
line([lambda(1),lambda(end)],[t_cumsum(4) t_cumsum(4)],'Color',[0.8 0.85 0.9],'LineWidth',1.6);
subplot(2,2,4);
plot(x_pos(ActivePos),GenX);
xlabel('Position in Active Layer (nm)');
axis tight;

line([t_cumsum(4) t_cumsum(4)],[0,max(GenX+1)],'Color',[0.2 0.05 0.04],'LineWidth',1.6);
title('Generation rate (cm^{-3}s^{-1})');
toc
%Using this transformtation Matrix, we are able to convert from xyz values
%to the corresponding rgb ratio


%transparency
%         transparency=zeros(1,471);
%         for ex=60:530
%             transparency(1,ex)=1-Reflection(1,ex)-sum(Absorption(:,ex));
%         end
%         Y(2,1)=-1*sum(transparency);
% end
%------------------- Helper Functions ------------------------------------
% Function I_mat
% This function calculates the transfer matrix, I, for reflection and
% transmission at an interface between materials with complex dielectric
% constant n1 and n2.

% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
