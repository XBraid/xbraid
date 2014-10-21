 % Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 % Produced at the Lawrence Livermore National Laboratory. Written by 
 % Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 % Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 % 
 % This file is part of XBraid. Email xbraid-support@llnl.gov for support.
 % 
 % This program is free software; you can redistribute it and/or modify it under
 % the terms of the GNU General Public License (as published by the Free Software
 % Foundation) version 2.1 dated February 1999.
 % 
 % This program is distributed in the hope that it will be useful, but WITHOUT ANY
 % WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 % PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 % License for more details.
 % 
 % You should have received a copy of the GNU Lesser General Public License along
 % with this program; if not, write to the Free Software Foundation, Inc., 59
 % Temple Place, Suite 330, Boston, MA 02111-1307 USA



% 
% Syntax:
% spatialcoarsen(cfl, alpha, order, intorder, zeta)
%
% Input:
%   cfl: dt/h on fine grid
%   alpha: artificial dissipation coefficient
%   order: spatial order of accuracy (2 and 4 are implemented) (default 2)
%   intorder: interpolation order of accuracy (2 and 4 are implemented) (default 2)
%   zeta: residual filter coefficient. Apply Laplace smoother to coarse grid residual
% Output:

function spatialcoarsen(cfl, alpha, order, intorder, zeta)

iu = 1.0i;

% scaling the AD coefficient on the coarse grid didn't work out
scalead = 0;

if (nargin < 5)
  zeta=0;
end

if (nargin < 4)
  intorder=4;
end

if (nargin < 3)
  order=2;
end

ximax=1.0;
eps = 1e-4;
nk = 201;

xi = ximax*(-0.25*pi:(pi)/(nk-1):0.75*pi); % shifted the wave numbers to be positive
%xi = eps:ximax/(nk-1):ximax;

% set up the vectors for the spectral radius od D1, D2, D1+D2
kappa1 = zeros(nk,1);
kappa2 = zeros(nk,1);
kappasum = zeros(nk,1);

% fs(q) is the array of the filter symbol
fs = (1 - 2*zeta*sin(2*xi).^2);

if (intorder == 2 || intorder == 4)

% second order case
if (order == 2)
la2 = cfl*(iu*sin(xi) - alpha*(2*sin(0.5*xi)).^2);

phi = 1 + la2 + 0.5*la2.^2 + 1/6*la2.^3 + 1/24*la2.^4;

phi2 = phi.^2;

% add pi to xi
lapi2 = cfl*(iu*sin(xi+pi) - alpha*(2*sin(0.5*(xi+pi))).^2);

phipi = 1 + lapi2 + 0.5*lapi2.^2 + 1/6*lapi2.^3 + 1/24*lapi2.^4;

phipi2 = phipi.^2;

% scale AD coefficient
if (scalead==1)
  alphac=alpha/2;
else
  alphac=alpha;
endif 

% coarse grid in time and space (keeping cfl constant), doubling xi
lac = cfl*(iu*sin(2*xi) - alphac*(2*sin(2*0.5*xi)).^2);

phic = 1+ lac + 0.5*(lac).^2 + 1/6*(lac).^3 + 1/24*(lac).^4;

sxi = abs(phi2.*(phi2-phic)) ./ (1-abs(phic));
sxipi = abs(phipi2.*(phipi2-phic)) ./ (1-abs(phic));

if (intorder == 2)
% linear interpolation
  betap = (1-sin(0.5*xi).^2);
else
% cubic interpolation
  betap = (1-sin(0.5*xi).^4.*(1+2*cos(0.5*xi).^2));
endif

t1 = (1-betap).*( abs(phi2) + sxipi );
t2 = (betap).*( abs(phipi2) + sxi );

% temp
%t1 = (sin(0.5*xi).^2).*( sxipi );
%t2 = (cos(0.5*xi).^2).*( sxi );

for q=1:nk
  D1 = [ (1-betap(q)*fs(q))*abs(phi2(q)),  betap(q)*fs(q)*abs(phipi2(q)) ; ...
         (1-betap(q))*fs(q)*abs(phi2(q)),  (1-fs(q)+fs(q)*betap(q))*abs(phipi2(q)) ];

  D2 = fs(q)*[ betap(q)*sxi(q),     betap(q)*sxipi(q) ; ...
               (1-betap(q))*sxi(q), (1-betap(q))*sxipi(q) ];

  kappa = eig(D1);
  kappa1(q) = max(kappa);
  
  if (finite(D2))
    lambda = eig(D2);
    kappa2(q) = max(lambda);
    lambda = eig(D1+D2);
    kappasum(q) = max(lambda);
  else
    kappa2(q) = -0.1;
    kappasum(q)=-0.1;
  endif

endfor

elseif (order == 4)

% 4th order case

la4 = cfl*iu*sin(xi).*(1 + 2/3*sin(0.5*xi).^2) - cfl*alpha*(2*sin(0.5*xi).^4);

phi = 1 + la4 + 0.5*la4.^2 + 1/6*la4.^3 + 1/24*la4.^4;

phi2 = phi.^2;

% add pi to xi
lapi4 = cfl*iu*sin(xi+pi).*(1 + 2/3*sin(0.5*(xi+pi)).^2) - cfl*alpha*(2*sin(0.5*(xi+pi)).^4);

phipi = 1 + lapi4 + 0.5*lapi4.^2 + 1/6*lapi4.^3 + 1/24*lapi4.^4;

phipi2 = phipi.^2;

% coarse grid in time and space (keeping cfl constant), doubling xi

% 2nd order on coarse grid
%lac = cfl*(iu*sin(2*xi) - alphac*(2*sin(2*0.5*xi)).^2);

% scale AD coefficient
if (scalead==1)
  alphac=alpha/8;
else
  alphac=alpha;
endif 

% this is the 4th order discretization
lac = cfl*iu*sin(2*xi).*(1 + 2/3*sin(2*0.5*xi).^2) - cfl*alphac*(2*sin(2*0.5*xi).^4);

phic = 1+ lac + 0.5*(lac).^2 + 1/6*(lac).^3 + 1/24*(lac).^4;

sxi = abs(phi2.*(phi2-phic)) ./ (1-abs(phic));
sxipi = abs(phipi2.*(phipi2-phic)) ./ (1-abs(phic));

if (intorder == 2)
% linear interpolation
  betap = (1-sin(0.5*xi).^2);
else
% cubic interpolation
  betap = (1-sin(0.5*xi).^4.*(1+2*cos(0.5*xi).^2));
endif

t1 = (1-betap).*( abs(phi2) + sxipi );
t2 = (betap).*( abs(phipi2) + sxi );

for q=1:nk
%  D1 = [ (1-betap(q))*abs(phi2(q)),  betap(q)*abs(phipi2(q)) ; ...
%         (1-betap(q))*abs(phi2(q)),  betap(q)*abs(phipi2(q)) ];

%  D2 = [ betap(q)*sxi(q),     betap(q)*sxipi(q) ; ...
%         (1-betap(q))*sxi(q), (1-betap(q))*sxipi(q) ];
% with residual filtering
  D1 = [ (1-betap(q)*fs(q))*abs(phi2(q)),  betap(q)*fs(q)*abs(phipi2(q)) ; ...
         (1-betap(q))*fs(q)*abs(phi2(q)),  (1-fs(q)+fs(q)*betap(q))*abs(phipi2(q)) ];

  D2 = fs(q)*[ betap(q)*sxi(q),     betap(q)*sxipi(q) ; ...
               (1-betap(q))*sxi(q), (1-betap(q))*sxipi(q) ];

  kappa = eig(D1);
  kappa1(q) = max(kappa);
  
  if (finite(D2))
    lambda = eig(D2);
    kappa2(q) = max(lambda);
    lambda = eig(D1+D2);
    kappasum(q) = max(lambda);
  else
    kappa2(q) = -0.1;
    kappasum(q)=-0.1;
  endif

endfor

endif

else

% Fourier interpolation / restriction

% 4th order case

la4 = cfl*iu*sin(xi).*(1 + 2/3*sin(0.5*xi).^2) - cfl*alpha*(2*sin(0.5*xi).^4);

phi = 1 + la4 + 0.5*la4.^2 + 1/6*la4.^3 + 1/24*la4.^4;

phi2 = phi.^2;

% add pi to xi
lapi4 = cfl*iu*sin(xi+pi).*(1 + 2/3*sin(0.5*(xi+pi)).^2) - cfl*alpha*(2*sin(0.5*(xi+pi)).^4);

phipi = 1 + lapi4 + 0.5*lapi4.^2 + 1/6*lapi4.^3 + 1/24*lapi4.^4;

phipi2 = phipi.^2;

% coarse grid in time and space (keeping cfl constant), doubling xi

% scale AD coefficient
if (scalead==1)
  alphac=alpha/8;
else
  alphac=alpha;
endif 

% this is the 4th order discretization with 3rd order AD
lac = cfl*iu*sin(2*xi).*(1 + 2/3*sin(2*0.5*xi).^2) - cfl*alphac*(2*sin(2*0.5*xi).^4);

% 2nd order on coarse grid
%lac = cfl*(iu*sin(2*xi) - alphac*(2*sin(2*0.5*xi)).^2);

% this is the 4th order discretization with 1st order AD
%lac = cfl*iu*sin(2*xi).*(1 + 2/3*sin(2*0.5*xi).^2) - cfl*alphac*(2*sin(2*0.5*xi).^2);

phic = 1 + lac + 0.5*(lac).^2 + 1/6*(lac).^3 + 1/24*(lac).^4;

sxi = abs(phi2.*(phi2-phic)) ./ (1-abs(phic));
sxipi = abs(phipi2.*(phipi2-phic)) ./ (1-abs(phic));

% explicitly using betap=1
t1 = abs(phipi2);
t2 = sxi;

for q=1:nk
  D1 = [ (1-fs(q))*abs(phi2(q)), 0              ; ...
         0,                      abs(phipi2(q)) ];

  D2 = fs(q)*[ sxi(q), 0 ; ...
               0,      0 ];

  kappa = eig(D1);
  kappa1(q) = max(kappa);
  
  if (finite(D2))
    lambda = eig(D2);
    kappa2(q) = max(lambda);
    lambda = eig(D1+D2);
    kappasum(q) = max(lambda);
  else
    kappa2(q) = -0.1;
    kappasum(q)=-0.1;
  endif

endfor

printf("Fourier interpolation/restriction\n");
%return;

endif

% plot the result
printf("Spatial order=%i, cfl=%e, alpha=%e, alphac=%e, zeta=%e, max spec-radius(D1+D2)=%e\n", order, ...
        cfl, alpha, alphac, zeta, max(kappasum));
h=plot(xi/pi, abs(phic), "g", xi/pi, fs, "r.", xi/pi, abs((phi.^2-phic).*phi.^2), "b", ...
     xi/pi, abs(phipi.^2), "r", xi/pi, kappa1, "c", xi/pi, kappasum, "k.");
set(h,"linewidth",2);
axis([-0.25 0.75 -.1 2]);
legend("abs(phic)","filter symbol", "abs((phi^2-phic)phi^2)", ...
       "abs(phipi^2)", "spec(D1)", "spec(D1+D2)"); 
set(gca,"fontsize",18);
%pause;
return;

% 6th order case (not updated)
la2 = cfl*iu*sin(xi).*(1 + 2/3*sin(0.5*xi).^2 + 8/15*sin(0.5*xi).^4) - cfl*alpha*(2*sin(0.5*xi).^6);

phi = 1 + la2 + 0.5*la2.^2 + 1/6*la2.^3 + 1/24*la2.^4;

phi2 = phi.^2;

% coarse grid in time P(2*la2)
phic = 1+ 2*la2 + 0.5*(2*la2).^2 + 1/6*(2*la2).^3 + 1/24*(2*la2).^4;

cf = abs((phi2-phic).*phi2)./(1-abs(phic));

% plot the 2nd order case:
plot(xi/pi, abs(phi2), "k", xi/pi, 1-abs(phic), "g", xi/pi, abs(phi2-phic), "b", xi/pi, cf, "r");
legend("abs(phi^2)","1-abs(phic)","abs(phi2 - phic)","conv fact"); 
return;

ka = eps:ximax/(nk-1):ximax;
pc = 1 + 2*iu.*ka - 2*ka.^2 - 8/6*iu*ka.^3 + 16/24*ka.^4;
pcb = conj(pc);

pc2 = real(pc.*pcb);
%loglog(ka,1- pc2, 'r');
%legend('1 - pc*pcb');

%pause;

maxphi = max([abs(phi),abs(phic)]);
minphi = min([abs(phi),abs(phic)]);

%plot(xi/pi, real(phi)-1, 'k', xi/pi, imag(phi) ,'b', xi/pi, real(phic)-1, 'r', xi/pi, imag(phic), 'm'); 
%xlabel('xi/pi')
%legend('real(phi)-1', 'imag(phi)', 'real(phic)-1', 'imag(phic)');

plot(xi/pi, abs(phi.^2), 'b')
xlabel('xi/pi')
legend('abs(phi^2)');

pause;

plot(xi/pi, real((phi.^2 - phic).*phi.^2), 'b', xi/pi, imag((phi.^2 - phic).*phi.^2), 'r', 
     xi/pi, abs((phi.^2 - phic).*phi.^2), 'k')
xlabel('xi/pi')
legend('real((phi^2 - phic)*phi^2)', 'imag((phi^2 - phic)*phi^2)', 'abs((phi^2 - phic)*phi^2)');

printf("Max |(phi^2 - phic)*phi^2)| = %e\n", max(abs((phi.^2 - phic).*phi.^2)));
