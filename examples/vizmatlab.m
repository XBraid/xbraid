function  vizmatlab(filename)

figure(1);
fname = strcat(filename, '.out.u.matlab');
U = importdata(fname);
surf(U);
shading interp;
title('Unknown (u)');
xlabel('Space');
ylabel('Time');

figure(2);
fname = strcat(filename, '.out.v.matlab');
V = importdata(fname);
surf(V);
shading interp;
title('Control (v)');
xlabel('Space');
ylabel('Time');

figure(3);
fname = strcat(filename, '.out.w.matlab');
W = importdata(fname);
surf(W);
title('Adjoint (w)');
shading interp;
xlabel('Space');
ylabel('Time');

commandwindow;

