
using DelimitedFiles
using PlotlyJS
#using Plots
#plotlyjs()

function vizjulia(filename::AbstractString)
    
    # Plot u
    fname = string(filename, ".out.u.julia")
    u = readdlm(fname, ',')'                    # transpose after reading
    (nx,nt) = size(u)
    x = range(0, stop=1, length=nx)
    t = range(0, stop=1, length=nt)
    layout = Layout(
        title = "Unknown (u)", autosize = false,
        width = 500, height = 500
    )
    pu = plot(surface(z=u, x=x, y=t), layout)
    
    # Plot v
    fname = string(filename, ".out.v.julia")
    v = readdlm(fname, ',')'                    # transpose after reading
    (nx,nt) = size(v)
    x = range(0, stop=1, length=nx)
    t = range(0, stop=1, length=nt)
    layout = Layout(
        title = "Control (v)", autosize = false,
        width = 500, height = 500
    )
    pv = plot(surface(z=v, x=x, y=t), layout)
    
    # Plot w
    fname = string(filename, ".out.w.julia")
    w = readdlm(fname, ',')'                    # transpose after reading
    (nx,nt) = size(w)
    x = range(0, stop=1, length=nx)
    t = range(0, stop=1, length=nt)
    layout = Layout(
        title = "Adjoint (w)", autosize = false,
        width = 500, height = 500
    )
    pw = plot(surface(z=w, x=x, y=t), layout)

    p = [pu pv pw]
    relayout!(p, height=500, width=1300)
    p
    
end



# # Matlab
# 
# figure(1);
# fname = strcat(filename, '.out.u.matlab');
# U = importdata(fname);
# surf(U);
# shading interp;
# title('Unknown (u)');
# xlabel('Space');
# ylabel('Time');
# 
# figure(2);
# fname = strcat(filename, '.out.v.matlab');
# V = importdata(fname);
# surf(V);
# shading interp;
# title('Control (v)');
# xlabel('Space');
# ylabel('Time');
# 
# figure(3);
# fname = strcat(filename, '.out.w.matlab');
# W = importdata(fname);
# surf(W);
# title('Adjoint (w)');
# shading interp;
# xlabel('Space');
# ylabel('Time');
# 
# commandwindow;

