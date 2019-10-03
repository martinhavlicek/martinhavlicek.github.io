##
using Plots

x = zeros(1,6);
f = copy(x);
u = zeros(1,400);

u[1,10:20] .= 1;


EE = -0.5;
IE = 0.5;
II = 0.1;
EI = copy(II);
C  = 1/16;
de1 = 0.6;
ga  = 1.5;
de2 = 0.6;
nr  = 3;
al  = 0.35;
tt = 2;
ve = 10;
dt = 0.1;

X = zeros(length(u),size(x,2));
for t = 1:length(u)

x[:,3:5]= exp.(x[:,3:5]);

f[:,1]  = f[:,1] .+ dt.*(EE.*x[:,1] .- IE.*x[:,6] .+ C*u[t]);
f[:,6]  = f[:,6] .+ dt.*(-II.*x[:,6] + EI.*x[:,1]);


# # m = (f-1 + nr)/nr - oxygen metabolism
# #--------------------------------------------------------------------------
m        = (x[:,3] .- 1 .+ nr)./nr;
# implement differential state equation y = dx/dt [hemodynamic]
#--------------------------------------------------------------------------
f[:,2]   = f[:,2] .+ dt.*(x[:,1] .- de1.*(x[:,2]));
dfin     = (ga.*x[:,2] .- de2.*(x[:,3] .- 1));
f[:,3]   = f[:,3] .+ dt.*(dfin./x[:,3]);
# test
#fv_de       = (tt.*x[:,4].^(1./al) + ve_de.*x[:,3])./(tt+ve);

#dv_de = (x[:,3] - fv_de)./tt;
#ve[dv_de<0]  = ve_de[dv_de<0];
# # Fout = f[v] - outflow
# #--------------------------------------------------------------------------
fv       = (tt.*x[:,4].^(1/al) .+ ve.*x[:,3])./(tt .+ ve);
f[:,4]   = f[:,4] .+ dt.*((x[:,3] .- fv)./(tt.*x[:,4]));
f[:,5]   = f[:,5] .+ dt.*((m - fv.*x[:,5]/x[:,4])./(tt.*x[:,5]));
x[:,:]   = f[:,:];
X[t,:]   = f[:];

end


plot(1:length(u),X[:,:])
