from secular_pt_fun import *
from secular_pt_param import *

data_dir = '/mnt/home/kvantilburg/ceph/GravitationalBasin/sec_pt/'

n_a = int(1e3)
a_min = 0.4
a_max = 10.
vec_a = np.arange(a_min,a_max,(a_max-a_min)/n_a)
prob_a = vec_a**(-2)/np.sum(vec_a**(-2)); # pdf that scales as f(a) \propto 1/a^2

n_e = int(1e3)
e_min = 0
e_max = 1
vec_e = np.arange(e_min,e_max,(e_max-e_min)/n_e)
prob_e = vec_e / np.sum(vec_e) # pdf scales as f(e) \propto e

n_inc = int(1e3)
inc_min = 0
inc_max = np.pi
vec_inc = np.arange(inc_min,inc_max,(inc_max-inc_min)/n_inc)
prob_inc = np.sin(vec_inc) / np.sum(np.sin(vec_inc))

n_sample = 2**10 # ideally power of 2 for FFT purposes
pts_xy_Earth_orbit = fn_gen_pts_xyz(df_planets[df_planets['name']=='Earth']['a'].to_numpy()[0],
                                    df_planets[df_planets['name']=='Earth']['e'].to_numpy()[0],
                                    n_points=n_sample)[:2,:]

n_orbits = int(1e6)
n_saves = int(1e3)
n_t_0 = int(1e3)
n_points = int(1e3)
n_bins_xy = int(1e2); n_bins_z = int(1e2);
lim_xy = 1.2; lim_z = 1.2;
bins_xy = np.linspace(-lim_xy,lim_xy,n_bins_xy);
bins_z = np.linspace(-lim_z,lim_z,n_bins_z);
t = 0

for l in range(n_orbits):
    a_0 = vec_a[np.random.choice(n_a,1,p=prob_a)]
    e_0 = vec_e[np.random.choice(n_e,1,p=prob_e)]
    inc_0 = vec_inc[np.random.choice(n_inc,1,p=prob_inc)]
    t_0 = np.random.random(n_t_0)*(-4.6e9*yr) #make these random to avoid aliasing
    omega_0 = np.random.random()*2*np.pi #random omega
    Omega_0 = np.random.random()*2*np.pi #random Omega
    pts_xyz_E = fn_gen_pts_xyz_E(e_0,inc_0,omega_0,Omega_0,a_0,t_0,t,df_planets,
                                 e_il,inc_il,g_l,f_l,beta_l,gamma_l,n_points)
    if np.mod(l,int(5))==0:
        print('l = '+str(l))
        
    if l == 0:
        hist = np.histogramdd(pts_xyz_E,bins=[bins_xy,bins_xy,bins_z])[0]
    elif np.mod(l,n_orbits//n_saves)==0:
        j = str(int(np.random.rand()*1e9))
        np.save(data_dir+'hist_strap_'+str(j)+'.npy',hist)
        print(j + ' saved!')
        hist = np.histogramdd(pts_xyz_E,bins=[bins_xy,bins_xy,bins_z])[0]
    else:
        hist += np.histogramdd(pts_xyz_E,bins=[bins_xy,bins_xy,bins_z])[0]
j = str(int(np.random.rand()*1e9))
np.save(data_dir+'hist_strap_'+str(j)+'.npy',hist)

print('Finished')