import sys
import math
import numpy as np
import glob
import toml
import json
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Use options [t=*] and [color_max=*]

# Parameters
target_t = 1000
color_max = None
SHOW_XYLABELS = True
MULTIREPLICA = True
interpolation_type = None # Other options: None, 'antialiased', 'bilinear', 'bicubic'
ndigits = 5
dump_num_digits = 4
nf = '{:.' + str(ndigits) + 'g}'
ofilen_base0 = 'spin_structure_factor_square'
n_comp = 1
ofile_ext = '.text'
figsize = (7.5,7.5)
m = 4  # Plot for [-m, m] * np
mskip_label = 2  # Label to 0, mskip_label * pi, 2 * mskip_label * pi, ...

# Extracting simulation parameters
def extract_parameter(dict, p):
    if p in dict:
        return str(dict[p])
    else:
        return '0'

pfn = 'parameters.toml'
dict = toml.load(open(pfn))
L = int(extract_parameter(dict, 'L'))

l = int(L/2)  # Use L/2 to avoid oscillations
nk = m * l
kmax = m * np.pi
rxy = 1
dkx = kmax / nk
dky = dkx * rxy
print('L =', L, ' dkx =', dkx, ' dky =', dky)
nkx = math.ceil(kmax / dkx)
nky = math.ceil(kmax / dky)
dt = 1
nsites_ucell = 1
eps = 1e-8

J1 = extract_parameter(dict, 'J1')
J2 = extract_parameter(dict, 'J2')
J3 = extract_parameter(dict, 'J3')
J4 = extract_parameter(dict, 'J4')
Jp = extract_parameter(dict, 'Jp')
kT = extract_parameter(dict, 'kT')
param_base = '-J1_'+J1+'-J2_'+J2+'-J3_'+J3+'-J4_'+J4+'-Jp_'+Jp+'-kT_'+kT+'-L_'+str(L)

a1 = np.array([1, 0])
a2 = np.array([0, 1])

argv = sys.argv
option_ = [option for option in argv if option.startswith('t=')]
if option_:
    option_str = option_[-1]
    start = option_str.find('=') + 1
    target_t = int(option_str[start:])
    argv.remove(option_str)

color_max = None
option_ = [option for option in argv if option.startswith('color_max=')]
if option_:
    option_str = option_[-1]
    start = option_str.find('=') + 1
    color_max = float(option_str[start:])
    print('color_max was set to', color_max)
    argv.remove(option_str)

# For plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'        
# plt.rcParams['font.family']= 'sans-serif'
# plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams["font.size"] = 20
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['axes.linewidth'] = 1.2
#plt.rcParams['axes.grid']=True
#plt.rcParams['grid.linestyle']='--'
#plt.rcParams['grid.linewidth'] = 0.3
plt.rcParams["legend.markerscale"] = 1.0
plt.rcParams["legend.labelspacing"] = 0.1
plt.rcParams["legend.handlelength"] = 1.2
#plt.rcParams["legend.fancybox"] = False
#plt.rcParams["legend.framealpha"] = 1
#plt.rcParams["legend.edgecolor"] = 'black'


def index_to_coordinate(ri, L):
    x = ri % L
    y = ri // L
    return a1 * x + a2 * y

def extract_spin(spin, ri):
    if n_comp == 3:
        return np.array([spin[3*ri], spin[3*ri+1], spin[3*ri+2]])
    elif n_comp == 1:
        return np.array([spin[ri]])
    else:
        exit('n_comp ==', n_comp, 'is not supported.')

# Naive calculation: the computational and memory complexities are O(NM^2) and O(M^2), where M is the number of kx (ky) points and N is the number of sites.
def Fourier_trans(spin):
    nsites = len(spin) // n_comp
    sf = np.zeros((2*nky+1, 2*nkx+1))

    for kyi in range(0, 2*nky+1):
        ky = dky * (kyi - nky)

        for kxi in range(0, 2*nkx+1):
            kx = dkx * (kxi - nkx)

            Sk = np.zeros(n_comp, dtype=np.cdouble)
            for ri in range(nsites):
                [x, y] = index_to_coordinate(ri, L)
                # print(ri, x, y)
                phase = np.exp(1j * (kx * x + ky * y))
                Sk += extract_spin(spin, ri) * phase

            Sk /= np.sqrt(nsites)
            sfk = np.real(np.einsum('i,i', np.conj(Sk), Sk))
            sf[kyi, kxi] = sfk
            # print(kx, ky, sfk)

    return sf

# Faster calculation using a moderate memory size ~ 2MN + M^2, where M is the number of kx (ky) points and N is the number of sites.
def Fourier_trans2(spin):
    nsites = len(spin) // n_comp
    sf = np.zeros((2*nky+1, 2*nkx+1))
    phase_kx_r = np.zeros((2*nkx+1, nsites), dtype=np.cdouble)
    phase_ky_r = np.zeros((2*nky+1, nsites), dtype=np.cdouble)

    for kxi in range(0, 2*nkx+1):
        kx = dkx * (kxi - nkx)
        for ri in range(nsites):
            [x, y] = index_to_coordinate(ri, L)
            phase = np.exp(1j * (kx * x))
            phase_kx_r[kxi, ri] = phase

    for kyi in range(0, 2*nky+1):
        ky = dky * (kyi - nky)
        for ri in range(nsites):
            [x, y] = index_to_coordinate(ri, L)
            phase = np.exp(1j * (ky * y))
            phase_ky_r[kyi, ri] = phase

    spin_r = np.reshape(spin, [nsites, n_comp])

    for kyi in range(0, 2*nky+1):
        phase_kyi = phase_ky_r[kyi,:]
        for kxi in range(0, 2*nkx+1):
            phase_kxi = phase_kx_r[kxi,:]
            Sk = np.einsum('i,i,ij->j', phase_kyi, phase_kxi, spin_r)
            Sk /= np.sqrt(nsites)
            sfk = np.real(np.einsum('j,j', np.conj(Sk), Sk))
            sf[kyi, kxi] = sfk
    
    return sf


def out_structure_factor(sf, ofilen, error=None):
    vf = '{:.8g}'
    with open(ofilen, 'w') as out:
        for kyi in range(0, 2*nky+1):
            ky = dky * (kyi - nky)    
            for kxi in range(0, 2*nkx+1):
                kx = dkx * (kxi - nkx)
                if error is not None:
                    out.write(str(vf.format(kx))+' '+str(vf.format(ky))+' '+str(vf.format(sf[kyi,kxi]))+' '+str(vf.format(error[kyi,kxi]))+'\n')
                else:
                    out.write(str(vf.format(kx))+' '+str(vf.format(ky))+' '+str(vf.format(sf[kyi,kxi]))+'\n')


def read_input(ifn):
    readfile = open(ifn, 'r')
    contents = ''
    for line in readfile:
        contents = contents + line
    readfile.close()
    return contents


def plot_structure_factor(sf, ofilen_base):
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'    
    colors = ['#0628bb', '#0029fe', '#00cbfe', '#00fe00', '#f2f200', '#febc00', '#ff0900']  # in terms of rgb
    cmap_name = 'my_list'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors)
    title = r'$S(q)$'+': '+r'$L=$'+str(L)+r'$, J_1=$'+str(J1)+r'$, J_2=$'+str(J2)+r'$, J_3=$'+str(J3)+r'$, J_4=$'+str(J4)+r'$, J_p=$'+str(Jp)+r'$, T=$'+str(kT)+r'$, t=$'+str(target_t)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1,1,1)
    ax.set_title(title, fontsize=16, pad=20)
    ticks = np.arange(-kmax, kmax+eps, np.pi)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)    
    if SHOW_XYLABELS:
        ax.set_xlabel(r'$q_x$')
        ax.set_ylabel(r'$q_y$')
        ax.xaxis.set_label_coords(0.5,-0.08)
        ax.yaxis.set_label_coords(-0.08,0.5)
        labels = []
        for tick in ticks:
            coef = int(np.round(tick / np.pi))
            if coef == 0:
                labels.append(r'$0$')
            elif coef % mskip_label == 0:
                labels.append(r'${}\pi$'.format(coef))
            else:
                labels.append('')

        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    cs = ax.imshow(sf, cmap=cmap, origin='lower', extent=[-kmax, kmax, -kmax, kmax], vmin=0, vmax=color_max, interpolation=interpolation_type)
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)    
    cax = fig.add_axes([0.83, 0.2, 0.03, 0.6])
    cb = fig.colorbar(cs, cax=cax)
    cb.ax.tick_params(which='both', direction='in')
    ofign = 'plot_'+ofilen_base+'.png'
    plt.savefig(ofign, format=ofign[-3:])


def read_dump_file(dfilen, verbose=False):
    if verbose:
        print('Decoding', dfilen)    
    contents = read_input(dfilen)
    decoded = json.loads(contents)
    time   = int(decoded['n_steps'] * dt)
    spin   = decoded['spin']
    nsites = len(spin) // n_comp
    nucells = nsites // nsites_ucell
    # Assume that the system is isotropic.
    assert(nucells == L * L)
    return spin, time

def calc_structure_factor(spin, ofilen=None, output=False):   
    # sf = Fourier_trans(spin)
    sf = Fourier_trans2(spin)
    if output:
        out_structure_factor(sf, ofilen)
    return sf

def calc_plot_structure_factor(basedir):
    dumps = glob.glob(basedir+'/dump*/')
    if dumps:
        print('Number of simulations:', len(dumps))
        sf = np.zeros((2*nky+1, 2*nkx+1))
        sf_sq = np.zeros((2*nky+1, 2*nkx+1))
        n_samples = 0      
        for dump in dumps:
            dumpfiles = glob.glob(dump+'dump*.json')
            for dfile in dumpfiles:
                spin, time = read_dump_file(dfile, verbose=False)
                if int(time - target_t) == 0:
                    print('Using', dfile)
                    sfi = calc_structure_factor(spin)
                    sf += sfi
                    sf_sq += np.square(sfi)
                    n_samples += 1

        ofilen_base = ofilen_base0+param_base+'-t_'+str(target_t)
        ofilen = ofilen_base+ofile_ext

        print('n_samples =', n_samples)
        if n_samples > 0:
            mu = sf / n_samples
            if n_samples > 1:
                sigma = np.sqrt((sf_sq / n_samples - np.square(mu)) / (n_samples-1))
                out_structure_factor(mu, ofilen, error=sigma)
            else:
                out_structure_factor(mu, ofilen)
            plot_structure_factor(mu, ofilen_base)
        else:
            print('No sample was found.')


def plot_spin_structure_factor_square():
    print('target_t =', target_t)
    if MULTIREPLICA:
        calc_plot_structure_factor('dump')
    else:
        calc_plot_structure_factor('.')


if __name__ == '__main__':
    plot_spin_structure_factor_square()