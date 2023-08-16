import os
import os.path
import sys
import numpy as np
import glob
import json
import toml
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

# Input
argv = sys.argv
argv.pop(0)
argc = len(argv)
sim_idx = []
if argc > 0:
    sim_idx = argv[0:]

# Plot parameters
output_format = 'png'
obsname = 'spin'
print('obsname =', obsname)
# combine_panels = False
combine_panels = True
ncols = 5
nrows = 6
npanels = ncols * nrows
dump_file_idx = range(0,npanels)
fig_scale = 30
font_size = fig_scale * 0.9
font_size_sup = fig_scale * 0.9
ndigits = 5
nf = '{:.' + str(ndigits) + 'g}'
ndigits2 = 3
nf2 = '{:.' + str(ndigits2) + 'f}'
pfn = '../../parameters.toml'
op_symbol = r'$\phi_{2 \times 1}$'

# Plot
plt.rcParams['font.family']= 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams["font.size"] = 16
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['axes.linewidth'] = 1.2
# plt.rcParams["legend.markerscale"] = 1.2

def set_title(time=None):
    params = toml.load(pfn)
    L = params['L']
    kT = params['kT']
    J1 = params['J1']
    J2 = params['J2']
    J3 = params['J3']
    Jp = params['Jp']
    title = 'Spin: '+r'$L=$'+str(L)+', '+r'$J_1=$'+str(J1)+', '+r'$J_2=$'+str(J2)+', '+r'$J_3=$'+str(J3)+', '+r'$J_p=$'+str(Jp)+', '+r'$T=$'+str(nf.format(kT))
    if time is not None:
        title += ', '+r'$t=$'+str(time)
    return title

def plot_snapshot(ifilen, fig=None, figidx=None):
    print('Processing: ', ifilen)
    ifile = open(ifilen, 'r')
    contents = ""
    for line in ifile:
      contents = contents + line
    ifile.close()
    decoded = json.loads(contents)
    dt = 1
    time   = int(decoded['n_steps'] * dt)
    if 'order_param' in decoded:
        op = decoded['order_param']
    else:
        op = None
    spin   = decoded[obsname]
    nsites = len(spin)
    print('time = ', time)

    params = toml.load(pfn)
    L = params['L']
    lx = L
    # ly = L

    a1 = np.array([1.0, 0])
    a2 = np.array([0, 1.0])
    def ucell_cod(i,j):
        return i*a1 + j*a2
    def s_cod(i):
        yi = i // lx
        xi = i % lx
        return ucell_cod(xi,yi)

    x      = np.empty([nsites])
    y      = np.empty([nsites])
    sz     = np.empty((L,L))

    for i in range(nsites):
        r = s_cod(i)        
        x[i] = r[0]
        y[i] = r[1]
        xi = i % lx
        yi = i // lx
        sz[yi, xi] = spin[i]

    # m_z = np.sum(sz) / nsites

    if combine_panels:
        ax = fig.add_subplot(nrows, ncols, figidx+1)
        title = r'$t=$'+str(time)
        if op is not None:
            title += '\n'+op_symbol+r'$=$'+str(nf2.format(op))
        ax.set_title(title, fontsize=font_size)        
    else:
        fig = plt.figure(figsize=(fig_scale, fig_scale))
        ax = fig.add_subplot()
        ax.set_title(set_title(time=time), fontsize=font_size)

    # ax.axis('off')    
    ax.set_xticklabels('')
    ax.set_xticks([])
    ax.set_yticklabels('')
    ax.set_yticks([]) 

    # sidx = np.arange(0, L, 1)
    # X, Y = np.meshgrid(sidx, sidx)
    colors = ['#0628bb', '#0029fe', '#00cbfe', '#00fe00', '#f2f200', '#febc00', '#fe0900']  # in terms of rgb         
    cmap_name = 'my_list'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors)    
    cs = ax.matshow(sz, cmap=cmap, vmin=-1, vmax=1)

    if not combine_panels:
        fign =obsname+'-snapshot'+'-t'+str(time)+'.'+output_format
        print('Output:', fign)
        plt.savefig(fign, bbox_inches='tight', format=output_format)

if combine_panels:
    fig = plt.figure(figsize=(fig_scale,fig_scale))
else:
    fig = None


dumps = glob.glob('dump/dump*/')
if dumps:
    if sim_idx:
        dumps = []        
        for sim_i in sim_idx:
            dump_i = 'dump/dump'+sim_i
            if os.path.exists(dump_i):
                dumps.append(dump_i)

        for dump in dumps:
            sim_idx = dump[9:]
            print('sim_idx =', sim_idx)
            os.chdir(dump)
            for figidx, dump_file_i in enumerate(dump_file_idx):
                i2 = str(dump_file_i).zfill(4)
                dfilen = 'dump'+str(i2)+'.json'
                plot_snapshot(dfilen, fig, figidx)

            if combine_panels:
                fign =obsname+'-snapshot'+'.'+output_format
                print('Output:', fign)
                plt.suptitle(set_title(), y=0.994, fontsize=font_size_sup)
                plt.tight_layout(pad=0.8)           
                plt.savefig(fign, bbox_inches='tight', format=output_format)
                plt.clf()

            os.chdir('../..')

plt.clf()
plt.close()