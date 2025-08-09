import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.animation as animation
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap, LogNorm


# Settings
wavefile = "wavefunction.bin"
paramfile = "parameters.txt"
potential_file = "potential.csv"
flux_file = "flux.bin"

output_mp4 = "double_slit.mp4"
save_mp4 = True  

# Rendering/quality options
display_interpolation = 'bilinear'  
dpi = 500         
fps = 30

# Color maps
psi_cmap = 'inferno'
colors = ['#48267700', '#453781FF', '#404788FF', '#39468CFF', '#33638DFF', '#2D708EFF', '#238A8DFF', '#1F968BFF', '#20A389FF', '#29AF7FFF', '#3CBB75FF', '#55C667FF', '#73D055FF', '#95D480FF', '#B8DE29FF', '#DCE319FF']
positions = [0.0, 0.1, 0.165, 0.23, 0.295, 0.36, 0.425, 0.49, 0.555, 0.62, 0.685, 0.75, 0.815, 0.88, 0.945, 1.0]
pot_cmap = LinearSegmentedColormap.from_list("custom_map", list(zip(positions, colors)))

detector_linewidth = 2.0
detector_color = 'white'

# Wavefunction shown threshold for LogNorm
psi_vmin = 1e-7 


# LaTeX settigns
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12

rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Latin Modern Roman']

plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = ["Latin Modern Roman"]
plt.rcParams['axes.titlepad'] = 10
plt.rcParams['axes.labelpad'] = 10
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=LEGEND_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)


# Parameters form parameters.txt file
def parse_parameters(fname):
    params = {}
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if ':' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    key = parts[0].strip().replace(' ', '_')
                    valstr = parts[1].strip()
                    tok = valstr.split()[0]
                    try:
                        val = int(tok)
                    except ValueError:
                        try:
                            val = float(tok)
                        except ValueError:
                            val = tok
                    params[key] = val
    return params

params = parse_parameters(paramfile)

Nx = int(params.get('Nx', params.get('Nx:', 0)))
Ny = int(params.get('Ny', params.get('Ny:', 0)))
Nt = int(params.get('Nt', params.get('Nt:', 0)))
xmax = float(params.get('xmax', 0.0))
ymax = float(params.get('ymax', 0.0))
dx = float(params.get('dx', (2 * xmax) / (Nx - 1) if Nx > 1 else 1.0))
dy = float(params.get('dy', (2 * ymax) / (Ny - 1) if Ny > 1 else 1.0))
dt = float(params.get('dt', 1.0))
x_flux = float(params.get('x_flux', params.get('x_flux:', 0.0)))
Nx_flux = int(params.get('Nx_flux', params.get('Nx_flux:', 0)))

print(f"Parsed parameters: Nx={Nx}, Ny={Ny}, Nt={Nt}, xmax={xmax}, ymax={ymax}, dx={dx}, dy={dy}, dt={dt}")
print(f"Detector (x_flux) = {x_flux} (Nx_flux={Nx_flux})")

x = (-xmax) + np.arange(Nx) * dx
y = (-ymax) + np.arange(Ny) * dy
t = np.arange(Nt) * dt

# Wavefucntion
data = np.fromfile(wavefile, dtype=np.complex128)
Nt = data.size // (Nx * Ny)

print(f"Loading wavefunction: total complex entries = {data.size}, inferred Nt = {Nt}")

psi_all = data.reshape((Nt, Ny, Nx))
psi_abs2 = np.abs(psi_all) ** 2
eps = 1e-12
psi_abs2[psi_abs2 < eps] = eps
vmax = psi_abs2.max()
vmin = psi_vmin if psi_vmin is not None else psi_abs2.min()

# Potential
pot_data = np.loadtxt(potential_file, skiprows=1)
pot_x = pot_data[:, 0]
pot_y = pot_data[:, 1]
pot_Re = pot_data[:, 2]
pot_Re_grid = pot_Re.reshape((Ny, Nx))

# Flux
flux_all = np.fromfile(flux_file, dtype=np.float64)
flux_vs_t_y = flux_all.reshape((Nt, Ny))

# Axis extents
extent = [x[0], x[-1], y[0], y[-1]]


# Figure size setup
fig = plt.figure(figsize=(12, 6))
ax_left = plt.subplot2grid((1, 10), (0, 0), colspan=7)
ax_right = plt.subplot2grid((1, 10), (0, 7), colspan=3)
plt.subplots_adjust(left=0.0)

# Potential plot
x_edges = np.linspace(x[0] - dx/2, x[-1] + dx/2, Nx * + 1)
y_edges = np.linspace(y[0] - dy/2, y[-1] + dy/2, Ny * + 1)
pcm = ax_left.pcolormesh(x_edges, y_edges, pot_Re_grid, shading='auto', cmap=pot_cmap)
cbar_pot = fig.colorbar(pcm, ax=ax_left, orientation='vertical', pad=0.02, fraction=0.046)
cbar_pot.ax.set_title(r"$V_\mathrm{ds}$", fontsize=SMALL_SIZE, pad=6, loc='left') 
cbar_pot.ax.xaxis.set_label_position('bottom')
cbar_pot.ax.tick_params(labelsize=SMALL_SIZE)

# Wavefunction plot
psi0_display = psi_abs2[0]
im = ax_left.imshow(psi0_display, origin='lower', extent=extent, cmap=psi_cmap, norm=LogNorm(vmin=vmin, vmax=vmax), interpolation=display_interpolation, animated=True)
cbar_psi = fig.colorbar(im, ax=ax_left, orientation='vertical', pad=0.1, fraction=0.046, location = 'left')
cbar_psi.ax.set_title(r"$|\psi|^2$", fontsize=SMALL_SIZE, pad=6, loc='left')
cbar_psi.ax.tick_params(labelsize=SMALL_SIZE)

# Detector
ax_left.axvline(x=x_flux, color=detector_color, linewidth=detector_linewidth, zorder=5)

# Graphs labels
ax_left.set_xlabel(r"$x$")
ax_left.set_ylabel(r"$y$")
ax_left.set_xlim(-xmax,xmax)
ax_left.set_ylim(-ymax,ymax)
ax_left.xaxis.labelpad = 0
ax_left.yaxis.labelpad = 0
ax_left.set_title("Double-slit experiment simulation")

# Right plot setup: time-integrated flux vs y
line_flux, = ax_right.plot([], [], lw=2, color='blue')
#ax_right.set_xlabel("Time-integrated flux")
ax_right.set_ylim(y[0], y[-1])
#ax_right.set_yticklabels([])
ax_right.set_title("Cumulative probability flux")
ax_right.grid(False)

cumulative_flux = np.cumsum(flux_vs_t_y[:Nt, :], axis=0) * dt 
flux_xmax = np.max(np.abs(cumulative_flux))
ax_right.set_xlim(0.0, flux_xmax*1.05)

# Time
time_text = ax_left.text(0.02, 0.95, '', transform=ax_left.transAxes,
                         color='white', fontsize=14)

# Animation
def update(frame):
    psi_frame_plot = psi_abs2[frame]
    im.set_array(psi_frame_plot)
    time_text.set_text(f"t = {frame*dt:.3f}")
    cumflux = cumulative_flux[frame]
    line_flux.set_data(cumflux, y)

    return [im, line_flux, time_text]

ani = animation.FuncAnimation(fig, update, frames=Nt, interval=1000/fps, blit=True)

# Save or show
if save_mp4:
    print(f"Animating...")
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='auto'), bitrate=2000)
    ani.save(output_mp4, writer=writer, dpi=dpi)
    print("Saved.")
else:
    plt.show()
