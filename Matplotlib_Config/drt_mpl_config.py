"""
To use copy below:

import urllib
response = urllib.request.urlretrieve('https://raw.githubusercontent.com/dtabuena/Resources/main/Matplotlib_Config/drt_mpl_config.py','drt_mpl_config.py')
%run drt_mpl_config.py

"""

version = 'v2.1'
print(f'dt_config {version}')



import matplotlib.font_manager as fm
from matplotlib import rcParams
from matplotlib import pyplot as plt
import matplotlib as mpl
import urllib
import importlib.util
import numpy as np
import os


try:
    tick_len_mm = 0.5
    tick_len_pts = tick_len_mm * 2.8346 # pts is the length unit 'points' not pixels
    global_markersize = 1.0
    
    fig_config = {
        # LINES
        "lines.linewidth": 0.5,
        "lines.markersize": 1.0,
        "lines.markersize": global_markersize,
        "lines.markeredgewidth": 0.0,     
    
        # PATCHES
        "patch.linewidth": 0.5,
    
        # BOXPLOT
        "boxplot.meanprops.markersize": global_markersize*.8,
          
    
        # FONT
        "font.size": 6,
        "font.family": "arial",
    
        # AXES
        "axes.linewidth": 0.5,
        "axes.labelsize": 6,
        "axes.titlesize": 6,
    
        # TICKS
        "xtick.labelsize": 5.25,
        "xtick.major.size": tick_len_pts,
        "xtick.major.width": 0.5,
        "xtick.minor.width": 0.5,
        "xtick.major.pad": 2,
        "ytick.labelsize": 5.25,
        "ytick.major.size": tick_len_pts,
        "ytick.major.width": 0.5,
        "ytick.minor.width": 0.5,
        "ytick.major.pad": 2,
    
        # GRIDS
        "grid.color": "grey",
        "grid.linestyle": "-",
        "grid.linewidth": 0.1,
    
        # LEGEND
        "legend.fontsize": 5.25,
        "legend.handlelength": 0.5,
        "legend.handleheight": 0.5,
        "legend.markerscale": 1,
        "legend.handletextpad": 0.25,
        "legend.borderaxespad": 0,
        "legend.borderpad": 0.2,
        "legend.labelspacing": 0.2,
        "legend.columnspacing": 0.5,
    
        # FIGURE
        "figure.titlesize": 6,
        "figure.figsize": [1.5, 1.5],
        "figure.dpi": 300,
    
        # SAVING FIGURES
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.transparent": True,
        "svg.fonttype": "none",
    }
    
    rcParams.update(fig_config)
    
    print('Matplotlib_config load success')
except:
    print('Matplotlib_config load failed')   

##################################################
####### Defining a Seurat-like Colorscheme #######
##################################################

# Check if colorspacious is installed
spec = importlib.util.find_spec("colorspacious")
if spec is not None:
    import colorspacious    
    from colorspacious import cspace_convert
    print(f"colorspacious is installed, version: {colorspacious.__version__}")
else:
    print("colorspacious is not installed; optional '!pip install colorspacious'")


def apply_farver_chroma_reduction(hue):
    """Apply farver chroma reduction model - R² = 0.9846"""
    
    # Fitted Gaussian parameters: center, width, depth
    c1, w1, d1 = 82.6, 21.2, 0.295   # Yellow region
    c2, w2, d2 = 188.3, 28.1, 0.461  # Cyan region  
    c3, w3, d3 = 270.1, 15.0, 0.058  # Blue region
    
    # Calculate circular distances
    def circular_distance(h1, h2):
        diff = abs(h1 - h2)
        return min(diff, 360 - diff)
    
    # Apply three Gaussians
    d1_dist = circular_distance(hue, c1)
    d2_dist = circular_distance(hue, c2) 
    d3_dist = circular_distance(hue, c3)
    
    gauss1 = d1 * np.exp(-(d1_dist**2) / (2 * w1**2))
    gauss2 = d2 * np.exp(-(d2_dist**2) / (2 * w2**2))
    gauss3 = d3 * np.exp(-(d3_dist**2) / (2 * w3**2))
    
    chroma_ratio = 1.0 - gauss1 - gauss2 - gauss3
    return np.clip(chroma_ratio, 0.0, 1.0)

def apply_farver_transform(hcl_colors):
    """Apply farver chroma reduction to HCL colors"""
    import warnings
    
    transformed_colors = []
    
    for h, c, l in hcl_colors:
        chroma_ratio = apply_farver_chroma_reduction(h)
        adjusted_c = c * chroma_ratio
        
        lch = np.array([[l, adjusted_c, h]])
        
        # Suppress colorspacious warnings about out-of-gamut colors
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rgb_out = cspace_convert(lch, "CIELCh", "sRGB1")[0]
        
        rgb_out = np.clip(rgb_out, 0, 1)
        hex_out = "#{:02x}{:02x}{:02x}".format(
            int(rgb_out[0]*255), int(rgb_out[1]*255), int(rgb_out[2]*255)
        )
        transformed_colors.append(hex_out)
    
    return transformed_colors

def hue_seurat(n_clusters, h_start=25,c=80,l=60):
    """Generate farver-corrected hue palette for scanpy"""
    
    if n_clusters == 0:
        raise ValueError("Must request at least one color")
    
    # Generate evenly spaced hues around color wheel
    hue_step = 360 / n_clusters
    hues = [(h_start + i * hue_step) % 360 for i in range(n_clusters)]
    
    # Create HCL colors and apply farver transform
    hcl_colors = [(h, c, l) for h in hues]
    hex_colors = apply_farver_transform(hcl_colors)
    
    return hex_colors


# Smoke test only runs if the dependency the functions need is present.
if spec is not None:
    _ = hue_seurat(10)
    print('hue_seurat load success')
else:
    print('hue_seurat defined but needs colorspacious to run')






# import matplotlib.colors as mpl.colors
def make_sequential_cmap(name, color_list, n_levels=256):
    if len(color_list) < 2:
        raise ValueError(f"A sequential ramp needs at least two colors, got {len(color_list)}.")
    positions = np.linspace(0.0, 1.0, n_levels)
    stops = np.linspace(0.0, 1.0, len(color_list))
    stop_colors = np.array([mpl.colors.to_rgba(entry) for entry in color_list])
    colors = np.empty((n_levels, 4))
    for channel in range(stop_colors.shape[1]):
        colors[:, channel] = np.interp(positions, stops, stop_colors[:, channel])
    result = mpl.colors.LinearSegmentedColormap.from_list(name, list(zip(positions, colors)), N=n_levels)
    return result



def make_neutral_cmap(name, cmap, neutral_pos, neutral='lightgrey', neutral_width=0.2):
    positions = np.linspace(0.0, 1.0, 256)
    colors = cmap(positions)
    neutral_rgba = np.array(mpl.colors.to_rgba(neutral))
    for i, p in enumerate(positions):
        dist = abs(p - neutral_pos)
        if dist <= neutral_width:
            t = 1.0 - (dist / neutral_width)
            colors[i] = np.clip((1 - t) * colors[i] + t * neutral_rgba, 0, 1)
    result = mpl.colors.LinearSegmentedColormap.from_list(name, list(zip(positions, colors)), N=256)
    return result

def desaturate_cmap(name, cmap, saturation_scale=0.5, value_boost=0.1):
    positions = np.linspace(0.0, 1.0, 256)
    rgba_colors = cmap(positions)
    rgb_colors = rgba_colors[:, :3]
    hsv_colors = mpl.colors.rgb_to_hsv(rgb_colors)
    scaled_saturation = hsv_colors[:, 1] * saturation_scale
    boosted_value = np.clip(hsv_colors[:, 2] + value_boost, 0, 1)
    hsv_colors[:, 1] = scaled_saturation
    hsv_colors[:, 2] = boosted_value
    muted_rgb_colors = mpl.colors.hsv_to_rgb(hsv_colors)
    muted_rgba_colors = np.concatenate([muted_rgb_colors, rgba_colors[:, 3:4]], axis=1)
    result = mpl.colors.LinearSegmentedColormap.from_list(name, list(zip(positions, muted_rgba_colors)), N=256)
    return result


def dt_save_fig(fig_obj, save_name_loc, formats=('jpeg','svg')):
    """
    Save one figure to the same path stem in several formats.

    fig_obj : matplotlib Figure to write.
    save_name_loc : path stem, directory plus base file name, without an extension.
        The directory must already exist; savefig does not create it.
    formats : file extensions to write, one file each.

    Returns the list of paths written, in the order of formats.
    """
    saved = []
    for f in formats:
        final_file = f'{save_name_loc}.{f}'
        fig_obj.savefig(final_file, bbox_inches='tight')
        saved.append(final_file)
    return saved

def show_legend_figure(source_figure, number_of_columns=None, legend_title=None):
    handle_by_label = {}
    for source_axis in source_figure.axes:
        axis_handles, axis_labels = source_axis.get_legend_handles_labels()
        for axis_handle, axis_label in zip(axis_handles, axis_labels):
            if axis_label not in handle_by_label:
                handle_by_label[axis_label] = axis_handle
    if len(handle_by_label) == 0:
        raise ValueError("The figure has no labeled artists; there is nothing to put in a legend.")
    if number_of_columns is None:
        number_of_columns = len(handle_by_label)

    legend_margin_factor = 1.02
    legend_figure = plt.figure(figsize=(1, 1))
    legend_object = legend_figure.legend(
        handles=list(handle_by_label.values()),
        labels=list(handle_by_label.keys()),
        loc="center",
        ncols=number_of_columns,
        title=legend_title,
    )
    legend_figure.canvas.draw()
    legend_extent = legend_object.get_window_extent()
    legend_width_inches = legend_margin_factor * legend_extent.width / legend_figure.dpi
    legend_height_inches = legend_margin_factor * legend_extent.height / legend_figure.dpi
    legend_figure.set_size_inches(legend_width_inches, legend_height_inches)
    plt.show()
    return legend_figure

"""
Some Fonts
"""
font_list = ['Helvetica','arial','Helvetica-Bold','Helvetica-BoldOblique','Helvetica-Oblique']
for font in font_list:
    try:
        _ = urllib.request.urlretrieve(f'https://github.com/dtabuena/Resources/raw/main/Fonts/{font}.ttf',f'{font}.ttf')
        fm.fontManager.addfont(f'./{font}.ttf')
    except:
        print(f"{font} failed")



# Tab60
tab60_colors = (list(plt.cm.tab20.colors)
                + list(plt.cm.tab20b.colors)
                + list(plt.cm.tab20c.colors))
tab60 = mpl.colors.ListedColormap(tab60_colors, name='tab60')

PiYG_grey    = make_neutral_cmap('PiYG_grey',    plt.cm.PiYG,    neutral_pos=0.5)
# rdbu_r_grey  = make_neutral_cmap('rdbu_r_grey',  plt.cm.RdBu_r,  neutral_pos=0.5)
# rdbu_grey    = make_neutral_cmap('rdbu_grey',    plt.cm.RdBu,    neutral_pos=0.5)
Purples_grey = make_neutral_cmap('Purples_grey', plt.cm.Purples, neutral_pos=0.0)
Oranges_grey = make_neutral_cmap('Oranges_grey', plt.cm.Oranges, neutral_pos=0.0)
turbo_muted = desaturate_cmap('turbo_muted', plt.cm.turbo, saturation_scale=0.75, value_boost=0.2)
fn_map = make_sequential_cmap('fn_map', ['#B3995E', '#FFFFFF', '#AA0000'])

rdbu_r_grey = 'rdbu_r_grey has been depricated, use coolwarm'
rdbu_grey = 'rdbu_grey has been depricated, use coolwarm'
piyg_grey = 'piyg_grey has been depricated, use PiYG_grey'

print('display_dt_customs() to display maps')


def display_dt_customs(custom_maps = [PiYG_grey,Purples_grey,Oranges_grey,turbo_muted,tab60,fn_map]):
    for cm in custom_maps:
        display(cm)
    return None

def plot_ruler(size_inches=(4, .75), text_sizes=(4.0, 5.25, 6.0, 8.0),
               marker_sizes=(1, 2, 3, 4, 6), scale_bar_inches=1.0,
               caption_fontsize=5.25, to_save=True,return_fig=False):
    """Size reference: text and marker samples, a one-inch bar, and the live rcParams.

    The border is on the figure boundary, so a display or save that crops or rescales
    the canvas shows it. Compare the bar against a real ruler to check the display.
    """
    figure, axes = plt.subplots(1, 2, figsize=size_inches, layout='constrained',
                                width_ratios=(len(text_sizes), len(marker_sizes)))
    for column, size_value in enumerate(text_sizes):
        axes[0].text((column + 0.5) / len(text_sizes), 0.80, f'Agj {size_value}',
                     fontsize=size_value, ha='center', va='center',
                     transform=axes[0].transAxes)
    axes[0].text(0.5, 0.46, 'font size (pt)', fontsize=caption_fontsize,
                 ha='center', va='center', transform=axes[0].transAxes)

    for column, marker_size in enumerate(marker_sizes):
        centre = (column + 0.5) / len(marker_sizes)
        axes[1].plot(centre, 0.84, marker='o', markersize=marker_size, color='0.2',
                     transform=axes[1].transAxes)
        axes[1].text(centre, 0.62, f'{marker_size}', fontsize=caption_fontsize,
                     ha='center', va='center', transform=axes[1].transAxes)
    axes[1].text(0.5, 0.46, 'markersize (pt)', fontsize=caption_fontsize,
                 ha='center', va='center', transform=axes[1].transAxes)
    for axis in axes:
        axis.set_axis_off()

    bar_fraction = scale_bar_inches / size_inches[0]
    figure.add_artist(mpl.lines.Line2D(
        [0.5 - bar_fraction / 2, 0.5 + bar_fraction / 2], [0.26, 0.26],
        transform=figure.transFigure, color='k', lw=0.8))
    figure.text(0.5, 0.29, f'{scale_bar_inches:g} inch', transform=figure.transFigure,
                ha='center', va='bottom', fontsize=caption_fontsize)
    figure.text(0.5, 0.06,
                f'figsize {size_inches[0]:g} x {size_inches[1]:g} in   '
                f'figure.dpi {mpl.rcParams["figure.dpi"]:g}   '
                f'savefig.dpi {mpl.rcParams["savefig.dpi"]:g}   '
                f'font.size {mpl.rcParams["font.size"]:g}',
                transform=figure.transFigure, ha='center', va='bottom',
                fontsize=caption_fontsize - 1.25)
    figure.add_artist(mpl.patches.Rectangle(
        (0.0, 0.0), 1.0, 1.0, transform=figure.transFigure, fill=False,
        edgecolor='k', linewidth=1.2, clip_on=False))
    
    if not to_save:
        dt_save_fig(figure,'ruler')
    if return_fig:
        return figure
    else:
        return None

