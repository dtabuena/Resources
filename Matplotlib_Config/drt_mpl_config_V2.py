""" To use copy below
import urllib
response = urllib.request.urlretrieve('https://raw.githubusercontent.com/dtabuena/Resources/main/Matplotlib_Config/drt_mpl_config_V2.py','drt_mpl_config_V2.py')
%run drt_mpl_config_V2.py
"""




version = 'v2.0'

import matplotlib.font_manager as fm
from matplotlib import rcParams
from matplotlib import pyplot as plt
import urllib
import importlib.util
import numpy as np
import os


try:
    TICK_LEN_PT = 0.5 * 72 / 25.4
    global_markersize = 1.0
    
    fig_config = {
        # LINES
        "lines.linewidth": 0.5,
        "lines.markersize": 1.0,
        "lines.markersize": global_markersize,
        "lines.markeredgewidth": 0.0,"      
    
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
        "xtick.major.size": TICK_LEN_PT,
        "xtick.major.width": 0.5,
        "xtick.minor.width": 0.5,
        "xtick.major.pad": 2,
        "ytick.labelsize": 5.25,
        "ytick.major.size": TICK_LEN_PT,
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
    print("colorspacious is not installed")
    print("try !pip install colorspacious") 

try:        

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
        
    # test #
    _ = hue_seurat(10)
   
    print('hue_seurat load success')
except:
    print('hue_seurat load failed')

def hue_seurat_cmap(h_start=25, c=80, l=60, n=256):
    hex_colors = hue_seurat(n, h_start=h_start, c=c, l=l)
    return mcolors.LinearSegmentedColormap.from_list('hue_seurat', hex_colors)



import matplotlib.colors as mcolors

def make_neutral_cmap(name, cmap, neutral_pos, neutral='lightgrey', neutral_width=0.2):
    positions = np.linspace(0.0, 1.0, 256)
    colors = cmap(positions)
    neutral_rgba = np.array(mcolors.to_rgba(neutral))
    for i, p in enumerate(positions):
        dist = abs(p - neutral_pos)
        if dist <= neutral_width:
            t = 1.0 - (dist / neutral_width)
            colors[i] = np.clip((1 - t) * colors[i] + t * neutral_rgba, 0, 1)
    result = mcolors.LinearSegmentedColormap.from_list(name, list(zip(positions, colors)), N=256)
    return result

def desaturate_cmap(name, cmap, saturation_scale=0.5, value_boost=0.1):
    positions = np.linspace(0.0, 1.0, 256)
    rgba_colors = cmap(positions)
    rgb_colors = rgba_colors[:, :3]
    hsv_colors = mcolors.rgb_to_hsv(rgb_colors)
    scaled_saturation = hsv_colors[:, 1] * saturation_scale
    boosted_value = np.clip(hsv_colors[:, 2] + value_boost, 0, 1)
    hsv_colors[:, 1] = scaled_saturation
    hsv_colors[:, 2] = boosted_value
    muted_rgb_colors = mcolors.hsv_to_rgb(hsv_colors)
    muted_rgba_colors = np.concatenate([muted_rgb_colors, rgba_colors[:, 3:4]], axis=1)
    result = mcolors.LinearSegmentedColormap.from_list(name, list(zip(positions, muted_rgba_colors)), N=256)
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
tab60 = mcolors.ListedColormap(tab60_colors, name='tab60')

piyg_grey    = make_neutral_cmap('piyg_grey',    plt.cm.PiYG,    neutral_pos=0.5)
rdbu_r_grey  = make_neutral_cmap('rdbu_r_grey',  plt.cm.RdBu_r,  neutral_pos=0.5)
rdbu_grey    = make_neutral_cmap('rdbu_grey',    plt.cm.RdBu,    neutral_pos=0.5)
Purples_grey = make_neutral_cmap('Purples_grey', plt.cm.Purples, neutral_pos=0.0)
Oranges_grey = make_neutral_cmap('Oranges_grey', plt.cm.Oranges, neutral_pos=0.0)
turbo_muted = desaturate_cmap('turbo_muted', plt.cm.turbo, saturation_scale=0.75, value_boost=0.2)

print('Custom Colors: piyg_grey, rdbu_r_grey, rdbu_grey, Purples_grey, Oranges_grey, tab60, turbo_muted')
print(f'dt_config {version}')


