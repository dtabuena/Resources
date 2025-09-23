""" To use copy below
import urllib
response = urllib.request.urlretrieve('https://raw.githubusercontent.com/dtabuena/Resources/main/Matplotlib_Config/Load_FS6.py','Load_FS6.py')
%run Load_FS6.py
"""



import matplotlib.font_manager as fm
from matplotlib import rcParams
import urllib
fig_config =    {
    "font.size": 6,
    "font.family": "arial",
    "axes.linewidth": 0.5,
    "lines.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "figure.titlesize": 6,
    "axes.titlesize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "axes.labelsize": 6,
    "legend.fontsize": 6,
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "grid.color": "black",
    "grid.linestyle": "-",
    "grid.linewidth": 0.1,
    "figure.figsize": [1.5,1.5],
    "svg.fonttype": "none",
    "savefig.bbox": 'tight',
    "savefig.transparent": True,
    'lines.markersize': 1.0,
    'boxplot.meanprops.markersize': 1.0,
    }
_ = urllib.request.urlretrieve('https://github.com/dtabuena/Resources/raw/main/Fonts/arial.ttf','arial.ttf')
fm.fontManager.addfont('./arial.ttf')
rcParams.update(fig_config)



##################################################
####### Defining a Seurat-like Colorscheme #######
##################################################

import numpy as np
from colorspacious import cspace_convert

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
    return np.clip(chroma_ratio, 0.5, 1.0)

def apply_farver_transform(colors, input_format='hex'):
    """Apply farver chroma reduction to a list of colors"""
    
    transformed_colors = []
    
    for color in colors:
        
        if input_format == 'hex':
            hex_clean = color.lstrip('#')
            rgb = [int(hex_clean[i:i+2], 16)/255.0 for i in (0, 2, 4)]
            hcl = cspace_convert([rgb], "sRGB1", "CIELCh")[0]
            h, c, l = hcl[1], hcl[0], hcl[2]
            
        elif input_format == 'hcl':
            h, c, l = color
            
        elif input_format == 'hsl':
            h_hsl, s, l_hsl = color
            import colorsys
            rgb = colorsys.hls_to_rgb(h_hsl/360, l_hsl/100, s/100)
            hcl = cspace_convert([rgb], "sRGB1", "CIELCh")[0]
            h, c, l = hcl[1], hcl[0], hcl[2]
            
        else:
            raise ValueError("input_format must be 'hex', 'hcl', or 'hsl'")
        
        chroma_ratio = apply_farver_chroma_reduction(h)
        adjusted_c = c * chroma_ratio
        
        lch = np.array([[l, adjusted_c, h]])
        rgb_out = cspace_convert(lch, "CIELCh", "sRGB1")[0]
        
        hex_out = "#{:02x}{:02x}{:02x}".format(
            int(rgb_out[0]*255), int(rgb_out[1]*255), int(rgb_out[2]*255)
        )
        transformed_colors.append(hex_out)
    
    return transformed_colors

def hue_seurat(n_clusters, h_start=10, c=100, l=65):
    """Generate farver-corrected hue palette for scanpy"""
    
    if n_clusters == 0:
        raise ValueError("Must request at least one color")
    
    # Generate evenly spaced hues around color wheel
    hue_step = 360 / n_clusters
    hues = [(h_start + i * hue_step) % 360 for i in range(n_clusters)]
    
    # Create HCL colors and apply farver transform
    hcl_colors = [(h, c, l) for h in hues]
    hex_colors = apply_farver_transform(hcl_colors, 'hcl')
    
    return hex_colors

