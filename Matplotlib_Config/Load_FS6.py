""" To use copy below
import urllib
response = urllib.request.urlretrieve('https://raw.githubusercontent.com/dtabuena/Resources/main/Matplotlib_Config/Load_FS6.py','Load_FS6.py')
%run Load_FS6.py
"""

import matplotlib.font_manager as fm
from matplotlib import rcParams
import urllib

try:
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
        'xtick.major.size': 2,
        'xtick.major.pad': 2,
        'ytick.major.size': 2,
        'ytick.major.pad': 2,
        'legend.handlelength': 1.0,
        'legend.handleheight': 0.7,
        'legend.markerscale': 0.5,
        }
    _ = urllib.request.urlretrieve('https://github.com/dtabuena/Resources/raw/main/Fonts/arial.ttf','arial.ttf')
    fm.fontManager.addfont('./arial.ttf')
    rcParams.update(fig_config)
    
    print('Matplotlib_config load success')
except:
    print('Matplotlib_config load failed')   

try:  
    ##################################################
    ####### Defining a Seurat-like Colorscheme #######
    ##################################################

   
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
        
    print('hue_seurat load success')
except:
    print('hue_seurat load failed')
