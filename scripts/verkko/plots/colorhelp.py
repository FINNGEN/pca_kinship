import colorsys
import numpy as np


"""
This modules contains helper functions for obtaining sets of distinct colors.

For more functionality and colorsets see the brewer2mpl python package.
"""


def get_distinct_colors(num_colors):
    """
    Get any number of distinct colors.
    
    If num_colors <= 12, the cbrewer palettes are used.
    Otherwise the colors are generated algorithmically based on 
    the HLS color model.
    """
    assert num_colors >= 0
    if num_colors <= 12:
        try:
            return get_qualitative_brewer_colors(num_colors)
        except:
            # in case cbrewer is not available.
            pass
    return get_arbitrary_n_of_distinct_colors(num_colors)
        
        
def get_arbitrary_n_of_distinct_colors(num_colors):
    """
    Helper function for getting an arbitrary number of distinct colors.
    The algorithm is pure heuristic.
    """
    color_list = []
    hueMax = 360.
    for j, i in enumerate(np.linspace(0., hueMax, num_colors)):
        hue = i / 360.
        lightness = (38 + j % 3 * 17) / 100.
        saturation = (80 + j % 2 * 20) / 100.
        color_list.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return np.array(color_list)

def get_qualitative_brewer_colors(num_colors):
    """
    Get colors according to the colorbrewer palette.
    """
    try:
        import brewer2mpl as br
    except:
        print "Could not import brewer2mpl, do you have it in your $PYTHONPATH?"
        raise
    
    assert num_colors <= 12, "max 12 colors are supported"

    if num_colors > 9:
        brmap = br.get_map('Paired', 'Qualitative', num_colors) 
    else:
        brmap = br.get_map('Set1', 'Qualitative', num_colors)
    return np.array(brmap.colors)/255
        

"""
See the original colorbrewer palettes at:
http://colorbrewer2.org/
"""

CBREWER11 = np.array([
                     (166, 206, 227),
                    (31, 120, 180),
                    (178, 223, 138),
                    (51, 160, 44),
                    (251, 154, 153),
                    (227, 26, 28),
                    (253, 191, 111),
                    (255, 127, 0),
                    (202, 178, 214),
                    (106, 61, 154),
                    (255, 255, 153)]
                     ) / 255.

CBREWER09 = np.array([[228, 26, 28],
                      [247, 129, 191],
                      [77, 175, 74],
                      [255, 255, 51],
                      [255, 127, 0],
                      [152, 78, 163],
                      [55, 126, 184],
                      [166, 86, 40],
                      [153, 153, 153]])
