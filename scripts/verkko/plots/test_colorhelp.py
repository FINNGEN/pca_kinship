import colorhelp as ch
from matplotlib.colors import ColorConverter


def test():
    n = 1
    colors = ch.get_distinct_colors(n)
    assert len(colors) is 1
    n = 10
    colors = ch.get_distinct_colors(n)
    assert len(colors) is 10
    colors = ch.get_arbitrary_n_of_distinct_colors(11)
    assert len(colors) is 11
    colors = ch.get_qualitative_brewer_colors(12)
    assert len(colors) is 12
    colors2 = ch.get_qualitative_brewer_colors(10)
    # same first colors should be returned with 10-12 colors used:
    assert (colors[9] == colors2[-1]).all()
    
    
    cc = ColorConverter()
    for color in colors:
        assert len(cc.to_rgba(color)) is 4
