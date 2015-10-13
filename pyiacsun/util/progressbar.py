__all__ = ['progressbar']

import sys

def progressbar(current, total, width=30, end=False):
    """Progress bar
    
    Args:
        current (float): current value of the bar
        total (float): total of the bar
        width (int, optional): number of spaces of the bar
        end (bool, optional): end character
    
    Returns:
        None: None
    """
    bar_width = width
    block = int(round(bar_width * current/total))
    text = "\rProgress: [{0}] {1} of {2}".\
        format("#"*block + "-"*(bar_width-block), current, total)
    if end:
        text = text +'\n'
    sys.stdout.write(text)
    sys.stdout.flush()
