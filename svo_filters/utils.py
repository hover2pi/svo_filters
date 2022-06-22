import astropy.units as q
import numpy as np


def _padded_differences(arr, pad_val=1e20):
    """

    Parameters
    ----------
    arr: array-like
    pad_val: astropy.units.quantity, float, int
        value for padding the first element of the difference array by when calculating 
        monotonically increasing (or not) differences. 
    """
    if isinstance(arr, q.Quantity):
        # Quantity units cannot be padded by non-quantity unless using this wrapper.
        diff = arr.ediff1d(to_begin=pad_val)
    else:
        diff = np.ediff1d(arr, to_begin=pad_val)
    return diff

def incremented_monotonic(arr, increment=None, increment_step=1000):
    """Returns input if monotonically increasing. Otherwise
    increment repeated elements by `increment`, which will be set to 1/`increment_step`
    of the smallest difference in array if `None`, the default.
    If not monotonically increasing (ignoring repeated elements), raises `ValueError`.

    Parameters
    ----------
    arr: array-like
        array to check for increment. Also handles astropy.units.quantity.Quantity arrays.
    increment: astropy.units.quantity, float, int (optional)
        value to increment repeated elements by. Set to 1/`increment_step` of the smallest difference 
        in array if `None`, the default. Unit conversion will be attempted if array is an instance
        of astropy.units.quantity.Quantity.
    increment_step: float, int
        The relative size difference between repeated elements if automatically determining. Only
        used if `increment = None`. 

    Returns
    -------
    array-like
        Input array if monotonically increasing else input array where repeat values
        have been incremented by `increment`.
    """   
    diff = _padded_differences(arr)

    if np.any(diff<0):
        raise ValueError("Input array must be monotonically increasing except for repeated values.")
    non_monotonic_mask = diff<=0
    
    # Exit early if monotonic
    if np.all(~non_monotonic_mask):
        return arr
    
    if increment is None:
        # a thousanth of the minimimum non-zero increment.
        increment = np.nanmin(diff[~non_monotonic_mask])/increment_step 
    # Try to help user with unit conversion, will fail if unconvertable.
    elif isinstance(arr, q.Quanity):
        increment = increment << arr.unit
        
    #non_monotonic_mask = non_monotonic_mask.reshape(arr.shape)
    repeats, multiples = get_multipliers(non_monotonic_mask)
    if isinstance(increment, q.Quantity):
        multiples = multiples << q.dimensionless_unscaled
    multiples *= increment

    flatarr = arr.flatten()
    flatarr[repeats] += multiples
    return flatarr.reshape(arr.shape)


def breadth_first(repeats, state, row):
    """somewhat convoluted 1D breadth first search 
    for getting multiples of repeated elements.
    """
    queue = [repeats[row]]
    index = 1
    additions = [index]
    state.append(row)
    while len(queue) > 0:
        idx = queue.pop(0)
        
        #print(idx, idx+1, repeats)
        if idx+1 in repeats:
            neighbor = idx+1
            state.append(row+index)
            index += 1
            queue.append(neighbor)
            additions.append(index)
            
    return additions, state


def get_multipliers(mask):
    """Get all repeats and their multiples using breadth first search

    Parameters
    ----------
    mask: array_like
        array of booleans representing repeated elements. True for repeated.
        
    Returns
    -------
    tuple 
        tuple of array_like of repeated indexes and their multiples for use in
        shifting multiple repeated indexes.
    """
    repeats = np.nonzero(mask.flatten())[0]
    if len(repeats)==0:
        raise ValueError("No repeats found. Input mask all False")
    groups = []
    state = []
    for j,val in enumerate(repeats):
        if j in state: 
            continue
        additions, state = breadth_first(repeats, state, j)
        groups.append(additions)
    groups = np.array([element for sublist in groups for element in sublist])

    return repeats, groups
   