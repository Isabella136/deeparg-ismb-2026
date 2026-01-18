import numpy as np
import warnings

def geom_mean(array: np.array) -> np.double:
    if len(array[array != 0]) == 0:
        return 1
    if len(array[array != 0]) == 1:
        return array[array != 0][0]
    return np.exp(np.log(array[array != 0]).mean())

def array_wise_clr_transform(array: np.array) -> np.array:
    warnings.simplefilter("error")
    try:
        return_array = array.copy()
        return_array[array == 0] = -1 * np.inf
        return_array[array != 0] = np.log(array[array != 0]/geom_mean(array))
        return return_array
    except:
        print(array)
        print(geom_mean(array))
        raise

def relative_abundance(array: np.array) -> np.array:
    return array / np.sum(array)