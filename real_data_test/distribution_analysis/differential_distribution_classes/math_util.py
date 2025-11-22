import numpy as np
import warnings

def geom_mean(array: np.array) -> np.double:
    if len(array[array != 0]) == 0:
        return 1
    if len(array[array != 0]) == 1:
        return array[array != 0][0]
    return np.exp(np.log(array[array != 0]).mean())

def clr_transform(array: np.array) -> np.array:
    warnings.simplefilter("error")
    try:
        new_array = array.copy()
        new_array[new_array == 0] = 0.01
        return np.log(new_array/geom_mean(array))
    except:
        print(array)
        print(geom_mean(array))
        raise