import numpy as np
cimport numpy as cnp

cnp.import_array()

DTYPE_t = np.float64
ctypedef cnp.float64_t DTYPE_ct

cpdef cython_predict_proba(DTYPE_ct[:, :] X_np,
                           DTYPE_ct[:, :] lp_prediction_np,
                           list reverse_combinations,  # Adjusted type here
                           unsigned long long int label_count):

    cdef unsigned long long int rs = X_np.shape[0]
    cdef unsigned long long int cs = label_count
    cdef DTYPE_ct[:, :] result = np.zeros((rs, cs), dtype=DTYPE_t)

    cdef int index, combination_id, label
    cdef cnp.ndarray label_array  # For accessing inner numpy arrays

    for index in range(X_np.shape[0]):
        for combination_id in range(lp_prediction_np.shape[1]):
            label_array = reverse_combinations[combination_id]
            for label in label_array:  # Iterate over inner numpy array
                result[index, label] += lp_prediction_np[index, combination_id]
    return np.asarray(result)
