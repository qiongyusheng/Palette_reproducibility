import numbers
from typing import Any, Mapping, Optional, TypeVar, Union

import anndata as ad
import h5py
import numpy as np
import scipy.sparse
import sklearn.neighbors

Array = Union[np.ndarray, scipy.sparse.spmatrix]
BackedArray = Union[h5py.Dataset, ad._core.sparse_dataset.SparseDataset]
AnyArray = Union[Array, BackedArray]
ArrayOrScalar = Union[np.ndarray, numbers.Number]
Kws = Optional[Mapping[str, Any]]
RandomState = Optional[Union[np.random.RandomState, int]]

T = TypeVar("T")  # Generic type var

def get_rs(x: RandomState = None) -> np.random.RandomState:
    r"""
    Get random state object

    Parameters
    ----------
    x
        Object that can be converted to a random state object

    Returns
    -------
    rs
        Random state object
    """
    if isinstance(x, int):
        return np.random.RandomState(x)
    if isinstance(x, np.random.RandomState):
        return x
    return np.random

def seurat_alignment_score(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01,
        n_repeats: int = 4, random_state: RandomState = None, **kwargs
) -> float:
    r"""
    Seurat alignment score

    Parameters
    ----------
    x
        Coordinates
    y
        Batch labels
    neighbor_frac
        Nearest neighbor fraction
    n_repeats
        Number of subsampling repeats
    random_state
        Random state
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    sas
        Seurat alignment score
    """
    rs = get_rs(random_state)
    idx_list = [np.where(y == u)[0] for u in np.unique(y)]
    min_size = min(idx.size for idx in idx_list)
    repeat_scores = []
    # ncell = np.sum([idx.size for idx in idx_list])
    # k = max(round(ncell * neighbor_frac), 1)
    k = 30
    # nn = sklearn.neighbors.NearestNeighbors(
            # n_neighbors=k + 1, **kwargs
        # ).fit(x)
    # nni = nn.kneighbors(x, return_distance=False)
    # same_y_hits = (
            # y[nni[:, 1:]] == np.expand_dims(y, axis=1)
        # ).sum(axis=1).mean()
    # repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
    
    
    for _ in range(n_repeats):
        subsample_idx = np.concatenate([
            rs.choice(idx, min_size, replace=False)
            for idx in idx_list
        ])
        subsample_x = x[subsample_idx]
        subsample_y = y[subsample_idx]
        # k = max(round(subsample_idx.size * neighbor_frac), 1)
        nn = sklearn.neighbors.NearestNeighbors(
            n_neighbors=k + 1, **kwargs
        ).fit(subsample_x)
        nni = nn.kneighbors(subsample_x, return_distance=False)
        same_y_hits = (
            subsample_y[nni[:, 1:]] == np.expand_dims(subsample_y, axis=1)
        ).sum(axis=1).mean()
        repeat_score = (k - same_y_hits) * len(idx_list) / (k * (len(idx_list) - 1))
        repeat_scores.append(min(repeat_score, 1))  # score may exceed 1, if same_y_hits is lower than expected by chance
    return np.mean(repeat_scores).item()
    # return repeat_score