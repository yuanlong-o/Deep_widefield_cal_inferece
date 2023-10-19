#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" functions that creates image from a video file

Primarily intended for plotting, returns correlation images ( local or max )

See Also:
------------

@author andrea giovannucci
"""

from builtins import range

import cv2
import logging
import numpy as np
from scipy.ndimage import convolve, generate_binary_structure
from scipy.sparse import coo_matrix
from typing import Any, List, Optional, Tuple
import numbers



def old_div(a, b):
    """
    Equivalent to ``a / b`` on Python 2 without ``from __future__ import
    division``.

    TODO: generalize this to other objects (like arrays etc.)
    """
    if isinstance(a, numbers.Integral) and isinstance(b, numbers.Integral):
        return a // b
    else:
        return a / b

__all__ = ['PY3', 'PY2', 'PYPY', 'with_metaclass', 'native', 'old_div']


def mean_psd(y, method='logmexp'):
    """
    Averaging the PSD

    Args:
        y: np.ndarray
             PSD values

        method: string
            method of averaging the noise.
            Choices:
             'mean': Mean
             'median': Median
             'logmexp': Exponential of the mean of the logarithm of PSD (default)

    Returns:
        mp: array
            mean psd
    """

    if method == 'mean':
        mp = np.sqrt(np.mean(old_div(y, 2), axis=-1))
    elif method == 'median':
        mp = np.sqrt(np.median(old_div(y, 2), axis=-1))
    else:
        mp = np.log(old_div((y + 1e-10), 2))
        mp = np.mean(mp, axis=-1)
        mp = np.exp(mp)
        mp = np.sqrt(mp)

    return mp


def get_noise_fft(Y, noise_range=[0.25, 0.5], noise_method='logmexp', max_num_samples_fft=3072,
                  opencv=True):
    """Estimate the noise level for each pixel by averaging the power spectral density.

    Args:
        Y: np.ndarray
            Input movie data with time in the last axis

        noise_range: np.ndarray [2 x 1] between 0 and 0.5
            Range of frequencies compared to Nyquist rate over which the power spectrum is averaged
            default: [0.25,0.5]

        noise method: string
            method of averaging the noise.
            Choices:
                'mean': Mean
                'median': Median
                'logmexp': Exponential of the mean of the logarithm of PSD (default)

    Returns:
        sn: np.ndarray
            Noise level for each pixel
    """
    T = Y.shape[-1]
    # Y=np.array(Y,dtype=np.float64)

    if T > max_num_samples_fft:
        Y = np.concatenate((Y[..., 1:max_num_samples_fft // 3 + 1],
                            Y[..., np.int(T // 2 - max_num_samples_fft / 3 / 2)
                                          :np.int(T // 2 + max_num_samples_fft / 3 / 2)],
                            Y[..., -max_num_samples_fft // 3:]), axis=-1)
        T = np.shape(Y)[-1]

    # we create a map of what is the noise on the FFT space
    ff = np.arange(0, 0.5 + 1. / T, 1. / T)
    ind1 = ff > noise_range[0]
    ind2 = ff <= noise_range[1]
    ind = np.logical_and(ind1, ind2)
    # we compute the mean of the noise spectral density s
    if Y.ndim > 1:
        if opencv:
            import cv2
            try:
                cv2.setNumThreads(0)
            except:
                pass
            psdx_list = []
            for y in Y.reshape(-1, T):
                dft = cv2.dft(y, flags=cv2.DFT_COMPLEX_OUTPUT).squeeze()[
                    :len(ind)][ind]
                psdx_list.append(np.sum(1. / T * dft * dft, 1))
            psdx = np.reshape(psdx_list, Y.shape[:-1] + (-1,))
        else:
            xdft = np.fft.rfft(Y, axis=-1)
            xdft = xdft[..., ind[:xdft.shape[-1]]]
            psdx = 1. / T * abs(xdft)**2
        psdx *= 2
        sn = mean_psd(psdx, method=noise_method)

    else:
        xdft = np.fliplr(np.fft.rfft(Y))
        psdx = 1. / T * (xdft**2)
        psdx[1:] *= 2
        sn = mean_psd(psdx[ind[:psdx.shape[0]]], method=noise_method)

    return sn, psdx


def max_correlation_image(Y, bin_size: int = 1000, eight_neighbours: bool = True, swap_dim: bool = True) -> np.ndarray:
    """Computes the max-correlation image for the input dataset Y with bin_size

    Args:
        Y:  np.ndarray (3D or 4D)
            Input movie data in 3D or 4D format

        bin_size: scalar (integer)
             Length of bin_size (if last bin is smaller than bin_size < 2 bin_size is increased to impose uniform bins)

        eight_neighbours: Boolean
            Use 8 neighbors if true, and 4 if false for 3D data (default = True)
            Use 6 neighbors for 4D data, irrespectively

        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front

    Returns:
        Cn: d1 x d2 [x d3] matrix,
            max correlation image
    """

    if swap_dim:
        Y = np.transpose(Y, tuple(np.hstack((Y.ndim - 1, list(range(Y.ndim))[:-1]))))

    T = Y.shape[0]
    if T <= bin_size:
        Cn_bins = local_correlations_fft(Y, eight_neighbours=eight_neighbours, swap_dim=False)
        return Cn_bins
    else:
        if T % bin_size < bin_size / 2.:
            bin_size = T // (T // bin_size)

        n_bins = T // bin_size
        Cn_bins = np.zeros(((n_bins,) + Y.shape[1:]))
        for i in range(n_bins):
            Cn_bins[i] = local_correlations_fft(Y[i * bin_size:(i + 1) * bin_size],
                                                eight_neighbours=eight_neighbours,
                                                swap_dim=False)
            logging.debug(i * bin_size)

        Cn = np.max(Cn_bins, axis=0)
        return Cn


#%%
def local_correlations_fft(Y,
                           eight_neighbours: bool = True,
                           swap_dim: bool = True,
                           opencv: bool = True,
                           rolling_window=None) -> np.ndarray:
    """Computes the correlation image for the input dataset Y using a faster FFT based method

    Args:
        Y:  np.ndarray (3D or 4D)
            Input movie data in 3D or 4D format
    
        eight_neighbours: Boolean
            Use 8 neighbors if true, and 4 if false for 3D data (default = True)
            Use 6 neighbors for 4D data, irrespectively
    
        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front
    
        opencv: Boolean
            If True process using open cv method

        rolling_window: (undocumented)

    Returns:
        Cn: d1 x d2 [x d3] matrix, cross-correlation with adjacent pixels
    """

    if swap_dim:
        Y = np.transpose(Y, tuple(np.hstack((Y.ndim - 1, list(range(Y.ndim))[:-1]))))

    Y = Y.astype('float32')
    if rolling_window is None:
        Y -= np.mean(Y, axis=0)
        Ystd = np.std(Y, axis=0)
        Ystd[Ystd == 0] = np.inf
        Y /= Ystd
    else:
        Ysum = np.cumsum(Y, axis=0)
        Yrm = (Ysum[rolling_window:] - Ysum[:-rolling_window]) / rolling_window
        Y[:rolling_window] -= Yrm[0]
        Y[rolling_window:] -= Yrm
        del Yrm, Ysum
        Ystd = np.cumsum(Y**2, axis=0)
        Yrst = np.sqrt((Ystd[rolling_window:] - Ystd[:-rolling_window]) / rolling_window)
        Yrst[Yrst == 0] = np.inf
        Y[:rolling_window] /= Yrst[0]
        Y[rolling_window:] /= Yrst
        del Ystd, Yrst

    if Y.ndim == 4:
        if eight_neighbours:
            sz = np.ones((3, 3, 3), dtype='float32')
            sz[1, 1, 1] = 0
        else:
            # yapf: disable
            sz = np.array([[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                           [[0, 1, 0], [1, 0, 1], [0, 1, 0]],
                           [[0, 0, 0], [0, 1, 0], [0, 0, 0]]],
                          dtype='float32')
            # yapf: enable
    else:
        if eight_neighbours:
            sz = np.ones((3, 3), dtype='float32')
            sz[1, 1] = 0
        else:
            sz = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype='float32')

    if opencv and Y.ndim == 3:
        Yconv = np.stack([cv2.filter2D(img, -1, sz, borderType=0) for img in Y])
        MASK = cv2.filter2D(np.ones(Y.shape[1:], dtype='float32'), -1, sz, borderType=0)
    else:
        Yconv = convolve(Y, sz[np.newaxis, :], mode='constant')
        MASK = convolve(np.ones(Y.shape[1:], dtype='float32'), sz, mode='constant')

    YYconv = Yconv * Y
    del Y, Yconv
    if rolling_window is None:
        Cn = np.mean(YYconv, axis=0) / MASK
    else:
        YYconv_cs = np.cumsum(YYconv, axis=0)
        del YYconv
        YYconv_rm = (YYconv_cs[rolling_window:] - YYconv_cs[:-rolling_window]) / rolling_window
        del YYconv_cs
        Cn = YYconv_rm / MASK

    return Cn


def local_correlations_multicolor(Y, swap_dim: bool = True) -> np.ndarray:
    """Computes the correlation image with color depending on orientation

    Args:
        Y:  np.ndarray (3D or 4D)
            Input movie data in 3D or 4D format

        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front

    Returns:
        rho: d1 x d2 [x d3] matrix, cross-correlation with adjacent pixels
    """
    if Y.ndim == 4:
        raise Exception('Not Implemented')

    if swap_dim:
        Y = np.transpose(Y, tuple(np.hstack((Y.ndim - 1, list(range(Y.ndim))[:-1]))))

    w_mov = (Y - np.mean(Y, axis=0)) / np.std(Y, axis=0)

    rho_h = np.mean(np.multiply(w_mov[:, :-1, :], w_mov[:, 1:, :]), axis=0)
    rho_w = np.mean(np.multiply(w_mov[:, :, :-1], w_mov[:, :, 1:]), axis=0)
    rho_d1 = np.mean(np.multiply(w_mov[:, 1:, :-1], w_mov[:, :-1, 1:,]), axis=0)
    rho_d2 = np.mean(np.multiply(w_mov[:, :-1, :-1], w_mov[:, 1:, 1:,]), axis=0)

    return np.dstack([rho_h[:, 1:] / 2, rho_d1 / 2, rho_d2 / 2])


def local_correlations(Y, eight_neighbours: bool = True, swap_dim: bool = True, order_mean=1) -> np.ndarray:
    """Computes the correlation image for the input dataset Y

    Args:
        Y:  np.ndarray (3D or 4D)
            Input movie data in 3D or 4D format
    
        eight_neighbours: Boolean
            Use 8 neighbors if true, and 4 if false for 3D data (default = True)
            Use 6 neighbors for 4D data, irrespectively
    
        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front

        order_mean: (undocumented)

    Returns:
        rho: d1 x d2 [x d3] matrix, cross-correlation with adjacent pixels
    """

    if swap_dim:
        Y = np.transpose(Y, tuple(np.hstack((Y.ndim - 1, list(range(Y.ndim))[:-1]))))

    rho = np.zeros(np.shape(Y)[1:])
    w_mov = (Y - np.mean(Y, axis=0)) / np.std(Y, axis=0)

    rho_h = np.mean(np.multiply(w_mov[:, :-1, :], w_mov[:, 1:, :]), axis=0)
    rho_w = np.mean(np.multiply(w_mov[:, :, :-1], w_mov[:, :, 1:]), axis=0)

    # yapf: disable
    if order_mean == 0:
        rho = np.ones(np.shape(Y)[1:])
        rho_h = rho_h
        rho_w = rho_w
        rho[:-1, :] = rho[:-1, :] * rho_h
        rho[1:,  :] = rho[1:,  :] * rho_h
        rho[:, :-1] = rho[:, :-1] * rho_w
        rho[:,  1:] = rho[:,  1:] * rho_w
    else:
        rho[:-1, :] = rho[:-1, :] + rho_h**(order_mean)
        rho[1:,  :] = rho[1:,  :] + rho_h**(order_mean)
        rho[:, :-1] = rho[:, :-1] + rho_w**(order_mean)
        rho[:,  1:] = rho[:,  1:] + rho_w**(order_mean)

    if Y.ndim == 4:
        rho_d = np.mean(np.multiply(w_mov[:, :, :, :-1], w_mov[:, :, :, 1:]), axis=0)
        rho[:, :, :-1] = rho[:, :, :-1] + rho_d
        rho[:, :, 1:] = rho[:, :, 1:] + rho_d

        neighbors = 6 * np.ones(np.shape(Y)[1:])
        neighbors[0]        = neighbors[0]        - 1
        neighbors[-1]       = neighbors[-1]       - 1
        neighbors[:,     0] = neighbors[:,     0] - 1
        neighbors[:,    -1] = neighbors[:,    -1] - 1
        neighbors[:,  :, 0] = neighbors[:,  :, 0] - 1
        neighbors[:, :, -1] = neighbors[:, :, -1] - 1

    else:
        if eight_neighbours:
            rho_d1 = np.mean(np.multiply(w_mov[:, 1:, :-1], w_mov[:, :-1, 1:,]), axis=0)
            rho_d2 = np.mean(np.multiply(w_mov[:, :-1, :-1], w_mov[:, 1:, 1:,]), axis=0)

            if order_mean == 0:
                rho_d1 = rho_d1
                rho_d2 = rho_d2
                rho[:-1, :-1] = rho[:-1, :-1] * rho_d2
                rho[1:,   1:] = rho[1:,   1:] * rho_d1
                rho[1:,  :-1] = rho[1:,  :-1] * rho_d1
                rho[:-1,  1:] = rho[:-1,  1:] * rho_d2
            else:
                rho[:-1, :-1] = rho[:-1, :-1] + rho_d2**(order_mean)
                rho[1:,   1:] = rho[1:,   1:] + rho_d1**(order_mean)
                rho[1:,  :-1] = rho[1:,  :-1] + rho_d1**(order_mean)
                rho[:-1,  1:] = rho[:-1,  1:] + rho_d2**(order_mean)

            neighbors = 8 * np.ones(np.shape(Y)[1:3])
            neighbors[0,   :] = neighbors[0,   :] - 3
            neighbors[-1,  :] = neighbors[-1,  :] - 3
            neighbors[:,   0] = neighbors[:,   0] - 3
            neighbors[:,  -1] = neighbors[:,  -1] - 3
            neighbors[0,   0] = neighbors[0,   0] + 1
            neighbors[-1, -1] = neighbors[-1, -1] + 1
            neighbors[-1,  0] = neighbors[-1,  0] + 1
            neighbors[0,  -1] = neighbors[0,  -1] + 1
        else:
            neighbors = 4 * np.ones(np.shape(Y)[1:3])
            neighbors[0,  :]  = neighbors[0,  :] - 1
            neighbors[-1, :]  = neighbors[-1, :] - 1
            neighbors[:,  0]  = neighbors[:,  0] - 1
            neighbors[:, -1]  = neighbors[:, -1] - 1

    # yapf: enable
    if order_mean == 0:
        rho = np.power(rho, 1. / neighbors)
    else:
        rho = np.power(np.divide(rho, neighbors), 1 / order_mean)

    return rho


def correlation_pnr(Y, gSig=None, center_psf: bool = True, swap_dim: bool = True,
                    background_filter: str = 'disk') -> Tuple[np.ndarray, np.ndarray]:
    """
    compute the correlation image and the peak-to-noise ratio (PNR) image.
    If gSig is provided, then spatially filtered the video.

    Args:
        Y:  np.ndarray (3D or 4D).
            Input movie data in 3D or 4D format
        gSig:  scalar or vector.
            gaussian width. If gSig == None, no spatial filtering
        center_psf: Boolean
            True indicates subtracting the mean of the filtering kernel
        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front
        background_filter: str
            (undocumented)

    Returns:
        cn: np.ndarray (2D or 3D).
            local correlation image of the spatially filtered (or not)
            data
        pnr: np.ndarray (2D or 3D).
            peak-to-noise ratios of all pixels/voxels

    """
    if swap_dim:
        Y = np.transpose(Y, tuple(np.hstack((Y.ndim - 1, list(range(Y.ndim))[:-1]))))

    # parameters
    _, d1, d2 = Y.shape
    data_raw = Y.reshape(-1, d1, d2).astype('float32')

    # filter data
    data_filtered = data_raw.copy()
    if gSig:
        if not isinstance(gSig, list):
            gSig = [gSig, gSig]
        ksize = tuple([int(2 * i) * 2 + 1 for i in gSig])

        if center_psf:
            if background_filter == 'box':
                for idx, img in enumerate(data_filtered):
                    data_filtered[idx, ] = cv2.GaussianBlur(
                        img, ksize=ksize, sigmaX=gSig[0], sigmaY=gSig[1], borderType=1) \
                        - cv2.boxFilter(img, ddepth=-1, ksize=ksize, borderType=1)
            else:
                psf = cv2.getGaussianKernel(ksize[0], gSig[0],
                                            cv2.CV_32F).dot(cv2.getGaussianKernel(ksize[1], gSig[1], cv2.CV_32F).T)
                ind_nonzero = psf >= psf[0].max()
                psf -= psf[ind_nonzero].mean()
                psf[~ind_nonzero] = 0
                for idx, img in enumerate(data_filtered):
                    data_filtered[idx,] = cv2.filter2D(img, -1, psf, borderType=1)

            # data_filtered[idx, ] = cv2.filter2D(img, -1, psf, borderType=1)
        else:
            for idx, img in enumerate(data_filtered):
                data_filtered[idx,] = cv2.GaussianBlur(img, ksize=ksize, sigmaX=gSig[0], sigmaY=gSig[1], borderType=1)

    # compute peak-to-noise ratio
    data_filtered -= data_filtered.mean(axis=0)
    data_max = np.max(data_filtered, axis=0)
    data_std = get_noise_fft(data_filtered.T, noise_method='mean')[0].T
    pnr = np.divide(data_max, data_std)
    pnr[pnr < 0] = 0

    # remove small values
    tmp_data = data_filtered.copy() / data_std
    tmp_data[tmp_data < 3] = 0

    # compute correlation image
    cn = local_correlations_fft(tmp_data, swap_dim=False)

    return cn, pnr


def iter_chunk_array(arr: np.array, chunk_size: int):
    if ((arr.shape[0] // chunk_size) - 1) > 0:
        for i in range((arr.shape[0] // chunk_size) - 1):
            yield arr[chunk_size * i:chunk_size * (i + 1)]
        yield arr[chunk_size * (i + 1):]
    else:
        yield arr





def prepare_local_correlations(Y, swap_dim: bool = False,
                               eight_neighbours: bool = False) -> Tuple[Any, Any, Any, Any, Any, Any, Any, Any]:
    """Computes the correlation image and some statistics to update it online

    Args:
        Y:  np.ndarray (3D or 4D)
            Input movie data in 3D or 4D format

        swap_dim: Boolean
            True indicates that time is listed in the last axis of Y (matlab format)
            and moves it in the front

        eight_neighbours: Boolean
            Use 8 neighbors if true, and 4 if false for 3D data
            Use 18 neighbors if true, and 6 if false for 4D data

    """
    # TODO: Tighten prototype above
    if swap_dim:
        Y = np.transpose(Y, (Y.ndim - 1,) + tuple(range(Y.ndim - 1)))

    T = len(Y)
    dims = Y.shape[1:]
    Yr = Y.T.reshape(-1, T)
    if Y.ndim == 4:
        d1, d2, d3 = dims
        sz = generate_binary_structure(3, 2 if eight_neighbours else 1)
        sz[1, 1, 1] = 0
    else:
        d1, d2 = dims
        if eight_neighbours:
            sz = np.ones((3, 3), dtype='uint8')
            sz[1, 1] = 0
        else:
            sz = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype='uint8')

    idx = [i - 1 for i in np.nonzero(sz)]

    def get_indices_of_neighbors(pixel):
        pixel = np.unravel_index(pixel, dims, order='F')
        x = pixel[0] + idx[0]
        y = pixel[1] + idx[1]
        if len(dims) == 3:
            z = pixel[2] + idx[2]
            inside = (x >= 0) * (x < d1) * (y >= 0) * (y < d2) * (z >= 0) * (z < d3)
            return np.ravel_multi_index((x[inside], y[inside], z[inside]), dims, order='F')
        else:
            inside = (x >= 0) * (x < d1) * (y >= 0) * (y < d2)
            return np.ravel_multi_index((x[inside], y[inside]), dims, order='F')

    # more compact but slower code
    # idx = np.asarray([i - 1 for i in np.nonzero(sz)])
    # def get_indices_of_neighbors(pixel):
    #     pixel = np.asarray(np.unravel_index(pixel, dims, order='F'))
    #     xyz = pixel[:, None] + idx
    #     inside = np.all([(x >= 0) * (x < d) for (x, d) in zip(xyz, dims)], 0)
    #     return np.ravel_multi_index(xyz[:, inside], dims, order='F')

    N = [get_indices_of_neighbors(p) for p in range(np.prod(dims))]
    col_ind = np.concatenate(N)
    row_ind = np.concatenate([[i] * len(k) for i, k in enumerate(N)])
    num_neigbors = np.concatenate([[len(k)] * len(k) for k in N]).astype(Yr.dtype)

    first_moment = Yr.mean(1)
    second_moment = (Yr**2).mean(1)
    crosscorr = np.mean(Yr[row_ind] * Yr[col_ind], 1)
    # slower for small T, less memory intensive, but memory not an issue:
    # crosscorr = np.array([Yr[r_].dot(Yr[c_])
    #                       for (r_, c_) in zip(row_ind, col_ind)]) / Yr.shape[1]
    sig = np.sqrt(second_moment - first_moment**2)

    M = coo_matrix(
        ((crosscorr - first_moment[row_ind] * first_moment[col_ind]) / (sig[row_ind] * sig[col_ind]) / num_neigbors,
         (row_ind, col_ind)),
        dtype=Yr.dtype)
    cn = M.dot(np.ones(M.shape[1], dtype=M.dtype)).reshape(dims, order='F')

    return first_moment, second_moment, crosscorr, col_ind, row_ind, num_neigbors, M, cn


def update_local_correlations(t,
                              frames,
                              first_moment,
                              second_moment,
                              crosscorr,
                              col_ind,
                              row_ind,
                              num_neigbors,
                              M,
                              del_frames=None) -> np.ndarray:
    """Updates sufficient statistics in place and returns correlation image"""
    dims = frames.shape[1:]
    stride = len(frames)
    if stride:
        frames = frames.reshape((stride, -1), order='F')
        if del_frames is None:
            tmp = 1 - float(stride) / t
            first_moment *= tmp
            second_moment *= tmp
            crosscorr *= tmp
        else:
            if stride > 10:
                del_frames = del_frames.reshape((stride, -1), order='F')
                first_moment -= del_frames.sum(0) / t
                second_moment -= (del_frames**2).sum(0) / t
                crosscorr -= np.sum(del_frames[:, row_ind] * del_frames[:, col_ind], 0) / t
            else:      # loop is faster
                for f in del_frames:
                    f = f.ravel(order='F')
                    first_moment -= f / t
                    second_moment -= (f**2) / t
                    crosscorr -= (f[row_ind] * f[col_ind]) / t
        if stride > 10:
            frames = frames.reshape((stride, -1), order='F')
            first_moment += frames.sum(0) / t
            second_moment += (frames**2).sum(0) / t
            crosscorr += np.sum(frames[:, row_ind] * frames[:, col_ind], 0) / t
        else:          # loop is faster
            for f in frames:
                f = f.ravel(order='F')
                first_moment += f / t
                second_moment += (f**2) / t
                crosscorr += (f[row_ind] * f[col_ind]) / t


#=======
#            del_frames = del_frames.reshape((stride, -1), order='F')
#            first_moment -= del_frames.sum(0) / t
#            second_moment -= (del_frames**2).sum(0) / t
#            crosscorr -= np.sum(del_frames[:, row_ind] * del_frames[:, col_ind], 0) / t
#        else:                                                                                               # loop is faster
#            for f in del_frames:
#                f = f.ravel(order='F')
#                first_moment -= f / t
#                second_moment -= (f**2) / t
#                crosscorr -= (f[row_ind] * f[col_ind]) / t
#    if stride > 10:
#        frames = frames.reshape((stride, -1), order='F')
#        first_moment += frames.sum(0) / t
#        second_moment += (frames**2).sum(0) / t
#        crosscorr += np.sum(frames[:, row_ind] * frames[:, col_ind], 0) / t
#    else:                                                                                                   # loop is faster
#        for f in frames:
#            f = f.ravel(order='F')
#            first_moment += f / t
#            second_moment += (f**2) / t
#            crosscorr += (f[row_ind] * f[col_ind]) / t
#>>>>>>> dev
    sig = np.sqrt(second_moment - first_moment**2)
    M.data = ((crosscorr - first_moment[row_ind] * first_moment[col_ind]) / (sig[row_ind] * sig[col_ind]) /
              num_neigbors)
    cn = M.dot(np.ones(M.shape[1], dtype=M.dtype)).reshape(dims, order='F')
    return cn




   