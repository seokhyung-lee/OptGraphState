import logging
import time
import numpy as np
import scipy.signal as scsig
import scipy.special as scsp
import decimal as dec
import math
import sympy as sp


def get_ftpdf_contracted(ftpdf1, ftpdf2, p_succ=0.5):
    """
    Compute the fourier-transformed probability distribution functions (FTPDF)
    of the resource count of a new node created by contracting a link of a
    fusion network from the FTPDFs of the two endpoints of the link.

    High-precision calculation with the `decimal` package is supported.

    Parameters
    ----------
    ftpdf1 : 1D array-like of float or decimal.Decimal
        FTPDF (with variable k) of an endpoint of the contracted link, which is
        the inverse of the polynomial of x=exp(2πik) where the coefficient for
        x^j is `ftpdf1[j]`.
    ftpdf2 : 1D array-like of float or decimal.Decimal
        FTPDF of the other endpoint expressed in the same way as above.
    p_succ : float (default: 0.5)
        Success probability of a single fusion.

    Returns
    -------
    ftpdf_cont : 1D numpy array of float or decimal.Decimal
        FTPDF of the new node created by contracting the link, expressed in
        the same way as `ftpdf1` and `ftpdf2`.
    """
    ftpdf1 = np.asanyarray(ftpdf1)
    ftpdf2 = np.asanyarray(ftpdf2)

    if isinstance(ftpdf1[0], dec.Decimal):
        p_succ = dec.Decimal(str(p_succ))

    ftpdf_cont = scsig.convolve(ftpdf1, ftpdf2) / p_succ
    ftpdf_cont[0] -= (1 - p_succ) / p_succ

    return ftpdf_cont


def recover_cmf_from_ftpdf(ftpdf,
                           cutoff,
                           raise_error_when_invalid_prob=True):
    """
    Recover the cumulative mass function (CMF) of the resource count for a node
    in a fusion network from its fourier-transformed probability distribution
    functions (FTPDF).

    High-precision calculation with the `decimal` package is supported.

    Parameters
    ----------
    ftpdf : 1D array-like of float or decimal.Decimal
        FTPDF (with variable k) of a node of a fusion network, which is
        the inverse of the polynomial of x=exp(2πik) where the coefficient for
        x^j is `ftpdf[j]`.
    cutoff : float (default: 0.95)
        Cutoff threshold between 0 and 1.
        The CMF is computed until it exceeds the given threshold.
    raise_error_when_invalid_prob : bool (default: True)
        Whether to raise an error when the calculated probability is invalid.

    Returns
    -------
    L : int
        First resource count that gives a nonzero probability.
    cmfs : 1D numpy array of float or decimal.Decimal
        CMF values for the resource counts starting from `m_start`.
        `cmfs[i]` is the value corresponding to the resource count of `m_start+i`.
    """
    assert 0 < cutoff < 1

    ftpdf = np.asanyarray(ftpdf)
    decimal = isinstance(ftpdf[0], dec.Decimal)

    L = ftpdf.size - 1
    rescale_factor = ftpdf[-1]

    ftpdf_diff = np.insert(np.diff(ftpdf), 0, ftpdf[0])
    ftpdf_diff_res = ftpdf_diff / rescale_factor

    one = dec.Decimal('1') if decimal else 1.
    cmfs = np.array([one, one])

    cutoff = dec.Decimal(str(cutoff)) if decimal else cutoff
    cutoff_rescaled = cutoff * rescale_factor

    summ = np.sum if decimal else math.fsum
    n = 2
    prev_cmf = 0
    while True:
        cmf = summ(ftpdf_diff_res[max(L - n + 1, 0):]
                   * cmfs[max(n - L - 1, 0):])

        if raise_error_when_invalid_prob \
                and (cmf < 0 or cmf > rescale_factor or cmf < prev_cmf):
            raise ValueError("Invalid probability. Use a higher precision.")

        cmfs = np.append(cmfs, cmf)
        if cmfs[-1] >= cutoff_rescaled:
            break

        n += 1
        prev_cmf = cmf

    cmfs /= rescale_factor

    return L, cmfs


def set_log(fname):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    file_handler = logging.FileHandler(fname)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def logging_time(original_fn):
    def wrapper_fn(*args, **kwargs):
        start_time = time.time()
        result = original_fn(*args, **kwargs)
        end_time = time.time()
        print("WorkingTime[{}]: {} sec".format(original_fn.__name__,
                                               end_time - start_time))
        return result

    return wrapper_fn
