import numpy as np
import astropy.units as u
import warnings
import scipy.optimize
import scipy.interpolate
import pkg_resources
import os
import pickle
import pandas


def calc_s_dMag(a, e, I, w, R_P, p, nu):  # noqa: E741
    r"""Compute the projected separation and \Delta{mag}

    Args:
        a (arraylike or Quantity):
            semi-major axis
        e (arraylike):
            eccentricity
        I (arraylike or Quantity):
            Inclination. If not a quantity, must be in radians
        w (arraylike or Quantity):
            Argument of periapse. If not a quantity, must be in radians
        R_P (arraylike or Quantity):
            Planet radius.  If not quantity, must be in same units as
            semi-major axis
        p (arraylike):
            Geometric albedo
        nu (arraylike or Quantity):
            True anomaly. If not a quantity, must be in radians

    Returns:
        tuple:
            s (arraylike or quantity):
                projected separation.  If quantity, in same units as
                semi-major axis input
            dMag (arraylike):
                \Delta{mag} values

    Notes:
        All inputs must be of the same size, or scalars.

    """

    r = a * (1 - e**2) / (1 + e * np.cos(nu))
    th = w + nu
    s = r * np.sqrt(1 - np.sin(I) ** 2 * np.sin(th) ** 2)
    beta = np.arccos(np.sin(I) * np.sin(th))
    Phi = Lambert_phase_function(beta)
    dMag = -2.5 * np.log10(p * Phi * (R_P / r) ** 2)

    return s, dMag


def Lambert_phase_function(beta):
    """Evaluate the Lambert phase function

    Args:
        beta (arraylike or Quantity):
            Phase angle.  If not quantity, must be in radians

    Returns:
        arraylike:
            Phase function values.
    """

    if isinstance(beta, u.Quantity):
        Phi = (np.sin(beta) + (np.pi - beta.to(u.rad).value) * np.cos(beta)) / np.pi
    else:
        Phi = (np.sin(beta) + (np.pi - beta) * np.cos(beta)) / np.pi

    return Phi


def forcendarray(x):
    """Convert any numerical value into 1-D ndarray

    Args:
        x (float or numpy.ndarray):
            Input
    Returns:
        numpy.ndarray:
            Same size as input but in ndarray form
    """

    return np.array(x, ndmin=1).astype(float).flatten()


def invKepler(
    M,
    e,
    tol=None,
    E0=None,
    maxIter=100,
    return_nu=False,
    convergence_error=True,
):
    """Finds eccentric/hyperbolic/parabolic anomaly from mean anomaly and eccentricity

    This method uses Newton-Raphson iteration to find the eccentric
    anomaly from mean anomaly and eccentricity, assuming a closed (0<e<1)
    orbit.

    Args:
        M (float or ndarray):
            mean anomaly (rad)
        e (float or ndarray):
            eccentricity (eccentricity may be a scalar if M is given as
            an array, but otherwise must match the size of M.)
        tolerance (float):
            Convergence of tolerance. Defaults to eps(2*pi)
        E0 (float or ndarray):
            Initial guess for iteration.  Defaults to Taylor-expansion based value for
            closed orbits and Vallado-derived heuristic for open orbits. If set, must
            match size of M.
        maxiter (int):
            Maximum numbr of iterations.  Optional, defaults to 100.
        return_nu (bool):
            Return true anomaly (defaults false)
        convergence_error (bool):
            Raise error on convergence failure. Defaults True.  If false, throws a
            warning.

    Returns:
        tuple:
            E (ndarray):
                eccentric/parabolic/hyperbolic anomaly (rad)
            numIter (ndarray):
                Number of iterations
            nu (ndarray):
                True anomaly (returned only if return_nu=True)

    Notes:
        If either M or e are scalar, and the other input is an array, the scalar input
        will be expanded to the same size array as the other input.  So, a scalar M
        and array e will result in the calculation of the eccentric anomaly for one
        mean anomaly at a variety of eccentricities, and a scalar e and array M input
        will result in the calculation of eccentric anomalies for one eccentricity at
        a variety of mean anomalies.  If both inputs are arrays then they are matched
        element by element.

    """

    # make sure M and e are of the correct format.
    # if either is scalar, expand to match sizes
    M = forcendarray(M)
    e = forcendarray(e)
    if e.size != M.size:
        if e.size == 1:
            e = np.array([e[0]] * len(M))
        if M.size == 1:
            M = np.array([M[0]] * len(e))

    assert e.shape == M.shape, "Incompatible inputs."
    assert np.all((e >= 0)), "e values below zero"

    if E0 is not None:
        E0 = forcendarray(E0)
        assert E0.shape == M.shape, "Incompatible inputs."

    # define output
    Eout = np.zeros(M.size)
    numIter = np.zeros(M.size, dtype=int)

    # circles
    cinds = e == 0
    Eout[cinds] = M[cinds]

    # ellipses
    einds = (e > 0) & (e < 1)
    if any(einds):

        Me = np.mod(M[einds], 2 * np.pi)
        ee = e[einds]

        # initial values for E
        if E0 is None:
            E = Me / (1 - ee)
            mask = ee * E**2 > 6 * (1 - ee)
            E[mask] = np.cbrt(6 * Me[mask] / ee[mask])
        else:
            E = E0[einds]

        # Newton-Raphson setup
        counter = np.ones(E.shape)
        err = np.ones(E.shape)

        # set tolerance is none provided
        if tol is None:
            etol = np.spacing(2 * np.pi)
        else:
            etol = tol

        while (np.max(err) > etol) and (np.max(counter) < maxIter):
            inds = err > etol
            E[inds] = E[inds] - (Me[inds] - E[inds] + ee[inds] * np.sin(E[inds])) / (
                ee[inds] * np.cos(E[inds]) - 1
            )
            err[inds] = np.abs(Me[inds] - (E[inds] - ee[inds] * np.sin(E[inds])))
            counter[inds] += 1

        if np.max(counter) == maxIter:
            if convergence_error:
                raise ValueError("Maximum number of iterations exceeded")
            else:
                warnings.warn("Maximum number of iterations exceeded")

        Eout[einds] = E
        numIter[einds] = counter

    # parabolae
    pinds = e == 1
    if np.any(pinds):
        q = 9 * M[pinds] / 6
        B = (q + np.sqrt(q**2 + 1)) ** (1.0 / 3.0) - (np.sqrt(q**2 + 1) - q) ** (
            1.0 / 3.0
        )
        Eout[pinds] = B

    # hyperbolae
    hinds = e > 1
    if np.any(hinds):
        Mh = M[hinds]
        eh = e[hinds]

        # initialize H
        if E0 is None:
            H = Mh / (eh - 1)
            mask = eh * H**2 > 6 * (eh - 1)
            H[mask] = np.cbrt(6 * Mh[mask] / eh[mask])
        else:
            H = E0[hinds]

        # Newton-Raphson setup
        counter = np.ones(H.shape)

        # set tolerance is none provided
        if tol is None:
            htol = 4 * np.spacing(np.abs(H))
        else:
            htol = np.ones(H.shape) * tol

        Hup = np.ones(len(H))
        # we will only update things that haven't hit their tolerance:
        inds = np.abs(Hup) > htol

        while np.any(inds) and np.all(counter < maxIter):
            Hup[inds] = (Mh[inds] - eh[inds] * np.sinh(H[inds]) + H[inds]) / (
                eh[inds] * np.cosh(H[inds]) - 1
            )
            H[inds] = H[inds] + Hup[inds]
            if tol is None:
                htol[inds] = 4 * np.spacing(np.abs(H[inds]))

            counter[inds] += 1
            inds = np.abs(Hup) > htol

        if np.max(counter) == maxIter:
            if convergence_error:
                raise ValueError("Maximum number of iterations exceeded")
            else:
                warnings.warn("Maximum number of iterations exceeded")

        Eout[hinds] = H
        numIter[hinds] = counter

    out = (Eout, numIter)
    if return_nu:
        nuout = np.zeros(M.size)

        # circles
        nuout[cinds] = M[cinds]

        # ellipses
        if np.any(einds):
            nuout[einds] = np.mod(
                2 * np.arctan(np.sqrt((1 + e[einds]) / (1 - e[einds])) * np.tan(E / 2)),
                2 * np.pi,
            )

        # parabolae
        if np.any(pinds):
            nuout[pinds] = 2 * np.arctan(B)

        # hyperbolae
        if np.any(hinds):
            nuout[hinds] = 2 * np.arctan(
                np.sqrt((e[hinds] + 1) / (e[hinds] - 1)) * np.tanh(H / 2)
            )

        out += (nuout,)

    return out


def observable_indices(s, dMag, d, IWA, OWA, dMaglim):
    """Compute when a planet is observable

    Args:
        s (arraylike or quantity):
            Projected separation.  In AU if no quantity.
        dMag (arraylike):
            Delta mag values
        d (float):
            Distance to star in parsecs
        IWA (float):
            Inner working angle in arcseconds
        OWA (float):
            Outer working angle in arcseconds
        dMaglim (float):
            Limiting Delta mag

    Returns:
        array:
            Boolean array that is true for indices of s and dMag where the planet is
            observable and false otherwise
    """

    if isinstance(s, u.Quantity):
        stmp = s.to(u.AU).value
    else:
        stmp = s

    return (stmp <= d * OWA) & (stmp >= d * IWA) & (dMag < dMaglim)


def split_observable_inds(observable_indices):
    """Split the observable indices into contiguous blocks"""

    inds = np.where(observable_indices)[0]
    j = 0
    sepinds = []
    for bind in np.where(np.diff(inds) > 1)[0]:
        sepinds.append(inds[j : bind + 1])  # noqa: E203
        j = bind + 1
    sepinds.append(inds[j:])

    return sepinds


def gen_Earthlike_values(N):
    """Generate a sample of vaugely Earth-like parameters

    Args:
        N (int):
            Number of samples to generate

    Returns:
        tuple:
            a (np.ndarray):
                semi-major axis values (AU)
            e (np.ndarray):
                eccentricity values
            w (np.ndarray):
                arg. of periapsis values (rad)
            I (np.ndarray):
                inclination values (rad)
            M (np.ndarray):
                Mean anomaly values (rad)
            R_P (float):
                1 Earth radius in AU
            p (float):
                Geometric albedo = 0.367
    """

    # generate an array of semi-major axis values
    # uniformly distributed between 0.7 and 1.5 AU:
    a = np.random.uniform(low=0.7, high=1.5, size=N)
    # generate an array of eccentricity values
    # uniformly distributed between 0 and 0.35
    e = np.random.uniform(low=0, high=0.35, size=N)
    # generate an array of argument of periapsis values
    # uniformly distributed between 0 and 2\pi
    w = np.random.uniform(low=0, high=2 * np.pi, size=N)
    # generate an array of inclination values
    # sinusoidally distributed between 0 and \pi
    I = np.arccos(np.random.uniform(low=-1, high=1, size=N))  # noqa: E741
    # Generate an array of mean anomaly values
    # uinformly distributed between 0 and 2\pi
    M = np.random.uniform(low=0, high=2 * np.pi, size=N)

    # these values are constant
    R_P = (1 * u.earthRad).to(u.AU).value  # Earth radius-equivalent planet in AU
    p = 0.367  # geometric albedo

    return a, e, w, I, M, R_P, p


def gen_Cpdf(bins, dMag, s, arange, erange, p, R_P):
    """Generate the 2D historgram of the joint imaging observable PDFs

    Args:
        bins (int):
            Number of bins to use in each dimension
        dMag (np.ndarray):
            Array of Delta mag values
        s (np.ndarray):
            Array of projected separation values (AU)
        arange (list):
            [min sma, max sma] in AU
        erange (list):
            [min e, max e]
        p (float):
            Geometric albedo
        R_P (float):
            Planet radius in AU

    Returns:
        tuple:
            Cpdf (np.ndarray):
                2D array representing the joint distribution
            sax (np.ndarray):
                1D array of values along the projected separation axis
            dMagax (np.ndarray):
                1D array of values along the dMag axis

    """

    N = len(dMag)  # number of samples

    # compute total orbital radius range
    rrange = [arange[0] * (1 - erange[1]), arange[1] * (1 + erange[1])]

    # define bins
    xedges = np.linspace(0.0, rrange[1], bins + 1)
    ymin = -2.5 * np.log10(p * (R_P / rrange[0]) ** 2)
    ymax = -2.5 * np.log10(p * (R_P / rrange[1]) ** 2 * 1e-11)
    yedges = np.linspace(ymin, ymax, bins + 1)

    # generate 2D histogram and normalize by bins
    Cpdf, yedges, xedges = np.histogram2d(
        dMag,
        s,
        bins=bins,
        range=[[yedges.min(), yedges.max()], [xedges.min(), xedges.max()]],
    )
    Cpdf /= N * (xedges[1] - xedges[0]) * (yedges[1] - yedges[0])

    # comppute bin centers
    xcent = 0.5 * (xedges[1:] + xedges[:-1])
    ycent = 0.5 * (yedges[1:] + yedges[:-1])

    # define final outputs
    sax = np.hstack((0.0, xcent, rrange[1]))
    dMagax = np.hstack((ymin, ycent, ymax))
    Cpdf = np.pad(Cpdf, 1, mode="constant")

    return Cpdf, sax, dMagax


def calc_Cpdf_limits(sax, arange, erange, p, R_P):
    """Compute limits on the joint imaging observable PDFs

    Args:
        sax (np.ndarray):
            Output of genCpdf
        arange (list):
            [min sma, max sma] in AU
        erange (list):
            [min e, max e]
        p (float):
            Geometric albedo
        R_P (float):
            Planet radius in AU

    Returns:
        tuple:
            dmagmax (np.ndarray):
                Maximum dmag values as function of sax
            dmagmin (np.ndarray):
                Minimum dmag values as function of sax
            dmag90 (np.ndarray):
                dmag values at quadrature as function of sax
    """
    # compute total orbital radius range
    rrange = [arange[0] * (1 - erange[1]), arange[1] * (1 + erange[1])]

    # compute dmag limits
    dmagmax = -2.5 * np.log10(
        p
        * ((R_P / rrange[1])) ** 2
        * Lambert_phase_function(np.pi - np.arcsin(sax / rrange[1]))
    )

    betastarexpr = (
        lambda beta: -(np.pi - beta) * np.sin(beta) ** 3 / np.pi
        + 2
        * ((np.pi - beta) * np.cos(beta) + np.sin(beta))
        * np.sin(beta)
        * np.cos(beta)
        / np.pi
    )
    betastar = scipy.optimize.fsolve(betastarexpr, 63 * np.pi / 180)[0]

    bp1 = rrange[0] * np.sin(betastar)
    bp2 = rrange[1] * np.sin(betastar)
    dmagmin = np.zeros(sax.size)
    dmagmin[sax < bp1] = -2.5 * np.log10(
        p
        * (R_P / rrange[0]) ** 2
        * Lambert_phase_function(np.arcsin(sax[sax < bp1] / rrange[0]))
    )
    inds = (sax >= bp1) & (sax < bp2)
    dmagmin[inds] = -2.5 * np.log10(
        p
        * (R_P / sax[inds]) ** 2
        * Lambert_phase_function(betastar)
        * np.sin(betastar) ** 2
    )
    dmagmin[sax >= bp2] = -2.5 * np.log10(
        p
        * (R_P / rrange[1]) ** 2
        * Lambert_phase_function((np.arcsin(sax[sax >= bp2] / rrange[1])))
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dmag90 = -2.5 * np.log10(
            p
            * (R_P / sax) ** 2
            * Lambert_phase_function(np.pi / 2)
            * np.sin(np.pi / 2) ** 2
        )

    return dmagmax, dmagmin, dmag90


def load_precomputed_completeness():
    """Load pre-computed completeness distribution

    Args:
        None

    Returns:
        tuple:
            Cpdf (np.ndarray):
                2D array representing the joint distribution
            sax (np.ndarray):
                1D array of values along the projected separation axis
            dMagax (np.ndarray):
                1D array of values along the dMag axis

    """

    fname = pkg_resources.resource_filename(
        "SSWYieldModelingTutorial", os.path.join("data", "precomputed_Cpdf.pkl")
    )

    with open(fname, "rb") as ff:
        Cpdf = pickle.load(ff)

    # compute total orbital radius range
    arange = [0.7, 1.5]
    erange = [0, 0.35]
    rrange = [arange[0] * (1 - erange[1]), arange[1] * (1 + erange[1])]
    R_P = (1 * u.earthRad).to(u.AU).value  # Earth radius-equivalent planet in AU
    p = 0.367  # geometric albedo
    bins = 1000

    # define bins
    xedges = np.linspace(0.0, rrange[1], bins + 1)
    ymin = -2.5 * np.log10(p * (R_P / rrange[0]) ** 2)
    ymax = -2.5 * np.log10(p * (R_P / rrange[1]) ** 2 * 1e-11)
    yedges = np.linspace(ymin, ymax, bins + 1)

    # comppute bin centers
    xcent = 0.5 * (xedges[1:] + xedges[:-1])
    ycent = 0.5 * (yedges[1:] + yedges[:-1])

    # define final outputs
    sax = np.hstack((0.0, xcent, rrange[1]))
    dMagax = np.hstack((ymin, ycent, ymax))
    Cpdf = np.pad(Cpdf, 1, mode="constant")

    return Cpdf, sax, dMagax


def calc_completeness(Cpdf, sax, dMagax, smin, smax, dMaglim, L=1):
    """Compute the completeness of an observation

    Args:
        Cpdf (np.ndarray):
            2D, normalized, joint probability histogram of s and dMag
        sax (np.ndarray):
            Projected separation axis of Cpdf (AU)
        dMagax (np.ndarray):
            dMag axis of Cpdf
        smin (arraylike):
            Minimum observable projected separation (projected IWA) (AU)
        smax (arraylike):
            Maximum observable projected separation (projected OWA) (AU)
        dMaglim (arraylike):
            Maximum observable Delta mag value
        L (arraylike):
            Stellar luminosity in solar luminosities. Defaults to 1.

    Returns:
        np.ndarray:
            Completeness values

    Notes:
        All arraylike inputs must have the same dimensionalities (or be scalars)

    """

    # define interpolant and summation function
    EVPOCpdf = scipy.interpolate.RectBivariateSpline(sax, dMagax, Cpdf.T)
    EVPOC = np.vectorize(EVPOCpdf.integral, otypes=[np.float64])

    # scale input by luminosity
    scaled_smin = np.array(smin, ndmin=1) / np.sqrt(L)
    scaled_smax = np.array(smax, ndmin=1) / np.sqrt(L)
    scaled_dMag = np.array(dMaglim, ndmin=1) - 2.5 * np.log10(L)

    # ensure that everything is the same length
    maxlen = np.max([scaled_smin.size, scaled_smax.size, scaled_dMag.size])
    if scaled_smin.size < maxlen:
        scaled_smin = np.repeat(scaled_smin, maxlen)
    if scaled_smax.size < maxlen:
        scaled_smax = np.repeat(scaled_smax, maxlen)
    if scaled_dMag.size < maxlen:
        scaled_dMag = np.repeat(scaled_dMag, maxlen)

    # we need to mask out values outside the bounds of the interpolant
    mask = (scaled_dMag > dMagax.min()) & (smin < sax.max())

    # compute completeness
    comp = np.zeros(maxlen)
    comp[mask] = EVPOC(scaled_smin[mask], scaled_smax[mask], 0.0, scaled_dMag[mask])
    # remove small values
    comp[comp < 1e-6] = 0.0

    return comp


def calc_intTime(C_p, C_b, M, SNR):
    """Find the integration time to reach a required SNR given the planet and
    background count rates as well as the optical system's noise floor.


    Args:
        C_p (arraylike Quantity):
            Planet count rate (1/time units)
        C_b (arraylike Quantity):
            Background count rate (1/time units)
        M (arraylike Quantity):
            Noise floor count rate (1/time units)
        SNR (float):
            Required signal to noise ratio

    Returns:
        ~astropy.units.Quantity(~numpy.ndarray(float)):
            Integration times

    .. note::

        All infeasible integration times should be returned as NaN values

    """

    intTime = (C_p + C_b) / ((C_p / SNR) ** 2 - M**2)

    # infinite and negative values are set to NAN
    intTime[np.isinf(intTime) | (intTime.value < 0.0)] = np.nan

    return intTime


def Cp_Cb_M(static_params, coronagraph, target, deltaMag):
    """Calculates electron count rates for planet signal, background noise,
    and noise floor for an observation.

    Args:
        static_params (dict):
            Dictionary of static parameters:
            lam (Quantity):
                Central wavelength of observing bandpass (length unit)
            deltaLam (Quantity):
                Bandpass (length unit)
            D (Quantity):
                Telescope aperture diameter (length unit)
            obsc (float):
                Fraction of the primary aperture that is obscured by secondary and
                secondary support structures
            tau (float):
                Optical system throughput excluding effects of starlight suppression
                system
            QE (float):
                Detector quantum efficiency
            F0 (Quantity):
                Spectral flux density of zero-magnitude star in observing band
                (1/length^2/length/time unit)
            pixelScale (Quantity):
                Instantaneous field of view of each detector pixel (angle unit)
            darkCurrent (Quantity):
                Dark current in counts/second/pixel (1/time unit)
            readNoise (float):
                Read noise in electrons/pixel/read
            texp (Quantity):
                Single readout exposure time (time unit)
            ppFac (float):
                Post-processing factor

        cornagraph (dict):
            Dictionary of coronagraph parameters for this observation:
            tau_core (float):
                Throughput of starlight suppression system for point sources
            tau_occ (float):
                Throughput of starlight suppression system for infinitely extended
                sources

        target (dict):
            Dictionary of target observing parameters:
            mag_star (float):
                Apparennt magnitude of target star in observing band
            zodi (Quantity)
                Surface brightness of local zodi in units of magnitude/arcsec^2.
                Input must have units of angle^{-2}.
            exozodi (Quantity)
                Surface brightness of exozodi in units of magnitude/arcsec^2.
                Input must have units of angle^{-2}.
        deltaMag (float):
                Assumed delta magnitude of planet in observing band


    Returns:
        tuple:
            C_p (~astropy.units.Quantity(~numpy.ndarray(float))):
                Planet signal electron count rate in units of 1/s
            C_b (~astropy.units.Quantity(~numpy.ndarray(float))):
                Background noise electron count rate in units of 1/s
            M (~astropy.units.Quantity(~numpy.ndarray(float))):
                Residual starlight spatial structure (systematic error)
                in units of 1/s
    """

    # Compute telescope collecting area:
    # this value should have units of length^2
    A = np.pi * (static_params["D"] / 2) ** 2 * (1 - static_params["obsc"])

    # compute the common factor of A*tau*deltaLam*QE
    # this value should have units of length^3
    eta = A * static_params["tau"] * static_params["deltaLam"] * static_params["QE"]

    # Compute the count rates of the planet
    # this should have units of 1/time
    C_star = static_params["F0"] * 10 ** (-0.4 * target["mag_star"])
    C_p = (C_star * 10 ** (-0.4 * deltaMag) * eta * coronagraph["tau_core"]).decompose()

    # Compute size of critically sampled photometric aperture
    # This should have units of angle^2
    Omega = np.pi * ((static_params["lam"] / 2 / static_params["D"]) ** 2).to(
        u.arcsec**2, equivalencies=u.dimensionless_angles()
    )

    # Compute the local and exozodi count rates
    # this should have units of 1/time
    C_zodi = (
        static_params["F0"]
        * 10 ** (-0.4 * target["zodi"])
        * Omega
        / u.arcsec**2
        * eta
        * coronagraph["tau_occ"]
    ).decompose()
    # Compute the local and exozodi count rates
    C_exozodi = (
        static_params["F0"]
        * 10 ** (-0.4 * target["exozodi"])
        * Omega
        / u.arcsec**2
        * eta
        * coronagraph["tau_occ"]
    ).decompose()

    # Compute the count rates of the starlight residual
    # this should have units of 1/time
    C_sr = (
        C_star * eta * coronagraph["contrast"] * coronagraph["tau_core"]
    ).decompose()

    # number of detector pixels in the photometric aperture = Omega / theta^2
    # this value should be unitless
    Npix = (Omega / static_params["pixelScale"] ** 2.0).decompose().value

    # Compute the dark current count rate
    # this should have units of 1/time
    C_dc = Npix * static_params["darkCurrent"]

    # Compute the read noise count rate
    # this should have units of 1/time
    C_rn = Npix * static_params["readNoise"] / static_params["texp"]

    # total background signal rate
    # this should have units of 1/time
    C_b = C_sr + C_zodi + C_exozodi + C_dc + C_rn

    # compute the noise floor rate
    # this should have units of 1/time
    M = C_sr * static_params["ppFac"]

    return C_p, C_b, M


def calc_dMag_from_intTime(static_params, coronagraph, target, SNR, intTime):
    """Find the achievable deltaMag for an observation for given SNR and integration
    time.

    Args:
        static_params (dict):
            Dictionary of static parameters. See Cp_Cb_M for details.

        cornagraph (dict):
            Dictionary of coronagraph parameters for this observation.
            See Cp_Cb_M for details.
        target (dict):
            Dictionary of target observing parameters. See Cp_Cb_M for details.
        SNR (float):
            Required signal to noise ratio
        intTime (Quantity):
            Integration time. Units of time.

    Returns:
        float:
            Achievable deltaMag

    """

    # Get values of C_b and M (first output is not needed)
    # Note that the deltaMag input here is entirely arbitrary, as we're not
    # using the C_p output.
    _, C_b, M = Cp_Cb_M(static_params, coronagraph, target, 0)

    # Compute C_p
    # This should have units of 1/time
    C_p = (SNR / 2 / intTime) * (
        SNR + np.sqrt(4 * M**2 * intTime**2 + 4 * C_b * intTime + SNR**2)
    )

    # compute the common factor of A*tau*deltaLam*QE
    # this value should have units of length^3
    A = np.pi * (static_params["D"] / 2) ** 2 * (1 - static_params["obsc"])
    eta = A * static_params["tau"] * static_params["deltaLam"] * static_params["QE"]

    # Compute the deltaMag corresponding to the inputs
    # This should be unitless
    deltaMag = -target["mag_star"] - 2.5 * np.log10(
        C_p / (static_params["F0"] * eta * coronagraph["tau_core"])
    )

    return deltaMag.value


def load_HWO_MissionStars():
    """Load HWO Mission Star List

    Args:
        None

    Returns:
        pandas.DataFrame:
            The target list
    """

    fname = pkg_resources.resource_filename(
        "SSWYieldModelingTutorial",
        os.path.join("data", "HWOStarList_20240526125022.pkl"),
    )

    data = pandas.read_pickle(fname)

    return data


def load_HWO_MissionStars_koMap():
    """Load HWO Mission Star List sample keepout map

    Args:
        None

    Returns:
        tuple:
            koTimes (~astropy.time.Time):
                Absolute MJD mission times from start to end in steps of 1 d
            koMap (~numpy.ndarray(bool)):
                True means a target unobstructed and observable, and False means a
                target unobservable due to obstructions in the keepout zone.
            targetNames (~numpy.ndarray(object)):
                Target names (for consistency checking)


    """

    fname = pkg_resources.resource_filename(
        "SSWYieldModelingTutorial",
        os.path.join("data", "HWOStarList_koMap.pkl"),
    )

    with open(fname, "rb") as f:
        data = pickle.load(f)

    return data["koTimes"], data["koMap"], data["targetNames"]


def check_target_availbility(koTimes, koMap, targetInds, obsStart, obsDurations):
    """Check whether targets are available (out of keepout) for proposed observations

    Args:
        koTimes (~astropy.time.Time):
            Absolute MJD mission times from start to end in steps of 1 d
        koMap (~numpy.ndarray(bool)):
            True means a target unobstructed and observable, and False means a
            target unobservable due to obstructions in the keepout zone.
        targetInds (arraylike):
            Indices of targets to check
        obsStart (~astropy.time.Time):
            Observation start time (absolute time)
        obsDurations (Quantity):
            Duration of observation(s). Time units.


    Returns:
        ~numpy.ndarray(bool):
            Boolean array of target availability. True means the target is available
            for observation and false means it is not. The size of this array should be
            the same as the size of targetInds.

    .. note::

        The size of targetInds and obsDurations must be the same.

    """

    # calculate observation end times:
    obsEnds = obsStart + obsDurations

    # find index  of start time in koTimes
    startInd = np.searchsorted(koTimes.value, obsStart.value)

    # find indices of end times in koTimes
    endInds = np.searchsorted(koTimes.value, obsEnds.value)

    # compute availbility
    avail = np.array(
        [
            np.all(koMap[sInd, startInd : endInds[j]])  # noqa
            for j, sInd in enumerate(targetInds)
        ]
    )

    return avail
