import numpy as np
import matplotlib.pyplot as plt


# Defining the 3 wavelengths of interest measured in microns
lambd_c = 590e-3
lambd_f = 532e-3
lambd_d = (lambd_c + lambd_f)/2


# Defining a function for the refractive index as a function of wavelength in microns
def n(x):
    n_temp = np.sqrt(1 + 0.6961663*x**2/(x**2 - 0.0684043**2) + 0.4079426*x **
                     2/(x**2 - 0.1162414**2) + 0.8974794*x**2/(x**2 - 9.896161**2))
    return n_temp

# Defining the function of the Abbé number for the refractive lens


def Vr(f, d, c):
    nf = n(f)
    nd = n(d)
    nc = n(c)
    Vr_temp = (nd - 1)/(nf - nc)
    return Vr_temp

# Defining the function of the Abbé number for the diffractive lens


def Vd(f, d, c):
    '''Takes the wavelengths of the wavespectrum in microns'''
    Vd_temp = (d)/(f - c)
    return Vd_temp

# Defining the function for the focal length of the refractive lens element in an achromat doublet


def fr(Vr, Vd, F):
    '''Takes the Abbé numbers of the refractive and diffractive lens elements and the total focal length of the achromat doublet'''
    fr_temp = F*(Vr - Vd)/Vr
    return fr_temp

# Defining the function for the focal length of the diffractive lens element in an achromat doublet


def fd(Vr, Vd, F):
    '''Takes the Abbé numbers of the refractive and diffractive lens elements and the total focal length of the achromat doublet'''
    fd_temp = F*(Vd - Vr)/Vd
    return fd_temp

#


def R1(fr, nd):
    '''Takes the focal length of the refractive lens element in an achromat doublet in meters and the refractive index for the lens element'''
    R1_temp = fr*(nd - 1)
    return R1_temp

# def lambda_min(d,na):


def lambda_min(R, d, f):
    '''Takes the units for the wavelength d in microns'''
    M_temp = M(R, d, f)
    lambda_min_temp = d*(f+M_temp*d)/R
    # lambda_min_temp = d*1e-6 / na
    return lambda_min_temp

# def M(R,d,na):


def M(R, d, f):
    '''Takes the units for the wavelength d in microns and f and R in meters'''
    # M_temp = R/(d*1e-6) * (1-(1-na**2)**0.5)/na
    M_temp = (np.sqrt(R**2 + f**2) - f)/d
    return M_temp


def na(R, F):
    '''Calculates the numerical aperture of the lens from it's focal length and radius'''
    na_temp = R/F
    return na_temp


def fr_calc(lam, R1):
    '''Calculates the focal length of the refractive lens element in an achromat doublet for use in testing of the chromatic aberration'''
    fr_temp = 1/((n(lam) - 1)/R1)
    return fr_temp


def fd_calc(lam, M, R, lam_min):
    '''Calculates the focal length of the diffractive lens element in an achromat doublet for use in testing of the chromatic aberration'''
    f_temp = R*lam_min/(lam*1e-6) - M*lam*1e-6
    return f_temp


def f_total(fd, fr):
    '''Calculates the total focal length of the achromat doublet'''
    f_temp = 1/(1/fd + 1/fr)
    # print(f'f_total = {f_temp}')
    return f_temp


class lens():
    '''Class for keeping track of the properties of a lens'''

    def __init__(self, f, R, achromat=True, lambd_f=lambd_f, lambd_d=lambd_d, lambd_c=lambd_c):
        '''Takes the focal length and radius of the lens in meters.
        The lens will be made from the theory of achromatic lens design unless achromat is set to False'''
        # Setting the properties of the lens
        self.f = f
        self.R = R
        self.na = na(R, f)
        # The values for the 3 wavelengths are required to calculate the Abbé numbers.
        # By default the values have been set as defined in the top of the document
        self.Vr = Vr(lambd_f, lambd_d, lambd_c)
        self.Vd = Vd(lambd_f, lambd_d, lambd_c)

        # If the lens is an achromat, the properties of the lens will be calculated using the achromatic lens design theory
        if achromat:
            self.fr = fr(self.Vr, self.Vd, self.f)
            self.R1 = R1(self.fr, n(lambd_d))
            self.fd = fd(self.Vr, self.Vd, self.f)
            # self.nad = na(self.R,self.fd)
            self.lambda_min = lambda_min(R=self.R, d=lambd_d, f=self.fd)
            self.M = M(R=self.R, d=lambd_d, f=self.fd)
        # If the lens is not set as an achromat, the parameters will be set from the basis that the radius of curvature for the refractive lens element
        # is 5 times the radius of the lens, so to keep to the thin lens approximation
        else:
            self.R1 = 5*R
            self.fr = fr_calc(lambd_d, self.R1)
            self.fd = 1/(1/f - 1/self.fr)
            # self.nad = na(self.R,self.fd)
            self.lambda_min = lambda_min(R=self.R, d=lambd_d, f=self.fd)
            self.M = M(R=self.R, d=lambd_d, f=self.fd)

    def print_vals(self):
        print(f'f = {self.f*1000} mm')
        print(f'R = {self.R*1000} mm')
        print(f'Vr = {self.Vr}')
        print(f'Vd = {self.Vd}')
        print(f'R1 = {self.R1*1000} mm')
        print(f'NA = {self.na}')
        print(f'fr = {self.fr*100} cm')
        print(f'fd = {self.fd*100} cm')
        # print(f'NA_d = {self.nad}')
        print(f'lambda_min = {self.lambda_min*1e6} um')
        print(f'M = {self.M}')

    def chromatic_focals(self, lambd_c=lambd_c, lambd_f=lambd_f):
        '''Returns the focal lengths at the two wavelengths'''
        frc = fr_calc(lambd_c, self.R1)
        fdc = fd_calc(lambd_c, self.M, self.R1, self.lambda_min)
        fc_total = f_total(frc, fdc)

        frf = fr_calc(lambd_f, self.R1)
        fdf = fd_calc(lambd_f, self.M, self.R1, self.lambda_min)
        ff_total = f_total(frf, fdf)

        return fc_total, ff_total

    def chromatic_aberration(self):
        '''Returns the chromatic aberration of the lens focal lengths at the lens' specified wavelengths'''
        fc_total, ff_total = self.chromatic_focals

        chromatic_aberration = fc_total - ff_total
        return chromatic_aberration
