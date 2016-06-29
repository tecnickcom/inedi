# iNEDI
*iNEDI (improved New Edge-Directed Interpolation) Image Zooming Algorithm*

[![Donate via PayPal](https://img.shields.io/badge/donate-paypal-87ceeb.svg)](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&currency_code=GBP&business=paypal@tecnick.com&item_name=donation%20for%20inedi%20project)
*Please consider supporting this project by making a donation via [PayPal](https://www.paypal.com/cgi-bin/webscr?cmd=_donations&currency_code=GBP&business=paypal@tecnick.com&item_name=donation%20for%20inedi%20project)*

* **category**    Application
* **package**     \Com\Tecnick\inedi
* **author**      Nicola Asuni <info@tecnick.com>
* **copyright**   2006-2016 Nicola Asuni - Tecnick.com LTD
* **license**     http://www.gnu.org/copyleft/lesser.html GNU-LGPL v3 (see LICENSE.TXT)
* **link**        https://github.com/tecnickcom/inedi
* **version**     1.6.2

## Description

This projects contains a Matlab implementation of the iNEDI (improved New Edge-Directed Interpolation) algorithm to enlarge (zoom) digital images.

This inedi.m function returns an enlarged image by a factor of 2^ZK for both horizontal and vertical directions.

This function implements the technicque described in:
N. Asuni, "iNEDI - Tecnica non lineare per l'interpolazione di immagini," Master's thesis - University of Cagliari, 2007.
Thesis [ITA]: http://www.tecnick.com/pagefiles/appunti/iNEDI_tesi_Nicola_Asuni.pdf
Presentation [ITA]: http://www.tecnick.com/pagefiles/appunti/iNEDI_presentazione_Nicola_Asuni.pdf



## USAGE
-
[EIM] = inedi(IM)
[EIM] = inedi(IM, ZK)
[EIM] = inedi(IM, ZK, MT)
[EIM] = inedi(IM, ZK, MT, ML)
[EIM] = inedi(IM, ZK, MT, ML, BT)
[EIM] = inedi(IM, ZK, MT, ML, BT, BS)
[EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ)
[EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ, VR)


## INPUT

IM : source image
ZK : power of the zoom factor (default = 1)
     the image enlargement on vertical and horizontal direction is
     2^ZK; the final image size will be (SIZE * 2^ZK) - (2^ZK - 1)
MT : maximum radius of the local circular window (default = 9)
ML : minimum radius of the local circular window (default = 3)
BT : brightness tollerance (default = 32)
     Y points are choosen from the ones inside a circular area around
     the pixel to be interpolated. These pixels must have a
     brightness between a range defined by
     the minimum and maximum brightness values of the 4 neighbors
     pixels of the pixel to interpolate. This range is extended by
     +/- BT value.
BS : maximum brigthess difference for smooth areas (default = 4)
     This parameter is compared to the difference between the maximum
     and minimum brightness value of 4 neighbors pixels to
     discriminate between smooth areas (bicubic interpolation) and
     edge areas (iNEDI interpolation).
SZ : number of image bits per layer (default = 8)
VR : verbose mode (default = false)
     if true prints some information during calculation


## OUTPUT

EIM : enlarged image


## Examples

Please check the inediexample.m file on how to use this function.


## Notes

This implementation is not intended to be used in a production environment. The main purpose of this script is to clearly show how this technique works. Better performaces could be obtained using a compiled version or rewriting this technique using a low-level programming language.

## KEYWORDS
iNEDI, NEDI, image, zooming, magnification, upsizing, resampling, resolution enhancement, interpolation, super-resolution, covariance-based adaptation, geometric regularity, matlab, octave.


