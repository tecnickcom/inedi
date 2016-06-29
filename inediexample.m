function [] = inediexample()
    % =====================================================================
    % File name   : inediexample.m
    % File Type   : m-file (script file for Matlab or Octave)
    % Requirements: Image Processing Toolbox for Matlab or 
    %               ImageMagick for Octave
    % Begin       : 2006-07-07
    % Last Update : 2007-08-04
    % Author      : Nicola Asuni
    % Description : Example script for inedi.m usage.
    % Copyright   : 2006-2016 Nicola Asuni - Tecnick.com LTD
    % License     : GNU GENERAL PUBLIC LICENSE v.3
    %               http://www.gnu.org/copyleft/gpl.html
    % Version     : 1.6.2
    % =====================================================================

    % DESCRIPTION
    % --------------------
    % This is an example script that shows you how to use the inedi.m function.
    %
    % The inedi.m function returns an enlarged image by a factor of 2^ZK for
    % both horizontal and vertical directions. 

    % USAGE
    % --------------------
    % inediexample

    % Example
    % --------------------
    % >> inediexample

    % ---------------------------------------------------------------------

    % Help
    disp('iNEDI function by Nicola Asuni.');
    disp('This is an example script that shows you how to use the inedi.m function.');
    disp('Please check the documentation inside this file for further information.');

    % --- set parameters ---

    % power of the zoom factor
    % the image enlargement on vertical and horizontal direction is 2^ZK 
    ZK = 1;
    
    % maximum radius of the local circular window
    MT = 6;
    
    % minimum radius of the local circular window
    ML = 3;
    
    % brightness tollerance
    %      Y points are choosen from the ones inside a circular area around
    %      the pixel to be interpolated. These pixels must have a
    %      brightness between a range defined by
    %      the minimum and maximum brightness values of the 4 neighbors
    %      pixels of the pixel to interpolate. This range is extended by
    %      +/- BT value.
    BT = 16;
    
    % maximum brigthess difference for smooth areas
    %      This parameter is compared to the difference between the maximum
    %      and minimum brightness value of 4 neighbors pixels to
    %      discriminate between smooth areas (linear interpolation) and
    %      edge areas (NEDI-based interpolation).
    BS = 8;
    
    % number of image bits per layer
    SZ = 8;
    
    % verbose mode
    %      if true prints some information during calculation
    VR = true;

    % --- end set parameters ---

    % load test image
    IM = imread('testimage.png'); 

    % display original image
    figure, imshow(IM);

    % get the enlarged image
    [EI] = inedi(IM, ZK, MT, ML, BT, BS, SZ, VR);

    % save image to disk
    %imwrite(EI,'testimage_inedi.tif','tif');

    % display enlarged image
    figure, imshow(EI);

% === EOF =================================================================
