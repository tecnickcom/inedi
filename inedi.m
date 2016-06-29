function [EI] = inedi(IM, ZK, MT, ML, BT, BS, SZ, VR)
    % =====================================================================
    % File name   : inedi.m
    % File Type   : m-file (script file for Matlab or Octave)
    % Begin       : 2006-07-07
    % Last Update : 2007-08-31
    % Author      : Nicola Asuni
    % Description : This function returns an enlarged image by a factor of
    %               2^ZK for both horizontal and vertical directions.
    % Copyright   : 2006-2016 Nicola Asuni - Tecnick.com LTD
    % License     : GNU GENERAL PUBLIC LICENSE v.3
    %               http://www.gnu.org/copyleft/gpl.html
    % Version     : 1.6.2
    % =====================================================================
    
    % DESCRIPTION
    % --------------------
    % This function returns an enlarged image by a factor of 2^ZK for both
    % horizontal and vertical directions.
    %
    % iNEDI is an acronym for "improved New Edge-Directed Interpolation".
    % This function implements the technicque described on:
    % N. Asuni, "iNEDI - Tecnica adattativa per l'interpolazione di 
    % immagini," Master's thesis - University of Cagliari, 2007.
    % Thesis [ITA]:  
    % http://www.tecnick.com/pagefiles/appunti/iNEDI_tesi_Nicola_Asuni.pdf
    % Presentation [ITA]: 
    % http://www.tecnick.com/pagefiles/appunti/iNEDI_presentazione_Nicola_Asuni.pdf
    % 
    % This technique is a deep improvement of the original NEDI technique
    % described on:
    % X. Li et. al., "new edge-directed interpolation". IEEE Trans. on
    % Image Processing, Vol.10, No.10, October 2001, pp.1521-1527
    
    
    % KEYWORDS
    % --------------------
    % iNEDI, NEDI, image, zooming, magnification, upsizing, resampling,
    % resolution enhancement, interpolation, super-resolution,
    % covariance-based adaptation, geometric regularity, matlab, octave.
    
    
    % WARNING
    % --------------------
    % This function is slow because of high computational complexity.
    
    
    % USAGE
    % --------------------
    % [EIM] = inedi(IM)
    % [EIM] = inedi(IM, ZK)
    % [EIM] = inedi(IM, ZK, MT)
    % [EIM] = inedi(IM, ZK, MT, ML)
    % [EIM] = inedi(IM, ZK, MT, ML, BT)
    % [EIM] = inedi(IM, ZK, MT, ML, BT, BS)
    % [EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ)
    % [EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ, VR)
    
    
    % INPUT
    % --------------------
    % IM : source image
    % ZK : power of the zoom factor (default = 1)
    %      the image enlargement on vertical and horizontal direction is
    %      2^ZK; the final image size will be (SIZE * 2^ZK) - (2^ZK - 1)
    % MT : maximum radius of the local circular window (default = 9)
    % ML : minimum radius of the local circular window (default = 3)
    % BT : brightness tollerance (default = 32)
    %      Y points are choosen from the ones inside a circular area around
    %      the pixel to be interpolated. These pixels must have a
    %      brightness between a range defined by
    %      the minimum and maximum brightness values of the 4 neighbors
    %      pixels of the pixel to interpolate. This range is extended by
    %      +/- BT value.
    % BS : maximum brigthess difference for smooth areas (default = 4)
    %      This parameter is compared to the difference between the maximum
    %      and minimum brightness value of 4 neighbors pixels to
    %      discriminate between smooth areas (bicubic interpolation) and
    %      edge areas (iNEDI interpolation).
    % SZ : number of image bits per layer (default = 8)
    % VR : verbose mode (default = false)
    %      if true prints some information during calculation
    
    
    % OUTPUT
    % --------------------
    % EIM : enlarged image
    
    
    % Examples
    % --------------------
    % Please check the inediexample.m file on how to use
    % this function.
    
    
    % Notes
    % --------------------
    % This implementation is not intended to be used in a production
    % environment. The main purpose of this script is to clearly show how
    % this technique works. Better performaces could be obtained using a
    % compiled version or rewriting this technique using a low-level
    % programming language.
    
    % ---------------------------------------------------------------------
    
    % Some initial tests on the input arguments
    
    if (nargin < 1)
        disp('iNEDI function by Nicola Asuni.');
        disp('This function returns an enlarged image.');
        disp('Usage:');
        disp('[EIM] = inedi(IM)');
        disp('[EIM] = inedi(IM, ZK)');
        disp('[EIM] = inedi(IM, ZK, MT)');
        disp('[EIM] = inedi(IM, ZK, MT, ML)');
        disp('[EIM] = inedi(IM, ZK, MT, ML, BT)');
        disp('[EIM] = inedi(IM, ZK, MT, ML, BT, BS)');
        disp('[EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ)');
        disp('[EIM] = inedi(IM, ZK, MT, ML, BT, BS, SZ, VR)');
        disp('Where:');
        disp('IM : source image');
        disp('ZK : power of the zoom factor [zoom factor = 2^ZK] (default = 1)');
        disp('MT : maximum radius of the local circular window (default = 9)');
        disp('ML : minimum radius of the local circular window (default = 3)');
        disp('BT : brightness tollerance (default = 32)');
        disp('BS : maximum brigthess difference for smooth areas (default = 4)');
        disp('SZ : number of image bits per layer (default = 8)');
        disp('VR : verbose mode (default = false)');
        disp('EI : enlarged image');
        EI = [];
        return;
    end
    
    % assign default values
    if (nargin > 8)
        error('Too many arguments');
    end
    if (nargin < 8)
        VR = false;
    end
    if (nargin < 7)
        SZ = 8;
    end
    if (nargin < 6)
        BS = 4;
    end
    if (nargin < 5)
        BT = 32;
    end
    if (nargin < 4)
        ML = 3;
    end
    if (nargin < 3)
        MT = 9;
    end
    if (nargin < 2)
        ZK = 1;
    end
    
    % ---------------------------------------------------------------------
    
    % check parameters
    
    if (ZK < 1)
        error('ZK must be a positive integer');
    end
    if (BS < 0)
        error('BS must be a positive integer');
    end
    if (BT < 0)
        error('BT must be a positive integer');
    end
    if (ML < 2)
        error('ML must be greater than or equal to 2');
    end 
    if (ML > MT)
        error('MT must be greater than or equal to ML');
    end
    
    % ---------------------------------------------------------------------
    
    if(VR)
        t1 = cputime;
    end
    
    % max intensity pixel value
    MAXPIX = (2^SZ) - 1;
    
    % get number of image colors
    IDIM = ndims(IM);
    if (IDIM == 3)
        % number of colors
        CL = size(IM,3);
    elseif (IDIM == 2)
        CL = 1;
    else
        error('Unrecognized image type, please use RGB or greyscale images');
    end
    
    % initialize some parameters used for bicubic interpolation
    v = [-1 0 1 2];
    V = zeros(16,16);
    k=1;
    for i=1:4
        for j=1:4
            % Vandermonde matrix
            V(k,:) = [((v(i)^3)*(v(j)^3)),((v(i)^3)*(v(j)^2)),((v(i)^3)*(v(j)^1)),((v(i)^3)*(v(j)^0)),((v(i)^2)*(v(j)^3)),((v(i)^2)*(v(j)^2)),((v(i)^2)*(v(j)^1)),((v(i)^2)*(v(j)^0)),((v(i)^1)*(v(j)^3)),((v(i)^1)*(v(j)^2)),((v(i)^1)*(v(j)^1)),((v(i)^1)*(v(j)^0)),((v(i)^0)*(v(j)^3)),((v(i)^0)*(v(j)^2)),((v(i)^0)*(v(j)^1)),((v(i)^0)*(v(j)^0))];
            k = k+1;
        end
    end
    % bicubic coefficients
    CXY = double(zeros(16,1));
    % 1/(1/2^i * 1/2^j)
    BXY = [64;32;16;8;32;16;8;4;16;8;4;2;8;4;2;1];
    
    % variable used to switch between iNEDI and Bicubic
    USEINEDI = false;
    
    % number of padding pixels
    MP = 2 + MT;
    
    % Circular kernel masks (KM) to be used for the selection of
    % Y points.
    % The use of a circular mask instead of the square window used in
    % original NEDI technique reduces the number of calculations and the
    % presence of artifacts in case of large windows.
    % The use of multiple masks with different sizes allows to account for
    % both lower and higher frequencies regions.
    KM = [];
    k = 1;
    for i = -MT+1:MT
        for j = -MT+1:MT
            % get pixel distance from center
            RAD = (sqrt(((i-0.5).^2) + ((j-0.5).^2)));
            if (RAD <= MT)
               KM(k,1) = i; % row coordinate (y)
               KM(k,2) = j; % column coordinate (x)
               KM(k,3) = ceil(RAD); % pixel distance from center
               k = k + 1;
            end
        end
    end
    
    % ---------------------------------------------------------------------
    
    % for each color
    for CID = 1:CL
        
        if(VR)
            fprintf('\n[%8.3f sec] LAYER: %02d\n', cputime-t1, CID);
        end
        
        % intensity image relative to the current color
        % (the image is traslated to reduce computational errors)
        IMG = double(IM(:,:,CID)) + MAXPIX;
        
        % the image is enlarged by scaling factor 2^ZK at each cycle
        for ZF = 1:ZK
            
            if(VR)
                fprintf('[%8.3f sec]    ZF: %02d\n', cputime-t1, ZF);
            end
            
            % image size
            [m,n] = size(IMG);

            % size of the expanded image (removing padding)
            fm = (2 * (m - 1)) + 1;
            fn = (2 * (n - 1)) + 1;

            % calculate interpolated points in two steps
            % by rotating the result image of each step by pi/4
            for s = 0:1
                
                if(VR)
                    fprintf('[%8.3f sec]        PHASE: %02d\n', cputime-t1, s);
                end
                
                % size of the expanded image rotated by pi/4
                mm = m + n - 1;
                nn = mm;
                
                % display current image
                %figure, imshow(uint8(round(IMG)));
                
                % creates a void expanded image
                IMGEXP = double(zeros(mm,nn));
                
                % starting coordinates to store expanded image on a rotated
                % matrix 
                SY = 1;
                SX = m;

                % padding the image
                IMG = [repmat(IMG(1,:),MP,1);IMG;repmat(IMG(m,:),MP,1)];
                IMG = [repmat(IMG(:,1),1,MP),IMG,repmat(IMG(:,n),1,MP)];
                
                % image size after padding
                [m,n] = size(IMG);

                for i = (1+MP):(m-MP) % for each row
                   
                    % coordinates for the new pi/4 rotated matrix
                    ii = SY;
                    jj = SX;
                    
                    for j = (1+MP):(n-MP) % for each column

                        % arrange pixel on expanded lattice
                        % (insert results on pi/4 rotated matrix)
                        IMGEXP(ii,jj) = IMG(i,j);
                        
                        % fill 4 empty triangles caused by image rotation
                        if (s == 0)
                            % top -> top-right
                            if (i == 1+MP) && (j > 1+MP)
                                for p = ii-1:-1:1
                                    IMGEXP(p,jj) = (IMGEXP(p,jj-1) + IMGEXP(p+1,jj)) / 2;
                                end
                            end
                            % left -> top-left
                            if (i > 1+MP) && (j == 1+MP)
                                for p = ii-1:-1:1
                                    IMGEXP(p,jj) = (IMGEXP(p,jj+1) + IMGEXP(p+1,jj)) / 2;
                                end
                            end
                            % bottom -> bottom-left
                            if (i == m-MP) && (j > 1+MP)
                                for p = jj-1:-1:1
                                    IMGEXP(ii,p) = (IMGEXP(ii-1,p) + IMGEXP(ii,p+1)) / 2;
                                end
                            end
                            % right -> bottom-right
                            if (i > 1+MP) && (j == n-MP)
                                for p = jj+1:nn
                                    IMGEXP(ii,p) = (IMGEXP(ii-1,p) + IMGEXP(ii,p-1)) / 2;
                                end
                            end
                        end
                        
                        % if we are not on the last row or last column before padding
                        if (i < m-MP) && (j < n-MP)
                            
                            % get min and max values of 4 neighbors
                            PMIN = min(min(IMG(i:i+1,j:j+1)));
                            PMAX = max(max(IMG(i:i+1,j:j+1)));
                            
                            % max intensity difference
                            PDIF = PMAX - PMIN;
                            
                            % NEDI technique is used only for pixels near
                            % edges
                            if  (PDIF > BS)
                                
                                % define an admissibility range for Y points
                                RNG = min(BT,PDIF);
                                LMIN = max(MAXPIX, PMIN - RNG);
                                LMAX = min(2*MAXPIX, PMAX + RNG);
                                
                                % build a mask to select only connected
                                % pixels to the 4 neighbors of the pixel 
                                % to interpolate
                                SELMASK = IMG(i-MT+1:i+MT,j-MT+1:j+MT);
                                SELMASK = logical((SELMASK >= LMIN) & (SELMASK <= LMAX));
                                %SELMASK = bwselect(SELMASK, [MT,MT+1,MT,MT+1], [MT,MT,MT+1,MT+1], 8);
                                SELMASKA = bwselect(SELMASK, MT, MT, 8);
                                SELMASKB = bwselect(SELMASK, MT+1, MT, 8);
                                SELMASKC = bwselect(SELMASK, MT, MT+1, 8);
                                SELMASKD = bwselect(SELMASK, MT+1, MT+1, 8);
                                SELMASK = logical(SELMASKA | SELMASKB | SELMASKC | SELMASKD); 
                                
                                % initialize previous MSE (Mean Square
                                % Error)
                                PREVMSE = Inf;
                                
                                % The use of multiple masks with different
                                % sizes allows to account for both lower
                                % and higher frequencies regions.
                                % We always start from the smallest mask
                                % (high frequencies).
                                for MRAD = ML:MT % mask radius
                                    
                                    % The Y vector contains the pixels around the
                                    % pixel to be interpolated. The 4 neighbors of
                                    % each of these pixels (see schema below) are
                                    % reported as columns of C matrix.
                                    %  c11   .   c12 
                                    %   .    y1   . 
                                    %  c13   .   c14
                                    Y = [];
                                    C = [];
                                    yid = 0; % index for Y vector and C matrix
                                    
                                    % get the mask IDs inside MRAD radius
                                    KMID = find(KM(:,3) <= MRAD);
                                                                        
                                    % number of pixels inside the circular mask
                                    KMSIZE = numel(KMID);
    
                                    % for each pixel inside the circular mask
                                    % centered on the pixel to interpolate
                                    for q = 1:KMSIZE
                                        
                                        % get pixels inside the circular area
                                        h = i + KM(KMID(q),1); % row
                                        k = j + KM(KMID(q),2); % column
                                        
                                        if (SELMASK(KM(KMID(q),1)+MT,KM(KMID(q),2)+MT))
                                            % four diagonal neighbors
                                            CTMP = [IMG(h-1,k-1) IMG(h-1,k+1) IMG(h+1,k-1) IMG(h+1,k+1)];
                                            % select only pixels that
                                            % belongs to an edge
                                            if ((max(CTMP) - min(CTMP)) > BS)
                                                yid = yid + 1;
                                                % Y vector
                                                Y(yid,1) = IMG(h,k);
                                                % C matrix
                                                C(yid,:) = CTMP;
                                            end
                                        end
                                    end
                                    
                                    % if Y vector has more than four points
                                    if (yid > 4)
                                        
                                        % calculates interpolation coefficients
                                        % (Moore-Penrose pseudoinverse is used
                                        % because the matix C'*C could be bad
                                        % conditioned.
                                        APP = pinv(C) * Y;
                                        
                                        % Mean Square Error
                                        MSE = sum(sum(((C * APP) - Y ).^2)) / yid;
                                        
                                        if (MSE >= PREVMSE)
                                            % the previous calculated
                                            % coefficients are better or
                                            % equal than current.
                                            break;
                                        end
                                        
                                        % store current MSE for next
                                        % comparison
                                        PREVMSE = MSE;
                                        
                                        % reshape the coefficients in a 2x2
                                        % matrix
                                        AP = [APP(1),APP(2);APP(3),APP(4)];
                                        USEINEDI = true;
                                    else
                                        % bicubic interpolation is used
                                        USEINEDI = false;
                                        % exit from the for mask loop
                                        break;
                                    end
                                end % end of for mask radius
                            else
                                % bicubic interpolation is used
                                USEINEDI = false;
                            end
                            
                            % calculates the missing point value and insert
                            % result on the pi/4 rotated matrix
                            if (USEINEDI)
                                IMGEXP(ii+1,jj) = sum(sum(AP .* IMG(i:i+1,j:j+1)));
                            else
                                % bicubic interpolation
                                CXY = V \ reshape(IMG(i-1:i+2,j-1:j+2)',16,1);
                                IMGEXP(ii+1,jj) = sum(CXY ./ BXY);
                            end
                            
                            % adjust the new value if it's out of the
                            % specified range
                            if (IMGEXP(ii+1,jj) < PMIN)
                                IMGEXP(ii+1,jj) = PMIN;
                            end
                            if (IMGEXP(ii+1,jj) > PMAX)
                                IMGEXP(ii+1,jj) = PMAX;
                            end
                        end
                        ii = ii + 1;
                        jj = jj + 1;
                    end
                    SY = SY + 1;
                    SX = SX - 1;
                end
                IMG = IMGEXP;
                [m,n] = size(IMG);
            end % end of (for s = 0:1)

            % crop image padding
            SR = ((fm - 1) / 2);
            SC = ((fn - 1) / 2);
            IMGEXP = IMGEXP(SR+1:mm-SR,SC+1:nn-SC);

            % restore the right image orientation
            IMG = rot90(IMGEXP);
            
        end % end of (for ZF = 1:ZK)
        
        % translate the image to original scale
        IMG = IMG - MAXPIX;
        
        % convert data to best integer type
        if (SZ > 32)
            EI(:,:,CID) = uint64(round(IMG));
        elseif (SZ > 16)
            EI(:,:,CID) = uint32(round(IMG));
        elseif (SZ > 8)
            EI(:,:,CID) = uint16(round(IMG));
        else
            EI(:,:,CID) = uint8(round(IMG));
        end
    end % end of (for CID = 1:CL)
    
    if(VR)
        fprintf('[%8.3f sec] END\n', cputime-t1);
    end

% === EOF =================================================================
