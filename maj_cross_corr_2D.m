% clear all;
clc;

Image_A = imread('image1_a.tif');
Image_B = imread('image1_b.tif');
 
Image_A = im2double(Image_A);
Image_B = im2double(Image_B);


win_size  = input('Enter the window size = ');
step_size = input('Enter the step size   = ');


l= 1; m = 1; l1 = 0; m1 = 0;
iter_lim = size(Image_A);
for i = 1: ((iter_lim(1)/step_size)-1)
    for j =1: ((iter_lim(2)/step_size)-1)
        
        A = Image_A(l:l1+win_size, m:m1+win_size);
        B = Image_B(l:l1+win_size, m:m1+win_size);
        % _________________________Cross Correlation ______________________________
        %--------------------------------------------------------------------------
        chk1 = max(max(A))/2;
        chk2 = max(max(B))/2;
        %..........................Mean of Image...................................
        avgA = mean(mean(A));
        avgB = mean(mean(B));
        if avgA > 0  
        %....................Standard Deviation of the Image.......................
        stdA = std(std(A));
        stdB = std(std(B));
        %....................Subtractiong Mean from Image .........................
        A = A - avgA;
        B = B - avgB;
        
        %_____________FFT based algorithm for cross correlation of window _________
        %--------------------------------------------------------------------------
        Ga = fft2(A);
        Gb = fft2(B);
        R = conj(Ga).*Gb;
        r = ifft2(R);
        
        CII = r/(stdA*stdB);
        CII = fftshift(CII);
      
        
        %..............Padding the Cross Corr Matrix...............................
        size_CII = size(CII);
        CII_pad = zeros(size_CII(1)+10, size_CII(2)+10);
        CII_pad(6:end-5, 6:end-5) = CII;
        CII = CII_pad;
        CII(isnan(CII)) = 100;
        size_CII = size(CII);
        
        %....................Peak detection of cross correlated Image..............
        
        %........................ Max Value Tracking block ........................
        
        for i_track = 1:size_CII(1)
            for j_track = 1:size_CII(2)
                if CII(i_track, j_track) == max(max(CII))
                    xa = i_track;
                    ya = j_track;
                end
            end
        end
        
        %*************Check Block *********
        
        
        %...........................Peak Centroid .................................
        %...Three point estimator for determining displacement.....................
        %......................from the correlation data...........................
        
        Xo_num   = (xa-1)*CII((xa-1), ya) + xa*CII(xa, ya) + (xa+1)*CII((xa+1), ya);
        Xo_denum = CII((xa-1), ya) + CII(xa, ya) + CII((xa+1), ya);
        Yo_num   = (ya-1)*CII(xa, (ya-1))+ ya*CII(xa, ya) + (ya+1)*CII(xa, (ya+1));
        Yo_denum = CII(xa, (ya-1))+ CII(xa, ya) + CII(xa, (ya+1));
        
        Xo = Xo_num/Xo_denum;
        Yo = Yo_num/Yo_denum;
        
        %..................Cordinate Allocation ..................................
        X_cord = Yo;
        Y_cord = Xo;
        
        
        %_______________________Displacement Block_________________________________
        %________________Positions of peaks on Ist Image___________________________
        %__________________Image A cordinates of X and Y___________________________
        %...............Determining the Centre of Correlation plane......
        
        X_cord_CII0 = size_CII(2)/2;
        Y_cord_CII0 = size_CII(1)/2;
        
        %..................Displacement of particle ..............................
        x_disp = (X_cord - X_cord_CII0-1);
        y_disp = (Y_cord - Y_cord_CII0-1);
        
        %##########################################################################
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SAVING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        u(i, j) = x_disp;
        v(i, j) = y_disp;
        
        %.................Changing steps............................................
        end
        m = j*step_size; m1 = m;
        
    end
    l = i*step_size; l1 = l; 
    m = 1; m1 = 0;
end
