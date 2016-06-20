#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imageproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <bitset>
#include <vector>
#include "2d_algo_fft.hpp"

/** This is "2D Particle image velocimetry, Version " **/

/** mhkhaniitb@gmail.com {PhD}, Guide **/
/** gowthamkuntumalla@gmail.com,krishnasandeepbhogaraju@gmail.com **/
/** Output data file consists of initial point(x1,y1) final point(x2,y2) **/
int main()
{
    // Import images
    int win_size;
    cout<<"Enter the number of images"<<endl;
    int num;
    cin>>num;
    cout<<"Enter the window size A x A (give only A which is a power of 2)"<<endl;//give only 16 or 32 or 64
    cin>>win_size;// window size 'A'
    for(int i=0; i<num; i++) /*all images*/
    {
        /**************Declarations**************/
        Mat image1 = imread("image1.tif",0);
        Mat image2 = imread("image2.tif",0);

        vector< vector <pair<int,int> > > max_coef_point;
        int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image.
        int totrow2= image2.rows,totcol2=image2.cols;
        int nHeight=totrow2;
        int nWidth=totcol2; //setting height of the img
        Mat img(Mat(nHeight, nWidth, CV_8U));//blank black image
        img=Scalar(0);//background color

        /**computations**/

        piv_2d_fft(image1,image2,max_coef_point,win_size,i);// main function

        /**Plot the displacements in a new image**/
       // cv::namedWindow( "Display window", CV_WINDOW_AUTOSIZE );

        for(int c=0; c<totcol1-win_size; c+=win_size/2) // according to win_size.
        {
            for(int r=0; r<totrow1-win_size; r+=win_size/2) // 50% overlap of windows
            {
                Point pt1 = Point(c+win_size/2,r+win_size/2);//here (x,y)=(col,row) point in cartesian coordiante system
                Point pt2 = Point(max_coef_point[r][c].second+win_size/2,max_coef_point[r][c].first+win_size/2);
                //syntax =  void arrowedLine(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int line_type=8, int shift=0, double tipLength=0.1)
                arrowedLine(img, pt1, pt2, Scalar(255), 1, 8, 0, .2);//white
            }
        }
        //cv::imshow( "Display window",img);
        imwrite("result_img_"+to_string(i)+".tif",img);
        //waitKey(0);// wait for a keystroke in the window
    }
    return 0;
}
