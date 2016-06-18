#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imageproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <string>
#include "2d_algo_fft.hpp"

/** This is "2D Particle image velocimetry, Version " **/

/** mhkhaniitb@gmail.com {PhD}, Guide **/
/** gowthamkuntumalla@gmail.com,krishnasandeepbhogaraju@gmail.com **/
/** Output data file consists of initial point(x1,y1) final point(x2,y2) **/
int main()
{
    // Import images
    int win_size;
    cout<<"Enter the window size A x A (give only A which is a power of 2)"<<endl;//give only 16 or 32 or 64
    cin>>win_size;// window size 'A'
    for(/*all images*/)
    {
        /**************Declarations**************/
        Mat image1 = imread("image1.tif",0);
        Mat image2 = imread("image2.tif",0);
        Mat img(Mat(nHeight, nWidth, CV_8U));//blank black image
        img=Scalar(0);//background color
        vector< vector <pair<int,int> > > max_coef_point;
        int totrows1= image1.rows,totcols1=image1.cols; //opencv functions to get the rows and columns of image.
        int totrows2= image2.rows,totcols2=image2.cols;

        /**computations**/

        piv_2d_fft(image1,image2,max_coef_point,win_size);// main function

        /**Plot the displacements in a new image**/
        cv::namedWindow( "Display window", CV_WINDOW_AUTOSIZE );
        for(int c=0; c</**/; c+=/**/)
        {
            for(int r=0; r</**/; r+=/**/)
            {
                Point pt1 = Point(c+/**/,r+/**/);//here (x,y)=(col,row) point in cartesian coordiante system
                Point pt2 = Point(max_coef_point[r][c].second+/**/,max_coef_point[r][c].first+/**/);
                //syntax =  void arrowedLine(Mat& img, Point pt1, Point pt2, const Scalar& color, int thickness=1, int line_type=8, int shift=0, double tipLength=0.1)
                arrowedLine(img, pt1, pt2, Scalar(255), 1, 8, 0, .2);//white
            }
        }
        //cv::imshow( "Display window",img);
        imwrite("result_img"+to_string(/**i**/)+".tif",img);
        //waitKey(0);// wait for a keystroke in the window
    }
    return 0;
}
