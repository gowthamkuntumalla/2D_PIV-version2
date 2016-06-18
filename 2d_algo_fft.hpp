#ifndef _2D_ALGO_FFT_HPP_
#define _2D_ALGO_FFT_HPP_
#include <cmath>
#include <bitset>
using namespace std;
using namespace cv;

/*********************** Auxiliary Functions ***********************/
double avg(const int x2,const int y2, Mat image,const int x1,const int y1)//average of sub matrix
{
    double sum=0.0;
    for (int i=x1; i<x2; i++) //window starting at (x,y)
    {
        for (int j=y1; j<y2; j++) // <i><x><col> equivalent and similar
        {
            double a = (int)image.at<uchar> (j,i);//(row,column)
            sum+=a;
        }
    }
    sum/=((x2-x1)*(y2-y1));
    return sum;
}
double sd(const int x2,const int y2, Mat image,const int x1,const int y1)//standard deviation of sub matrix
{
    double var= 0.0;//variance
    double aver = avg(x2,y2,image,x1,y1);
    //myfile<<aver<<endl;
    for (int i=x1; i<x2; i++) //window starting at (x,y)
    {
        for (int j=y1; j<y2; j++)
        {
            double a = (int)image.at<uchar> (j,i);//(row,column)
            var+=((a-aver)*(a-aver));
        }
    }
    var=sqrt(var)/((x2-x1)*(y2-y1));
    return var;
}
void max_coef(vector< vector<double> > t,const int wincol,const int winrow,int& max_x,int& max_y )
{
    int a=-1;// all intensities are greater than this
    int m=max_x,n=max_y;
    for(int i=m; i<wincol+m; i++)
    {
        for(int j=n; j<winrow+n; j++)
        {
            if(t[j][i]>a)//(row,column)
            {
                a=t[j][i];//(row,column)
                max_x=i;
                max_y=j;
            }
        }
    }
    return;
}
/** Fast Fourier Transform **/
void fft(Mat img,vector< vector <pair<double,double> > > &ft_img,int c,int r)//image to fourier transformed 2D vector
{
    bitset</** 4 for 16, 5 for 32 etc**/> b1{/**number**/};
}
void inv_fft(vector< vector <pair<double,double> > > &ft_img,vector< vector <int> > &inv_ft_img)//fourier transformed 2D vector to image(2D vector)
{

}
void piv_2d_fft(Mat image1, Mat image2,vector< vector <pair<int,int> > > &max_points,int win_size)
{
    /******This is valid for displacements of maximum A/3 where A is window_Size  ******/
    ofstream myfile; // like 'cout' it outputs to a file
    int initial_value = 0;
    double avg1=0,avg2=0;
    double sd1=0, sd2=0;
    int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image.
    int totrow2= image2.rows,totcol2=image2.cols;

    vector< vector <pair<double,double> > > /**/; //Fourier transform For image1
    vector< vector <pair<double,double> > > /**/; //For image2
    vector< vector <pair<double,double> > > /**/; //For product (I1)*(I2')
    /**/ = vector<vector<double> >(totrow1, vector<double>(totcol1,initial_value));// initialization
    /** Loop over entire image**/
    for(int c=0; c</**/; c+=/**/) // according to win_size
    {
        for(int r=0; r</**/; r+=/**/)
        {
            /**do FFT**/
            //bit reversal int b=a>>1;
            /**calculate i1*(i2')**/
            /**do inverse FFT**/
            /**Save the data in a "data_i.txt" file**/
        }
    }
    return;
}

/******* Declarations *******/

#endif // _2D_ALGO_FFT_HPP_


