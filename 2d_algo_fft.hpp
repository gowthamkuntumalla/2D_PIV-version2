/** This Header file consists a total of 11 Functions **/
#ifndef _2D_ALGO_FFT_HPP_
#define _2D_ALGO_FFT_HPP_

#define PI 3.14159265
using namespace std;
using namespace cv;


/**access in vector is [row][column]**/
/** <i><x><col> are equivalent and similarly <j><y><row> are equivalent**/
/** convention : (col,row), (x,y) order is followed**/
/** cartesian coordinates (c,r) i.e., the top left corner represents the address of subwindow**/

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
void max_coef(Mat t,const int win_size,int& max_x,int& max_y )
{
    double a=-1;// all intensities are greater than this
    for(int i=0; i<win_size; i++)
    {
        for(int j=0; j<win_size; j++)
        {
            double mag= t.at<uchar>(j,i);
            if(mag>a)//access by t(row,column)
            {
                a=mag;//(row,column)
                max_x=i;
                max_y=j;
            }
        }
    }
    return;
}
/***********************Main Function************************/
void piv_2d_fft(Mat image1, Mat image2,vector< vector <pair<int,int> > > &max_coef_point,int win_size,int step_size,int i) // i= image number
{
    ofstream myfile; // like 'cout' it outputs to a file
    double sd1=0,sd2=0;
    double avg1,avg2;

    /**This is valid for displacements of maximum A/3 where A is window_Size **/
    int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image1.
    Mat win1=Mat::zeros(win_size,win_size,CV_8U);
    Mat win2=Mat::zeros(win_size,win_size,CV_8U);//interrogation area
    Mat fft_win1,fft_win2;//Fourier transformed interrogation areas
    Mat fft_mat_prod;//multiplication of i1,i2'
    Mat mat_prod;
    myfile.open ("data_"+to_string(i)+".txt");// output to a text file

    /** Loop over entire image**/
    for(int c=0; c<=(totcol1-win_size); c+=step_size)
    {
        for(int r=0; r<=(totrow1-win_size); r+=step_size)
        {
            int m=0,n=0;//displacement of the max coefficient point
            myfile<<(c+win_size/2)<<","<<(r+win_size/2)<<", ";//initial center point (x,y) of the interrogation spot
            sd1=sd(win_size+c,win_size+r,image1,c,r);//standard deviation
            sd2=sd(win_size+c,win_size+r,image2,c,r);
            avg1=avg(win_size+c,win_size+r,image1,c,r);//averages
            avg2=avg(win_size+c,win_size+r,image2,c,r);
            for(int ra=0; ra!=win_size; ra++)
            {
                for(int ca=0; ca!=win_size; ca++)
                {
                    //copy elements
                    win1.at<uchar>(ra,ca)=((int)image1.at<uchar>(r+ra,c+ca)-avg1)/sd1;
                    win2.at<uchar>(ra,ca)=((int)image2.at<uchar>(r+ra,c+ca)-avg2)/sd2;
                }
            }
            Mat planes1[] = {Mat_<float>(win1), Mat::zeros(win1.size(), CV_32F)};
            Mat planes2[] = {Mat_<float>(win2), Mat::zeros(win2.size(), CV_32F)};
            merge(planes1, 2,fft_win1 );
            merge(planes2, 2,fft_win2 );
            /**do FFT for image1 & image2**/
            dft(fft_win1,fft_win1);
            dft(fft_win2,fft_win2);

            /**calculate i1*(i2')**/
            mulSpectrums(fft_win1, fft_win2, fft_mat_prod,0,true);

            /**do inverse FFT**/
            idft(fft_mat_prod,mat_prod);
            /** find maximum correlation coefficient point**/
            split(mat_prod, planes1);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
            magnitude(planes1[0], planes1[1], planes1[0]);// planes[0] = magnitude
            Mat magI = planes1[0];
            max_coef(magI,win_size,m,n);

            /**Save the data in a "data_i.txt" file**/
            max_coef_point[r][c].first=n+r;//row index
            max_coef_point[r][c].second=m+c;//column index
            myfile<<c+m<<","<<n+r<<endl;//cartesian displacements (m,n)
            //myfile << "Writing this to a file.\n";
        }
    }
    myfile.close();
    return;
}
#endif // _2D_ALGO_FFT_HPP_
