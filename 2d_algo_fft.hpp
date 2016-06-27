/** This Header file consists of Functions **/
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
            double a = image.at<double> (j,i);//(row,column)
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
            double a = image.at<double> (j,i);//(row,column)
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
            double mag= t.at<double>(j,i);
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
void piv_2d_fft(Mat image1, Mat image2,vector< vector <pair<double,double> > > &max_coef_point,int win_size,int step_size,int i) // i= image number
{
    ofstream myfile; // like 'cout' it outputs to a file
    double sd1=0,sd2=0;
    double avg1=0,avg2=0;

    /**This is valid for displacements of maximum A/3 where A is window_Size **/
    int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image1.
    Mat win1=Mat::zeros(win_size,win_size,CV_64F);
    Mat win2=Mat::zeros(win_size,win_size,CV_64F);//interrogation area
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
                    win1.at<double>(ra,ca)=(image1.at<double>(r+ra,c+ca)-avg1)/sd1;
                    win2.at<double>(ra,ca)=(image2.at<double>(r+ra,c+ca)-avg2)/sd2;
                }
            }
            Mat ab1,ab2;
            Mat planes1[] = {Mat_<double>(win1), Mat::zeros(win1.size(), CV_64F)};
            Mat planes2[] = {Mat_<double>(win2), Mat::zeros(win2.size(), CV_64F)};
            merge(planes1, 2,ab1 );
            merge(planes2, 2,ab2 );

            /**do FFT for image1 & image2**/
            dft(ab1,fft_win1);
            dft(ab2,fft_win2);

            /**calculate i1*(i2')**/
            mulSpectrums(fft_win1, fft_win2, fft_mat_prod,0,true);

            /**do inverse FFT**/
            dft(fft_mat_prod,mat_prod,cv::DFT_INVERSE);

            /**Patching up**/
            split(mat_prod, planes2);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
            magnitude(planes2[0], planes2[1], planes2[0]);// planes[0] = magnitude. Planes(number) can be either 1 or 2
            Mat magI_nopad = planes2[0];// no padding here
            Mat magI;
            double cx = magI_nopad.cols/2;
            double cy = magI_nopad.rows/2;

            Mat q0(magI_nopad, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
            Mat q1(magI_nopad, Rect(cx, 0, cx, cy));  // Top-Right
            Mat q2(magI_nopad, Rect(0, cy, cx, cy));  // Bottom-Left
            Mat q3(magI_nopad, Rect(cx, cy, cx, cy)); // Bottom-Right

            Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
            q0.copyTo(tmp);
            q3.copyTo(q0);
            tmp.copyTo(q3);

            q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
            q2.copyTo(q1);
            tmp.copyTo(q2);
            // normalize(magI_nopad, magI_nopad, 0, 1, CV_MINMAX);
            /** Padding **/
            copyMakeBorder(magI_nopad, magI, 5, 5, 5, 5, cv::BORDER_CONSTANT, Scalar::all(0));// 5 pixels extra on all four sides
            //magI=magI_nopad;
            /** Find maximum correlation coefficient point**/
            max_coef(magI,win_size+10,m,n);

            /** Three point peak estimation **/
            //Peak Centroid Method
            double Xo_num   = (m-1)*magI.at<double>(n,(m-1)) + (m)*magI.at<double>(n,m) + (m+1)*magI.at<double>(n,(m+1));//in matrix (row,col) is to be used ::(y,x)
            double Xo_denum = magI.at<double>(n,(m-1)) + magI.at<double>(n,m) + magI.at<double>(n,(m+1));
            double Yo_num   = (n-1)*magI.at<double>((n-1),m)+(n)*magI.at<double>(n,m) + (n+1)*magI.at<double>((n+1),m);
            double Yo_denum = magI.at<double>((n-1),m)+ magI.at<double>(n,m) + magI.at<double>((n+1),m);

            double Xo = Xo_num/Xo_denum;
            double Yo = Yo_num/Yo_denum;
            Xo=Xo-5;// shifting back because of padding
            Yo=Yo-5;

            /**Save the data in a "data_i.txt" file**/
            max_coef_point[r][c].first=Yo+r;//row index/**Yo**/
            max_coef_point[r][c].second=Xo+c;//column index/**Xo**/
            myfile<<Xo+c<<","<<Yo+r<<endl;//cartesian displacements (Xo,Yo)
            //myfile << "Writing this to a file.\n";
        }
    }
    myfile.close();
    return;
}
#endif // _2D_ALGO_FFT_HPP_
