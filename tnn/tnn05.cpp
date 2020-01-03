/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/////////////////////// Test 01 fitting quadratic polynomial

//multivariable regression

// MechSys
#include <mechsys/nn/Domain.h>
#define PI 3.1415926

int main(int argc, char **argv) try
{
    
 size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    size_t ni = 2;  //Number of inputs
    size_t no = 2;  //Number of outputs
    size_t nh = 1;  //Number of hidden layers
    size_t nn = 200;//Number of neurons per hidden layer
    size_t nt = 220; //Number of training sets

    //Preparing input and output data for training sets
    Array<Array <double> > In;
    Array<Array <double> > Out;
    In .Resize(nt);
    Out.Resize(nt);
    
    //Build a random input
    std::default_random_engine e;
    std::uniform_real_distribution<float> u(0.0, 1.0);

    double dx1 = u(e);
    double dx2 = u(e);

    for (size_t i=0;i<nt;i++) //For each training set
    {
        In [i].Resize(ni);
        Out[i].Resize(no);


        In [i][0] = dx1;
        In [i][1] = dx2;

        Out[i][0] = dx1*dx1 + dx2*dx2;
        Out[i][1] = sin(PI*(dx1*dx1 + dx2*dx2)/2);

        dx1 = u(e);
        dx2 = u(e);
    }

    NN::Domain dom(In,Out,nh,nn);
    dom.Initialize();
    dom.Alpha = 0.1; //learning rate
    dom.Type  = 5; //choose the activation function 
    //dom.Train(40000/*epochs*/,Nproc/*number of cores*/);

    //modification1
    dom.Train(40000/*epochs*/,Nproc/*number of cores*/);
    
    //Try the trained NN
    In[0][0] = 0.4;
    In[0][1] = 0.4;

    dom.Predict(In[0],Out[0],1);

    //The caculation of error
    std::cout<< "---------- The results ----------"<<std::endl;
    std::cout << "The output1: " << Out[0][0] << std::endl;
    std::cout << "The error1: " << fabs((Out[0][0] - (In[0][0]*In[0][0] + In[0][1]*In[0][1]))/(In[0][0]*In[0][0] + In[0][1]*In[0][1])) << std::endl;


    std::cout << "The output1: " << Out[0][1] << std::endl;
    std::cout << "The error2: " << fabs((Out[0][1] - (sin(PI*(In[0][0]*In[0][0] + In[0][1]*In[0][1])/2)))/(sin(PI*(In[0][0]*In[0][0] + In[0][1]*In[0][1])/2))) << std::endl;


    //Saving the NN
    dom.Save("tnn05");
    return 0;

}
MECHSYS_CATCH
