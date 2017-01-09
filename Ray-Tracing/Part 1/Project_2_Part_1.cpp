//
//  main.cpp
//  Project 2
//
//  Created by test on 15/4/2016.
//  Copyright © 2016年 gilbertan. All rights reserved.
//

#include <random>
#include <iostream>
#include <math.h>

double findIbar(int size, double sample[]){
    double sum = 0.0;
    for (int n = 0; n < size; ++n)
            sum += sample[n];
    return sum/double(size);
}

double findSampleSD(double Ibar, double sample[], int size){
    double sum = 0.0;
    for (int n = 0; n < size; ++n) {
        sum += pow((sample[n]-Ibar ),2);
    }
    return pow(sum/double(size-1),0.5);
}

double findSbar(double s, int size){
    return (2*s)/pow(double(size),0.5);
}

double findI(double p, double funct){
    return funct/p;
}

bool ray(double x[], double w[]){
    double centerS[3] = { 1.0,1.0,5.0 }, radius = 1.0, translatedS[3];
    double L1 = centerS[2]/cos(w[0]), L2 = L1*sin(w[0]) ,L3 = L2*sin(w[1]), L4 = L2*cos(w[1]);
    translatedS[0] = 1.0 - L3;
    translatedS[1] = 1.0 - L4;
    translatedS[2] = 0.0;
    double checkIn = pow(pow(x[0]-translatedS[0], 2.0) + pow(x[1]-translatedS[1], 2.0),0.5);
    if ( checkIn <= pow(radius,2.0))
        return true;
    else
        return false;
}

int main(int argc, const char * argv[]) {
    // Task 1-1
    const double p1 = 0.25;
    const int size = 100000;
    int i = 0;
     std::cout << "------------------------Task 1-1------------------------" << std::endl;
    while ( i < 10){
        double sample[size];
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-2, 2);
        for (int n = 0; n < size; ++n)
            sample[n] = findI(p1,exp(-pow(dis(gen),2.0)/2.0));

        double Ibar = findIbar(size, sample);
        double s = findSampleSD(Ibar, sample, size);
        double sbar = findSbar(s, size);
        
        std::cout << "Round " << i+1 << ":" << std::endl;
        std::cout << "Ibar = " << Ibar << std::endl;
        std::cout << "Sbar = "  << sbar <<std::endl;
        std::cout << "Ibar - Sbar = " << Ibar-sbar << std::endl;
        std::cout << "Ibar + Sbar = " << Ibar+sbar << std::endl;
        std::cout << "\n";
        i++;
    }
    
    //Task1-2
    int k = 0;
    const double lambda = 1.0;
    std::cout << "------------------------Task 1-2------------------------" << std::endl;
    while ( k < 10){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);
        double sample[size], tempX, tempP;
        for (int n = 0; n < size; ++n){
            tempX = 1 - (log(dis(gen))/lambda);
            tempP = lambda*exp(-(lambda*(tempX - 1)));
            sample[n] = findI(tempP, exp(-pow(tempX,2.0)/2.0));
        }
        double Ibar = findIbar(size, sample);
        double s = findSampleSD(Ibar, sample, size);
        double sbar = findSbar(s, size);
    
        std::cout << "Round " << k+1 << ":" << std::endl;
        std::cout << "Ibar = " << Ibar << std::endl;
        std::cout << "Sbar = "  << sbar <<std::endl;
        std::cout << "Ibar - Sbar = " << Ibar-sbar << std::endl;
        std::cout << "Ibar + Sbar = " << Ibar+sbar << std::endl;
        std::cout << "\n";
        k++;
    }
    
    //Task2-2
    std::cout << "------------------------Task 2-1------------------------" << std::endl;
    int g = 0;
    const double p21 = 1.0/(M_PI*M_PI);
    while ( g < 10){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);
        double sample[size], tempX;
        for (int n = 0; n < size; ++n){
            tempX = sin((M_PI/2.0)*dis(gen));
            sample[n] = findI(p21, tempX);
        }
        double Ibar = findIbar(size, sample);
        double s = findSampleSD(Ibar, sample, size);
        double sbar = findSbar(s, size);
        
        std::cout << "Round " << g+1 << ":" << std::endl;
        std::cout << "Ibar = " << Ibar << std::endl;
        std::cout << "Sbar = "  << sbar <<std::endl;
        std::cout << "Ibar - Sbar = " << Ibar-sbar << std::endl;
        std::cout << "Ibar + Sbar = " << Ibar+sbar << std::endl;
        std::cout << "\n";
        g++;
    }
    
    std::cout << "------------------------Task 2-2------------------------" << std::endl;
    int j = 0;
    const double p22 = 1.0/(2.0*M_PI);
    while ( j < 10){
        std::random_device rd1;
        std::mt19937 gen1(rd1());
        std::random_device rd2;
        std::mt19937 gen2(rd2());
        std::uniform_real_distribution<> dis(0, 1);
        
        double theta, phi;
        
        double sample[size], funct;
        for (int n = 0; n < size; ++n){
            theta = acos(dis(gen1));
            phi = 2*M_PI*dis(gen2);
            funct = pow(cos(theta)*sin(theta)*cos(phi),2.0)*sin(M_PI/2);
            sample[n] = findI(p22, funct);
        }
        double Ibar = findIbar(size, sample);
        double s = findSampleSD(Ibar, sample, size);
        double sbar = findSbar(s, size);
        
        std::cout << "Round " << j+1 << ":" << std::endl;
        std::cout << "Ibar = " << Ibar << std::endl;
        std::cout << "Sbar = "  << sbar <<std::endl;
        std::cout << "Ibar - Sbar = " << Ibar-sbar << std::endl;
        std::cout << "Ibar + Sbar = " << Ibar+sbar << std::endl;
        std::cout << "\n";
        j++;
    }
    
    std::cout << "------------------------Task 3------------------------" << std::endl;
    int f = 0, size2 = 500000;
    const double p3 = 1.0/(2.0*M_PI);
    while ( f < 10){
        std::random_device rd;
        std::mt19937 gen1(rd());
        std::mt19937 gen2(rd());
        std::mt19937 gen3(rd());
        std::mt19937 gen4(rd());
        std::uniform_real_distribution<> dis1(0, 1);
        std::uniform_real_distribution<> dis2(-0.5, 0.5);

        double x[2], w[2];
        double theta, phi, L = 100.0;
        double sample[size2];
        for (int n = 0; n < size2; ++n){
            theta = acos(dis1(gen1));
            phi = 2*M_PI*dis1(gen2);
            w[0] = theta;
            w[1] = phi;
            x[0] = dis2(gen3);
            x[1] = dis2(gen4);
            double funct, indicaFunct = 1;
            if (ray(x,w))
                indicaFunct = 0;
    
            funct = L*indicaFunct*cos(theta)*sin(M_PI/2);
            sample[n] = findI(p3, funct);
        }
        double Ibar = findIbar(size2, sample);
        double s = findSampleSD(Ibar, sample, size2);
        double sbar = findSbar(s, size2);
        
        std::cout << "Round " << f+1 << ":" << std::endl;
        std::cout << "Ibar = " << Ibar << std::endl;
        std::cout << "Sbar = "  << sbar <<std::endl;
        std::cout << "Ibar - Sbar = " << Ibar-sbar << std::endl;
        std::cout << "Ibar + Sbar = " << Ibar+sbar << std::endl;
        std::cout << "\n";
        f++;
    }
    return 0;
}


