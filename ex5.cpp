#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>

double f(double x) {
    return exp(-x) * pow(x, 3);
}

// Метод Неймана (метод браковки)
double neymanMethod(int n) {
    double a = 0;
    double b = 1;
    double h = (b - a) / n;
    double integral = 0.0;

    for (int i = 0; i < n; i++) {
        double x1 = a + i * h;
        double x2 = a + (i + 1) * h;
        integral += (f(x1) + f(x2)) * h / 2;
    }
    
    return integral;
}

// Метод выделения главной части
double extractionMethod(int n) {
    double a = 0;
    double b = 1;
    double integral = 0.0;
    double mainPartIntegral = 0.0;
    for (int i = 0; i < n; i++) {
        double x = a + (b - a) * i / n;
        mainPartIntegral += f(x);
    }
    mainPartIntegral *= (b - a) / n;
    for (int i = 0; i < n; i++) {
        double x = a + (b - a) * i / n;
        integral += (f(x) - mainPartIntegral);
    }

    return mainPartIntegral + integral * (b - a) / n;
}

// Метод Монте-Карло
double monteCarloMethod(int n) {
    TRandom3 randGen(0); 
    double a = 0;
    double b = 1;
    double integral = 0.0;

    for (int i = 0; i < n; i++) {
        double x = randGen.Uniform(a, b);
        integral += f(x);
    }

    return integral * (b - a) / n;
}

void calculateIntegrals() {
    int n = 10000; 


    double integralNeyman = neymanMethod(n);
    double integralExtraction = extractionMethod(n);
    double integralMonteCarlo = monteCarloMethod(n);

    std::cout << "Integral using Neyman Method: " << integralNeyman << std::endl;
    std::cout << "Integral using Extraction Method: " << integralExtraction << std::endl;
    std::cout << "Integral using Monte Carlo Method: " << integralMonteCarlo << std::endl;
}