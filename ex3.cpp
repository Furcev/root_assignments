#include <iostream>
#include <cmath>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMinuit.h>

//Для построения графика используйте run_variational_method()
//для запуска алгоритма минимизации find_minimum()

double potential(double x) {

if (x<-10 && x>10){
    return 0 ;
}

if (x>-10*pow(10, -10)){
    return -0.5;
}

}

double wave_function(double x, double a) {
    return pow((2.0/3.14), 0.25)*(1/ sqrt(a)) * exp(-x * x / (a * a));
}
double energy(double a) {
    const int n_points = 1000; 
    const double x_min = -10; 
    const double x_max = 10; 
    const double dx = (x_max - x_min) / n_points;

    double integral = 0.0;
    for (int i = 0; i < n_points; ++i) {
        double x = x_min + i * dx;
        double psi = wave_function(x, a);
        integral += psi*psi*(-3.8*(-2*a*a+4*x*x)/(a*a*a*a)+potential(x)) * dx;
    }
    
    return integral;
}

void energy_function(int &npar, double *gin, double &f, double *par, int iflag) {
    f = energy(par[0]);
}
void find_minimum() {
    TMinuit minuit(1);
    
    minuit.SetFCN(energy_function);

    double arglist[10];
    Double_t ierflg = 0;
    Double_t initial_a = 1.0; 
    minuit.DefineParameter(0, "a", initial_a, 0.1, 0, 0);

    arglist[0] = 1; 
    minuit.Command("MIGRAD");

    double min_a, min_energy;
    minuit.GetParameter(0, min_a, ierflg);
    min_energy = energy(min_a);

    std::cout << "Минимальное значение a: " << min_a << std::endl;
    std::cout << "Минимальное значение энергии: " << min_energy << std::endl;
}


void run_variational_method() {

    const int n_graph_points = 1000;
    TGraph *graph = new TGraph(n_graph_points);
    double a_value = 7;
    for (int i = 0; i < n_graph_points; ++i) {
        a_value = a_value + 0.01 ;
        double energy_value = energy(a_value); 
        graph->SetPoint(i, a_value, energy_value);
    }
    printf("%f",energy(0.02e-12) );
    graph->SetTitle("E(a)");
    
    TCanvas *canvas = new TCanvas("canvas", "Energy vs Parameter a", 800, 600);
    graph->Draw("AL"); 
    canvas->Update();
}

int main() {
    run_variational_method();
    return 0;
}