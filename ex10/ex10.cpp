#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>

// Определяем индексы параметров для удобства
enum { BackgroundParam, SignalScaleParam, SignalMeanParam, SignalSigmaParam };

// Объявление функций для расчета фона, сигнала и функции подгонки
double CalculateBackground(double *x, double *params); 
double CalculateSignal(double *x, double *params);
double CalculateFitFunction(double *x, double *params);

// Функция для вычисления значения хи-квадрат
void ComputeChiSquared(int &npar, double *gin, double &f, double *params, int iflag);

// Функция для заполнения гистограммы данными из файла
void PopulateHistogram(TH1D *histogram, TString filePath);


void Analysis() {
    // Создаем две гистограммы для хранения данных
    TH1D *firstHistogram = new TH1D("firstHistogram", "data_1.dat", 100, 500, 600);
    TH1D *secondHistogram = new TH1D("secondHistogram", "data_2.dat", 100, 500, 600);

    // Заполняем гистограммы данными из файлов
    PopulateHistogram(firstHistogram, "data_1.dat");
    PopulateHistogram(secondHistogram, "data_2.dat");

    // Инициализируем массив параметров для подгонки
    double fitParams[4] = {0}; 
    fitParams[BackgroundParam] = 1.0; // Начальное значение фона
    fitParams[SignalMeanParam] = 550;  // Начальное значение среднего сигнала
    fitParams[SignalSigmaParam] = 10;   // Начальное значение сигма сигнала
    fitParams[SignalScaleParam] = 100;  // Начальное значение масштаба сигнала

    // Создаем функцию подгонки с заданными параметрами
    TF1 *fittingFunction = new TF1("fittingFunction", CalculateFitFunction, 500, 600, 4);
    fittingFunction->SetParameters(fitParams);

    // Выполняем подгонку первой гистограммы с использованием функции подгонки
    firstHistogram->Fit(fittingFunction, "R");
    
    // Создаем холст для отображения результатов
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);
    firstHistogram->SetLineColor(kBlue);   
    secondHistogram->SetLineColor(kRed);  
    
    // Рисуем гистограммы и функцию подгонки на одном холсте
    firstHistogram->Draw();
    secondHistogram->Draw("SAME");
    fittingFunction->Draw("SAME");
    
    // Создаем легенду для обозначения графиков
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(firstHistogram, "Data 1", "l");
    legend->AddEntry(secondHistogram, "Data 2", "l");
    legend->AddEntry(fittingFunction, "Fitting Function", "l");
    legend->Draw();
    
    // Получаем значение хи-квадрат из функции подгонки
    double chiSquareValue = fittingFunction->GetChisquare();
    
    // Записываем результаты в текстовый файл
    std::ofstream outputFile("fitting_results.txt");
    outputFile << "Chi-squared: " << chiSquareValue << std::endl;
    outputFile << "Parameters:" << std::endl;
    outputFile << "Background: " << fitParams[BackgroundParam] << std::endl;
    outputFile << "Signal Scale: " << fitParams[SignalScaleParam] << std::endl;
    outputFile << "Signal Mean: " << fitParams[SignalMeanParam] << std::endl;
    outputFile << "Signal Sigma: " << fitParams[SignalSigmaParam] << std::endl;

    outputFile.close(); 
}


double CalculateBackground(double *x, double *params) {
    return params[BackgroundParam]; 
}


double CalculateSignal(double *x, double *params) {
    return params[SignalScaleParam] * TMath::Gaus(x[0], params[SignalMeanParam], params[SignalSigmaParam], true);
}


double CalculateFitFunction(double *x, double *params) {
    return CalculateBackground(x, params) + CalculateSignal(x, params);
}


void ComputeChiSquared(int &npar, double *gin, double &f, double *params, int iflag) {
    TH1D *firstHistogram = nullptr; 
    TH1D *secondHistogram = nullptr;


    gROOT->GetObject("firstHistogram", firstHistogram);
    gROOT->GetObject("secondHistogram", secondHistogram);
    
    const int numBinsFirst = firstHistogram->GetNbinsX(); // Количество бинов в первой гистограмме
    const int numBinsSecond = secondHistogram->GetNbinsX(); // Количество бинов во второй гистограмме
    
    double chiSquared = 0; // Инициализация значения хи-квадрат
    double difference, count, xValue;

    // Цикл по всем бинам первой гистограммы для вычисления хи-квадрат
    for (int i = 1; i <= numBinsFirst; ++i) {
        count = firstHistogram->GetBinContent(i); // Получаем количество событий в текущем бине
        if (count != 0) { // Игнорируем пустые бины
            xValue = firstHistogram->GetBinCenter(i); // Получаем центр бина
            difference = (count - CalculateFitFunction(&xValue, params)) / TMath::Sqrt(count); // Вычисляем разницу
            chiSquared += difference * difference; // Обновляем значение хи-квадрат
        }
    }

    // Цикл по всем бинам второй гистограммы для вычисления хи-квадрат
    for (int i = 1; i <= numBinsSecond; ++i) {
        count = secondHistogram->GetBinContent(i); 
        if (count != 0) { 
            xValue = secondHistogram->GetBinCenter(i); 
            difference = (count - CalculateBackground(&xValue, params)) / TMath::Sqrt(count); 
            chiSquared += difference * difference; 
        }
    }
    f = chiSquared;
}
void PopulateHistogram(TH1D *histogram, TString filePath) {
    std::ifstream inputFile(filePath.Data(), std::ios::in); 
    if (inputFile.is_open()) {
        double value;
        while (inputFile >> value) {
            histogram->Fill(value);
        }
    } 
}