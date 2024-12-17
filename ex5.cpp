#include <TH1D.h>
#include <ROOT/TSeq.hxx>
#include <iostream>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>

//График
void Hist(TRandom3 *randomGenerator)
{
    // Создание холста для отображения гистограммы
    TCanvas *canvas = new TCanvas("c1", "Histogram for one formula", 800, 600);
    TH1D *histogram = new TH1D("histogram", "Distribution", 100, 0, 2.5);

    // Количество выборок
    int totalSamples = 100000;

    // Заполнение гистограммы значениями
    for (auto i : ROOT::TSeqI(totalSamples))
        histogram->Fill(TMath::Sqrt(-TMath::Log(randomGenerator->Rndm())));

    // Отрисовка гистограммы 
    histogram->Draw();
}

// Функция для вычисления значения функции
double CalculateFunction(double x) 
{ 
    return TMath::Exp(-x) * x * x * x; 
}

// Метод выбраковки Неймана для интегрирования
double RejectionSampling(TRandom3 *randomGenerator, const unsigned int totalSamples)
{
    double lowerBound = 0, upperBound = 1, constant = 0.9; // Параметры интегрирования
    double acceptedSamples = 0; // Счетчик принятых образцов

    // Генерация 
    for (auto i : ROOT::TSeqI(totalSamples))
    {
        double randomX = (upperBound - lowerBound) * randomGenerator->Rndm() + lowerBound; 
        double randomY = constant * randomGenerator->Rndm();
        if (randomY <= CalculateFunction(randomX))
            acceptedSamples++;
    }
    

    return acceptedSamples / totalSamples * (upperBound - lowerBound) * constant;
}

// Метод среднего для интегрирования
double AverageSampling(TRandom3 *randomGenerator, const unsigned int totalSamples)
{
    double lowerBound = 0, upperBound = 1; // Параметры интегрирования
    double sum = 0; // Сумма значений функции

    // Генерация образцов методом среднего
    for (auto i : ROOT::TSeqI(totalSamples))
        sum += (upperBound - lowerBound) * CalculateFunction((upperBound - lowerBound) * randomGenerator->Rndm() + lowerBound);

    // Возвращаем среднее значение интеграла
    return sum / totalSamples;
}

// Метод важной выборки для интегрирования
double ImportanceSampling(TRandom3 *randomGenerator, const unsigned int totalSamples)
{
    double lowerBound = 0, upperBound = 1; // Параметры интегрирования
    double sum = 0; // Сумма значений функции

    // Генерация образцов методом важной выборки
    for (auto i : ROOT::TSeqI(totalSamples))
    {
        double randomSample = TMath::Sqrt(TMath::Sqrt(randomGenerator->Rndm()));
        sum += CalculateFunction(randomSample) / (4 * randomSample * randomSample * randomSample);
    }
    
    // Возвращаем среднее значение интеграла
    return sum / totalSamples;
}

// Основная функция для выполнения задачи
void task5(int sampleCount = 1000, int histogramCount = 1000)
{
    TRandom3 *randomGen = new TRandom3(); // Создаем генератор случайных чисел
    randomGen->SetSeed(time(NULL)); // Устанавливаем начальное значение генератора

    Hist(randomGen); // Генерируем и отображаем гистограмму

    double exactValue = 6 - 16 / TMath::E(); // Точное значение интеграла
    double rejectionResult   = RejectionSampling(randomGen, sampleCount); // Результат метода выбраковки
    double importanceResult   = ImportanceSampling(randomGen, sampleCount); // Результат метода важной выборки
    double averageResult      = AverageSampling(randomGen, sampleCount); // Результат метода среднего

    std::cout << std::fixed;
    std::cout.precision(7);

    std::cout << "Exactly:  " << exactValue << "\t  " << std::endl;
    std::cout << "Rejection:  " << rejectionResult << "\t  "  << std::endl;
    std::cout << "Average:  " << averageResult << "\t  " << std::endl;
    std::cout << "Importance:  " << importanceResult << "\t  " << std::endl;
}