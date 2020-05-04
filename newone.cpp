#include <iostream>
#include "matrix.h"


int main(){

    std::vector<std::vector<double>> teste;
    std::vector<std::vector<double>> teste2;
    teste.push_back({1, 2, 3});
    //teste.push_back({9, 10});
    teste2.push_back({1});
    teste2.push_back({0});
    teste2.push_back({3});

    Matrix Mtx1(teste);
    Matrix Mtx2(teste2);

    Mtx1.ShowContent();
    Mtx2.ShowContent();


    Matrix Mtx3 = Mtx2+Mtx1;
    Mtx3.ShowContent();

    std::cout<<"ARRUME O ERRO DE MULTIPLICAÇÃO DE MATRIZES"<<std::endl;
return 0;
}