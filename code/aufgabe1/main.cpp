#include <iostream>
#include <cstdlib>
#include<stdlib.h>
#include<time.h>
//#include <sparsematrix.cpp>
#include <eigen-3.3.7/Eigen/Dense>


int main() {

    for(int i = 1; i <= 10; i++)
        std::cout << rand() %10 << std::endl;
    return 0;
}
