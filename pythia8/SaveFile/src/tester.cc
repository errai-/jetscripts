#include <iostream>
#include <fstream>
#include <limits>
#include <string>

#include "MinimalEvent.h"

  using std::cout;
int main(void){
  std::ofstream output;
  output.open("tester.txt", std::ofstream::out | std::ofstream::app);

  output << "sis on sas\n";
  output << "\n";
  output << "\n";
  output.close();

  std::ifstream input;
  input.open("pythia8data.txt", std::ifstream::in);
  std::string jaa;
  MinimalEvent ses;
  for (size_t i=0; i!=100; ++i){
    ses.Read(&input);
  }
  std::cout << ses.particles << std::endl;
  for (size_t i=0; i!=ses.particles; ++i){
    std::cout << ses.px[i] << " " << ses.py[i] << " " << ses.pz[i] << " " << ses.e[i] << " " << ses.id[i] << std::endl;
  }

  std::cout << std::numeric_limits<double>::digits << " " <<  std::numeric_limits<double>::digits10 << std::endl;
  std::cout << std::numeric_limits<uint>::digits << " " <<  std::numeric_limits<uint>::digits10 << std::endl;
  std::cout << std::numeric_limits<int>::digits << " " <<  std::numeric_limits<int>::digits10 << std::endl;
  std::cout << std::numeric_limits<size_t>::digits << " " <<  std::numeric_limits<size_t>::digits10 << std::endl;
  
  return 0;
}
