#ifndef SIGMA_H
#define SIGMA_H

class SIGMAOPTS{

 public:

  // input options
  bool sigma_is_on;    // if we calculate sigma (true/1) or not(false/0)
  bool sigma_static;   // if true/1, do static sigma
  bool sigma_dynamic;  // if true/1, do dynamic sigma. If false, not performing dynamic sigma

  int npt;             // total number of data points that will be evaluated (see sigmadata.in file)
  double **eval;       // n, n', and \omega  (TODO: n and n' should be an integer)

  void readInputFile(char*);
  void readInputData(char*);
  void state_class_out(){
    FILE *fp; fp=fopen("SIGMAOPTS_class.out","w");
    fprintf(fp,"Sigma calculation is on/off (1/0): %d\n",sigma_is_on);
    fprintf(fp,"Static sigma is on/off (1/0): %d\n",sigma_static);
    fprintf(fp,"Dynamic sigma is on/off (1/0): %d\n",sigma_dynamic);
  }
  
  // default constructor
  SIGMAOPTS(){
  
    char fileName[256]
    sprintf(fileName,"sigma.in");

    // read sigma.in file. 
    readInputFile(fileName);

    sprintf(fileName,"sigmadata.in");

    // read sigmadata.in file 
    // this file includes n, n', and \omega values
    readInputData(fileName);

    // write out class variables
    state_class_out();
  };
  
};

#endif
