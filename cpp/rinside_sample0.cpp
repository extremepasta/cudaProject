#include <RInside.h>         

extern "C"
int quantMain(double ** matrix, int nR, int nC);

extern "C"
int secondMain(double ** matrix, int nR, int nC, double tau);

           // for the embedded R via RInside

Rcpp::NumericMatrix createMatrix(const int n) {
    Rcpp::NumericMatrix M(n,n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            M(i,j) = i*10 + j; 
        }
    }
    return(M);
}

int main(int argc, char *argv[]) {

    RInside R(argc, argv);              // create an embedded R instance 

    std::string txt = "x <- matrix(rnorm(10), nrow=5, ncol=5)";
    Rcpp::NumericMatrix xx = R.parseEval(txt);

    double** mat = new double*[5];

    for(int i = 0; i < 5; ++i)
        mat[i] = new double[5];

    for(int r = 0; r<5; r++)
        for(int c = 0; c<5; c++)
            mat[r][c] = (double)xx(r,c);

    //ACQUISIZIONE DATI
    std::cout << "Insert your matrix: " << std::endl;
    //R["x"] = createMatrixFromData(R);
    
    R["x"] = xx;
    //CALCOLO QUANTILE
    //R["tau"] = 
    double tau = (double)R.parseEval("library(ldfun);"
        "tau <- quantile.localdepth.functional(x, probs=0.3);");
    
    std::cout << "Qua " << std::endl;
    //R["ld"] = 

    int x = quantMain(mat, 5,5);

    std::cout << "Qui " << std::endl;

    int j = secondMain(mat, 5, 5, tau);
    
    exit(0);
}




