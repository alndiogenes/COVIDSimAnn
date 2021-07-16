#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <conio.h>
#include "math.h"
#include "randomc.h"
#include <time.h>
#include <string.h>


using namespace std;

typedef vector<float> paramVec;
typedef vector<paramVec> pImage;

int readCorr( const char* fileName, paramVec &cOriginal );
int calcMedia( int i, paramVec &cOriginal, paramVec &media );
int probability( double diff, double T, TRandomMersenne &rg );

int main(int argc,char **argv) {
    
    // programa para tentar ajustar uma curva de dist de tam de poro 3D com uma log-normal
    
    if ( argc<2 ) {
	    printf("Uso: curveFitting arqEntrada\n" );
	    printf("Ex: curveFitting dOriginal.txt\n" );
	    return 0;
	}
	
	TRandomMersenne rg ( time (NULL) );
	
	// Curva Log-Normal
	
	paramVec dOriginal;
	paramVec dGerada3D;
	paramVec mediaMovel;
	
	readCorr( argv[1], dOriginal );
	
	cout << "aqui" << endl;
	for (int i=0; i<dOriginal.size(); i++) {
		calcMedia(i, dOriginal, mediaMovel);
	}
	
	
	// Log Logistic CDF
	
	float m=1;
	
	if (argc == 3 ) {
		m = atof(argv[2]);
	}
	
	double alpha = 25;
	double beta = 4;
	double mult = m*dOriginal[dOriginal.size()-1];
	double div = 1.5;
	
	double pAlpha = 0.5;
	double pBeta = 0.05;
	double pMult = mult/100;
	double pDiv = 0.5;
	
	
	ofstream os ("dGerada3D.txt");
	
	dGerada3D.push_back( 0 );
	
	os << "Raio \t dGerada3D \t dOriginal" << endl;
	
	for ( int i=1; i<dOriginal.size(); i++ ) {
		double val = mult/pow(1+pow(i/(alpha*div), -beta),2.0);
		dGerada3D.push_back( val );
		os << i << " \t" << dGerada3D[i]  << " \t" << dOriginal[i]  << endl;
	}
	
	os.close();
	
	//return 0;
	
	double E=0;
	double T=10;
	double decrement = 0.99;
	int it=0;
	int dec=0;

	ofstream ts ("graficoSaida.txt");
	
	ts << "Energia\talpha\tbeta\tmult" << endl;
	
	for ( int i=1; i<dOriginal.size(); i++ ) {
		E += pow( dOriginal[i] - dGerada3D[i], 2 );
	}
	
	printf("E inicial %.5e\n", E );
    
	
	int var=1;
	int varT=1;
	
	int cOut=0;
	
	double bsfE = E;
	double bsfMult = mult;
	double bsfAlpha = alpha;
	double bsfBeta = beta;
	double bsfDiv = div;
	
	while ( ( E>1e-10 )&&( T>0.0000001 ) ) {
		double alphaA = alpha;
		double multA = mult;
		double betaA = beta;
		double divA = div;

	  	if ( kbhit() ) {
		   	printf("Parar a execucao? Pressione \'s\' para parar: ");
		   	char tx[80];
			gets(tx);
		   	if ( strcmp( tx, "s" )==0 ) {
			    break;
			}
		}
		
		double EA = 0;
		
		if ( varT%1000==0 ) {
		    dec=1;
		    varT=1;
		    ts << E << "\t" << alpha*div << "\t" << beta << "\t" << mult << endl;
		}
		
		if ( bsfE > E ) {
		    bsfE = E;
		    varT = 1;
		    cOut=0;
			bsfMult = mult;
			bsfAlpha = alpha;
			bsfBeta = beta;
			bsfDiv = div;
		} else {
		  	varT++;
		  	cOut++;
		}
		
		if ( cOut>100000 ) {
		    break;
		}
		
	    switch( var ) {
			//Log-Logistic
			//var = 4;
			case 1:
				if ( ( rg.Random()>0.5 )&&( beta>=2 ) ) {
				    beta -= pBeta;
				} else {
				  	if ( beta < 2 ) {
				  	    beta += pBeta;
					}
				}
				//var++;
				break;
			case 2:
				if ( ( rg.Random()>0.5 )&&( alpha>5) ) {
				    alpha -= pAlpha;
				} else {
				  	alpha += pAlpha;
				}
				//var++;
				break;
			case 3:
				if ( ( rg.Random()>0.5 )&&( div>0) ) {
				    div -= pDiv;
				} else {
				  	div += pDiv;
				}
				//var++;
				break;
			case 4:
				if ( ( rg.Random()>0.5 )&&( mult>2000 ) ) {
				    mult -= pMult;
				} else {
				  	mult += pMult;
				}
				//var = 1;
				break;
				
		}
		
		var++;
		if (var==5) {
			var = 1;
		}
		
		for ( int i=dOriginal.size()-1; i>=1; i-- ) {
			double val = mult/pow(1+pow(i/(alpha*div), -beta),2.0);
			dGerada3D[i] = val;
		}
		
		for ( int i=dOriginal.size()-1; i>=1; i-- ) {
			EA += pow( dOriginal[i] - dGerada3D[i], 2 );
		}
			
	    if ( it==0 ) {
		    if ( EA>E ) {
			    T = -(EA-E)/log(0.1);
			} else {
			  	if ( EA-E==0 ) {
				    T = 0.1;
				} else {
				  	T = -(E-EA)/log(0.1);
				}
			}
			//T=10;
			printf("T inicial %.5lf\n", T );
		} else {
		  	if ( dec==1 ) {
				T = decrement*T;    // geometrica
				dec = 0;
			}
		}
		if ( E<EA ) {
			if ( probability( EA-E, T, rg )==0 ) {
			    alpha = alphaA;
			    mult = multA;
			    beta = betaA;
			    div = divA;
			    
			} else {
			  	E = EA;
			}
		} else {
		  	E = EA;
		  	//varT = 1;
		}
		it++;
		printf("it %.5d E %.1e bsfE %.1e T %.1e varT %.3d mult %.2e alpha %.2f beta %.2f div %.1f \r", it, E, bsfE, T, varT, mult, alpha, beta, div );
	}
	
	mult = bsfMult;
	alpha = bsfAlpha;
	beta = bsfBeta;
	div = bsfDiv;

	
	cout << endl << "Parametros bsfE " << bsfE << " mult " << mult << " alpha " << alpha << " beta " << beta << " div " << div << endl; 
	
	ts << E << "\t" << alpha*div << "\t" << beta << "\t" << mult << endl;
	ts.close();

	
	printf("\nPrograma Encerrado\n");
	
	ofstream cs ("curvaSaida.txt");
	
	cs << "Raio \tdGerada3D \tOriginal3D \t mult \t" << mult << endl;
	
	for ( int i=1; i<dOriginal.size(); i++ ) {
		double val = mult/pow(1+pow(i/(alpha*div), -beta),2.0);
		cs << i << " \t" << val  << " \t" << dOriginal[i]  << endl;
	}
	cs.close();
    
    return 1;
}


int readCorr( const char* fileName, paramVec &cOriginal ) {
    
    ifstream is( fileName );
    if ( !is.is_open() ) {
	    return 0;
	}
    float c;
	float lixo;
	is >> lixo;
	cOriginal.push_back( 0 );
	while ( is >> c ) {
		cOriginal.push_back( c );
		cout << c << endl;
		is >> lixo;
		//++i;
	}
	is.close();
	return 1;
}

int probability( double diff, double T, TRandomMersenne &rg ) {
	float p =  rg.Random();

	if ( exp( -diff/T ) > p ) {
	  	 return 1;
    } else {
	  	 return 0;
    }
    return 1;
}

int calcMedia( int i, paramVec &cOriginal, paramVec &media ) {
	
	double medTemp=0;
	int k;
	double size;
	if (i<6) {
		k=-1;
		size=i+1;
	} else {
		k=i-7;
		size=7;
	}
	//cout << i << " " << k << " " << size << endl;
	for (int j=i; j>k; j--) {
		medTemp += cOriginal[j];
		//cout << " " << medTemp << " " << cOriginal[j] << endl;
	}
	medTemp = medTemp/size;
	media.push_back(medTemp);
	//cout << " " << medTemp << endl;
	return 1;	
}





