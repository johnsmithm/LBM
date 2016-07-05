#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>
#include <string>



#include "LBM.h"
#include "Timer.h"


using namespace std;




void run_test(){
	

}

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>0);
	string par = argv[1];

	vi param;
	// sizeXCylinder, sizeCylinderY, ballCenterX, ballCenterY, ballDiameter
	vi cylinder({0.06,0.02,0.02,0.008,0.005});
	int nr=0;
	if(par == "scenario1"){
		param = getParameters({0.000001,0.01,30},cylinder);
		nr = 3./param[0];
	}else{
		param = getParameters({0.000001,0.0016,60},cylinder);
		nr = 5./param[0];
	}
	LatticeB simulator(param);
	//nr = 500;
	cerr<<"nr:"<<nr<<"\n";
	simulator.run(nr);
	simulator.getResult();
}
