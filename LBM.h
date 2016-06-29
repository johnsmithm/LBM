#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <fstream>
#include <list>


using namespace std;



typedef vector<vector<vector<double> > > ma;
typedef vector<double > vi;

#define x first
#define y second
#define pb push_back
#define mp make_pair
#define forn(i,n) for(size_t i=0;i<n;++i)

/**
 * get the lattice parameter from real parameters
 *
 * @param given {viscosity,accceleration,resolution of cylinder diameter}
 * @param cylinder{sizeX, sizeY, ballCenterX, ballCenterY, ballDiameter}
 * return {delta t, delta x, acceleration, number cell X, number cell Y, ballCenterCellX, ballCenterCellY, ballDiameter}
 */
vi getParameters(vi given, vi cylinder){
	
	double dt, dx, w, acc, vis;
	dx = cylinder[1] / given[2];
	dt = 0.0000001;
	vis = given[0] * (dt/(dx*dx));
	w = 1./(3.*vis+0.5);
	cerr<<"w0:"<<w<<" vis0:"<<vis<<"\n";
	while(w<=1.8){
		dt *= 10.;		
		vis = given[0] * (dt/(dx*dx));
		w = 1./(3.*vis+0.5);
	}
	acc = given[1] * ((dt*dt)/dx);
	int nrCellsX = cylinder[0] / dx;
	int ballCenterX = cylinder[2] / dx;
	int ballCenterY = cylinder[3] / dx;
	int ballDiameter = cylinder[4] / dx;
	vi result({dt, dx, acc, (double)nrCellsX, given[2], (double)ballCenterX, (double)ballCenterY, (double)ballDiameter });
	return result;
}

//Lattice Boltzmann

class LatticeB{
	public:
		LatticeB(vi param){
			forn(i, param.size())
			 cerr<<param[i]<<"\n";
		}

		/**
		* run the simulation 'nr' time
		*/
		void run(int nr){

		}

		/**
		* get the x component of velocity
		*/
		void getResult(){

		}

	private:


		int nrCellsY, nrCellsX, ballCenterX, ballCenterY, ballDiameter;
		double dx, dt, acc ;
		ma domain;
};