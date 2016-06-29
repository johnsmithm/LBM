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
 */
vi getParameters(vi given){
	vi result;

	return result;
}

//Lattice Boltzmann

class LatticeB{
	public:
		LatticeB(vi par){

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
		int nrCellsY, nrCellX;
		double dx, dt, acc ;
		ma domain;
};