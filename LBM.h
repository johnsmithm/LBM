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
#define forn(i,n) for(size_t i=0;i<(size_t)(n);++i)
#define forr(i,a,n) for(size_t i=a;i<=(size_t)n;++i)

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
			 dt=param[0];
			 dx=param[1];
			 acc=param[2];
			 nrCellsX=param[3]+1;
			 nrCellsY=param[4]+1;
			 ballCenterX=param[5];//??
			 ballCenterY=param[6];
			 ballDiameter=param[7];
			 neighbours = vector<pair<int,int> >({{1,-1},{1,0},{1,1},{0,1},{-1,1},{1,0},{-1,-1},{0,-1},{0,0}});
			 w = vector<double>({1./9, 1./36.,1./9, 1./36.,1./9, 1./36.,4./9.});
			 feq = vector<double>(9,0.);
		}

		/**
		* run the simulation 'nr' time
		*/
		void run(int nr){
			return;
			initiate();
			forn(i,nr){
				periodicBoundaryHandler();
				noSlipBoundaryHandler();
				streamStep();
				collideStep();
			}
		}

		/**
		* get the x component of velocity
		*/
		void getResult(){

		}

	private:
		void initiate(){
			domain.assign(nrCellsY+1, vector<vector<double> >(nrCellsX+1, vector<double>(9,0.)));
			domainHelper.assign(nrCellsY+1, vector<vector<double> >(nrCellsX+1, vector<double>(9,0.)));
			forr(i,1,nrCellsY)
				forr(j,1,nrCellsX){
					for(int k=0;k<9;++k)						 
						domain[i][j][k] = w[k];
					// 7 0 1
					// 6 8 2
					// 5 4 3
				}

		}
		void periodicBoundaryHandler(){
			forr(i,1,nrCellsY-1)
				forn(j,9){
					domain[i][0][j] = domain[i][nrCellsX-1][j];
					domain[i][nrCellsX][j] = domain[i][1][j];
				}
		}
		void noSlipBoundaryHandler(){
			forr(i,1,nrCellsX-1)
				forn(j,9){
					domain[0][i][j] = domain[nrCellsY-1][i][j];
					domain[nrCellsY][i][j] = domain[1][i][j];
				}

			//ToDo ballCell

		}
		void streamStep(){
			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue;
					forn(q,9)
						domainHelper[i][j][q] = domainHelper[i+neighbours[q].x][j+neighbours[q].y][q];					
				}

			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue;
					forn(q,9)
						domain[i][j][q] = domainHelper[i][j][q];
				}
		}
		void collideStep(){
			double q;
			pair<double , double> u = mp(0.,0.);

			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue
					q = 0.;
					forn(k,9)
						q+=domain[i][j][k];
					forn(k,9){
						u.x += domain[i][j][k]*neighbours[k].x;						
						u.y += domain[i][j][k]*neighbours[k].y;
						feq[k] = 0.;
					}
					u.x /= q;
					u.y /= q;
					forn(k,9){
						double pr = neighbours[k].x*u.x+neighbours[k].y*u.y;
						double u2 = u.x*u.x+u.y*u.y;
						feq[k] = w[k]*q*(1.+ 3*pr + ((9.*pr*pr)/2.)- (3.*u2)/2.);
						domain[i][j][k] -= w[k]*(domain[i][j][k]-feq[k])+3.*w[k]*q*(neighbours[k].x*acc+neighbours[k].y*acc);//??
					}
				}
		}

		int nrCellsY, nrCellsX, ballCenterX, ballCenterY, ballDiameter;
		double dx, dt, acc ;
		vector<pair<int,int> > neighbours;
		vi w;
		vi feq;
		ma domain, domainHelper;
};