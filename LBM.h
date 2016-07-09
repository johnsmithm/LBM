#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <fstream>
#include <list>
#include <set>
#include <iomanip>

#include "imageClass/GrayScaleImage.h"

using namespace std;



typedef vector<vector<vector<double> > > ma;
typedef vector<double > vi;

#define x first
#define y second
#define pb push_back
#define mp make_pair
#define forn(i,n) for(size_t i=0;i<(size_t)(n);++i)
#define forr(i,a,n) for(size_t i=a;i<=(size_t)n;++i)

#define DEBUG

/**
 * get the lattice parameter from real parameters
 *
 * @param given {viscosity,accceleration,resolution of cylinder diameter}
 * @param cylinder{sizeX, sizeY, ballCenterX, ballCenterY, ballDiameter}
 * return {delta t, delta x, acceleration, number cell X, number cell Y, ballCenterCellX, ballCenterCellY, ballDiameter}
 */
vi getParameters(vi given, vi cylinder){
	
	double dt, dx, w, acc, vis;
	dx = cylinder[4] / given[2];
	dt = 0.000000001;
	vis = given[0] * (dt/(dx*dx));
	w = 1./(3.*vis+0.5);
	//cerr<<"w0:"<<w<<" vis0:"<<vis<<"\n";
	while(w>=1.9800){
		dt *= 10.;		
		vis = given[0] * (dt/(dx*dx));
		w = 1./(3.*vis+0.5);
		//cerr<<"w<"<<w<<"\n";
	}

	#ifdef DEBUG
		cerr<<"w:"<<w<<" vis:"<<vis<<"\n";
		cerr<<"dt:"<<dt<<" dx:"<<dx<<"\n";
	#endif

	acc = given[1] * ((dt*dt)/dx);
	int nrCellsX = cylinder[0] / dx;
	int nrCellsY = cylinder[1] / dx;
	int ballCenterX = cylinder[2] / dx;
	int ballCenterY = cylinder[3] / dx;
	int ballDiameter = cylinder[4] / dx;

	#ifdef DEBUG
		cerr<<"ballDiameter:"<<ballDiameter<<" x:"<<ballCenterX<<" y:"<<ballCenterY<<"\n";
		cerr<<"domain:"<<nrCellsX<<"X"<<nrCellsY<<"\n";
	#endif

	vi result({dt, dx, acc, (double)nrCellsX, (double)nrCellsY, (double)ballCenterX, (double)ballCenterY, (double)ballDiameter, w });
	return result;
}

//Lattice Boltzmann

class LatticeB{
	public:
		LatticeB(vi param){
			//forn(i, param.size())
				//cerr<<param[i]<<"\n";
			 dt=param[0];
			 dx=param[1];
			 acc=param[2];

			 #ifdef DEBUG
			 	cerr<<"acc:"<<acc<<"\n";
			 #endif

			 nrCellsX=param[3]+1;
			 nrCellsY=param[4]+1;
			 ballCenterX=param[5];//??
			 ballCenterY=param[6];
			 ballDiameter=param[7];		
			 W = param[8];
			 cerr<<"W:"<<W<<"\n";
			
			
		}

		/**
		* run the simulation 'nr' time
		*/
		void run(int nr1){
			
			initiate();
			//nr1 = 1000;
			for(int i1=0;i1<nr1;++i1){
				// periodic boundary handler
				int iline, nrCellsX19 = (nrCellsX-1)*9,nrCellsX9=nrCellsX*9;
				for(int i=1;i<nrCellsY;++i)
					for(int j=0;j<9;++j){
						iline = i*line;
						domain[j + (iline)] = domain[j + (iline) + (nrCellsX19)];
						domain[j + (iline) + (nrCellsX9)] = domain[j + 9 + (iline)];
					}
		
				// no slip boundary handler 
				int i9, nrCellsYline = nrCellsY*line;
				for(int i=1;i<nrCellsX;++i)
					for(int j=0;j<3;++j){
						i9 = i*9;
						domain[(2 - j) + (i9+neighbours[j*2]*9)] = domain[(6 - j) + line + (i9)];
						domain[j + 4 + (nrCellsYline) + (i9+neighbours[j*2]*9)] = domain[j + (nrCellsYline-line) + (i9)];
					}

				/*for(auto p : ballCells){// px - Y, py - X
					forn(k,9)
						if(ballCells.count(mp(p.x+neighbours[1+2*k],p.y+neighbours[2*k]))==0){
							domain[k + (p.x*line) + (p.y*9)] = domain[((k+4)%8) + ((p.x + neighbours[1+2*k])*line)+ ((p.y+neighbours[2*k])*9)];					
						}
				}*/
				for(int i=0;i<nk;i+=2)
					domain[rel[i]] = domain[rell[i+1]];

			    //{{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{0,0}}
				double q,feq1,uy,ux;
				int nn ;

				for(int i=1;i<nrCellsY;++i)
					for(int j=1;j<nrCellsX;++j){
						if(ballCellsB[i*nrCellsX+j])continue;
						q = ux = uy = 0.;
					    nn = (i*line)+(j*9);
						for(int k=0;k<9;++k){
							// stream step
							feq1 = domain[k+(line*(i-neighbours[k*2+1]))+(9*(j-neighbours[k*2]))];
							domainHelper[k+nn] = feq1;
							q   += feq1;
							ux += feq1*neighbours[k*2];						
							uy += feq1*neighbours[1+k*2];
						}
						ux /= q;
						uy /= q;
						// collide step
						for(int k=0;k<9;++k){
							double pr = neighbours[k*2]*ux+neighbours[1+k*2]*uy;
							double u2 = ux*ux+uy*uy;
							feq1 = w[k]*q*(1.+ 3.*pr + (4.5*pr*pr)- (1.5*u2));
							double pr1 = domainHelper[k+nn];
							domain[k+nn] = pr1*(1.-W) + W*feq1+ 3.*w[k]*q*acc*neighbours[k*2];
						}				
					}

				if(i1%500==0)
				  cerr<<"iteration "<<i1<<"\n";
			}
		}

		/**
		* get the x component of velocity
		*/
		void getResult(){
			ofstream cerrr("result1.out");
			GrayScaleImage img(nrCellsX-1,nrCellsY-1);
			cerr<<"img:"<<img.width()<<"x"<<img.height()<<"\n";
			pair<double , double> u = mp(0.,0.);
			double mi = 0.00000, ma = 0.0000010 ;
			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX -1){
					double q = 0.;
					u = mp(0.,0.);
					forn(k,9)
						q+=domain[k+i*line+9*j];
					forn(k,9){
						u.x += domain[k+i*line+9*j]*neighbours[2*k];//??				
						u.y += domain[k+i*line+9*j]*neighbours[1+2*k];
					}
					u.x /= q;
					if(ballCells.count(mp(i,j))==0){
						cerrr<<setw(13)<<left<<u.x<<" ";
						if(u.x > ma)ma = u.x;
						if(u.x < mi)mi = u.x;
					}
					else {
						cerrr<<setw(13)<<left<<"0.0"<<" ";
						
					}

					//img.setElement(j-1,i-1,(u.x+0.1));
				}
				cerrr<<"\n";
			}
			cerr<<"max:"<<ma<<" min:"<<mi<<"\n";
			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX -1){
					double q = 0.;
					u = mp(0.,0.);
					forn(k,9)
						q+=domain[k+i*line+9*j];
					forn(k,9){
						u.x += domain[k+i*line+9*j]*neighbours[2*k];//??				
						u.y += domain[k+i*line+9*j]*neighbours[1+2*k];
					}
					u.x /= q;
					int val = 0;
					if(ballCells.count(mp(i,j))==0){
						
					}
					else {
						u.x = 0.;
					}
					val = (255 - ((u.x - mi)/(ma-mi))*255);
					img.getElement(j-1,i-1) = val;
				}
				cerrr<<"\n";
			}

			img.save((ballCenterY==48?"scenario12.png":"scenario22.png"));
			cerrr.close();					
		}	

	private:
		void initiate(){
			domain       = new double[(nrCellsY+1)*(nrCellsX+1)*9];
			domainHelper = new double[(1+nrCellsY)*(nrCellsX+1)*9];
			ballCellsB   = new bool[(1+nrCellsY)*(nrCellsX+1)];
			line = (nrCellsX+1)*9;

			size_t raza = (ballDiameter / 2) ;
			forr(i,ballCenterY-raza, ballCenterY+raza)
				forr(j,ballCenterX - raza, ballCenterX + raza)
					{
						assert(i>=1 && i<(size_t)nrCellsY && j >=1 && j<(size_t)nrCellsX);
						size_t ds = (i-ballCenterY)*(i-ballCenterY) + (j-ballCenterX)*(j-ballCenterX);
						if(raza*raza >= ds){
							ballCells.insert(mp(i,j));
							ballCellsB[j + i*nrCellsX] = true;
						}
					}

			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX-1){
					for(int k=0;k<9;++k)						 
						domain[k + i*line + j*9] = w[k];
						// 0 1 2
						// 7 8 3
						// 6 5 4
				}
			}

			for(auto p : ballCells){// px - Y, py - X
					forn(k,9)
						if(ballCells.count(mp(p.x+neighbours[1+2*k],p.y+neighbours[2*k]))==0){
							//domain[k + (p.x*line) + (p.y*9)] = domain[((k+4)%8) + ((p.x + neighbours[1+2*k])*line)+ ((p.y+neighbours[2*k])*9)];					
							rel[nk] = k + (p.x*line) + (p.y*9);
							rell[nk+1] = ((k+4)%8) + ((p.x + neighbours[1+2*k])*line)+ ((p.y+neighbours[2*k])*9);
							nk += 2;
						}
				}
		}
		

		int nrCellsY, nrCellsX, ballCenterX, ballCenterY, ballDiameter, line, nk=0;
		double dx, dt, acc , W; int rel[5000], rell[5000];
		int neighbours[18] =  {-1,1,0,1,1,1,1,0,1,-1,0,-1,-1,-1,-1,0,0,0};

		double  w[9] = {1./36.,1./9., 1./36.,1./9., 1./36.,1./9., 1./36.,1./9., 4./9.};
		
		double * domain, * domainHelper;
		bool * ballCellsB;
		set<pair<int,int> > ballCells;
};