#include <math.h>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <fstream>
#include <list>
#include <set>
#include <iomanip>

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
	while(w>=1.8){
		dt *= 10.;		
		vis = given[0] * (dt/(dx*dx));
		w = 1./(3.*vis+0.5);
		//<<"w<"<<w<<"\n";
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
			 W = param[8];	                	// lt 0    ct1  rt2   3rc   4rb     5cb    6lb    7 lc    cc
			 cerr<<"W:"<<W<<"\n";
			 neighbours = vector<pair<int,int> >({{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{0,0}});// x,y
			 w = vector<double>({1./36.,1./9, 1./36.,1./9, 1./36.,1./9, 1./36.,1./9, 4./9.});
			 feq = vector<double>(9,0.);
			
		}

		/**
		* run the simulation 'nr' time
		*/
		void run(int nr){
			
			initiate();
			//return;
			forn(i,nr){
				periodicBoundaryHandler();
				
				noSlipBoundaryHandler();//
				//return;
				#ifdef DEBUG
					if(i==1)//debug
						showTest("it2",0,0,5);
				#endif
				streamStep();//
			//	return;
				#ifdef DEBUG
					if(i==1)//debuf
						showTest("it2a",0,0,5);
				#endif
				collideStep();
				if(i%100==0)
				  cerr<<"iteration "<<i<<"\n";
				//return;
			}
		}

		/**
		* get the x component of velocity
		*/
		void getResult(){
			ofstream cerrr("result.out");

			pair<double , double> u = mp(0.,0.);
			//map<double,int> dd;
			//int nr = 0;

			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX -1){
					// if ball continue
					double q = 0.;
					u = mp(0.,0.);
					forn(k,9)
						q+=domain[i][j][k];
					//assert(q>0.5 && q<1.5);
					//cerr<<"q="<<q<<"\n";
					forn(k,9){

						//cerr<<k<<" "<<domain[i][j][k]<<" n1:"<<neighbours[k].x<<" "<<neighbours[k].y<<"\n";
						u.x += domain[i][j][k]*neighbours[k].x;//??				
						u.y += domain[i][j][k]*neighbours[k].y;
					}

					u.x /= q;
					//cerrr<<i<<" "<<j<<" ";
					if(ballCells.count(mp(i,j))==0)
					cerrr<<setw(13)<<left<<u.x<<" ";
					else cerrr<<setw(13)<<left<<"0.0"<<" ";
					//cerrr<<"\n";
					//if(ballCells.count(mp(i,j))!=0) continue;

					//if(dd.count(u.x)==0)dd[u.x]= nr++;

					
					//cerrr<<dd[u.x]<<" ";else cerrr<< "--";//plot with gnuplot

					

					//return;
				}
				cerrr<<"\n";
			}

			//for(auto p : dd)
				//cerrr<<p.y<<"-"<<p.x<<"|";

/*
			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX-1){
					double q = 0.;
					u = mp(0.,0.);
					forn(k,9)
						q+=domain[i][j][k];
					//assert(q>0.5 && q<1.5);
					//cerr<<"q="<<q<<"\n";
					forn(k,9){

						//cerr<<k<<" "<<domain[i][j][k]<<" n1:"<<neighbours[k].x<<" "<<neighbours[k].y<<"\n";
						u.x += domain[i][j][k]*neighbours[k].x;//??				
						u.y += domain[i][j][k]*neighbours[k].y;
					}
					u.x /= q;

					if(dd.count(u.x)==0)
						cerr<<"-";
					else
						cerr<<((int)(((double)dd[u.x]/nr)*10));
				}
				cerr<<"\n";
			}*/

			cerrr.close();
					
		}

		void showTest(string name, int bY, int bX, int eX){
			ofstream out(name);
			out << setprecision(6) << fixed;
			forr(i,5,10)
			{
				forn(k,9){
					out<<k<<"| ";

					forr(j,25,35){
						out<<domain[i][j][k];
						if(ballCells.count(mp(i,j))!=0)out<<"* ";
						else out<<"  ";
					}
					//out<<"---";
				//	forr(j,nrCellsX-5,nrCellsX-1)out<<domain[i][j][k]<<" ";
					out<<"|"<<k<<"\n";
				}
				out<<"\n";
			}
			out<<"----------------------------------------------------\n";
			forr(i,nrCellsY-4,nrCellsY)
			{
				forn(k,9){
					out<<k<<"| ";
					forr(j,bX,eX)out<<domain[i][j][k]<<" ";
					out<<"---";
					forr(j,nrCellsX-5,nrCellsX-1)out<<domain[i][j][k]<<" ";
					out<<"|"<<k<<"\n";
				}
				out<<"\n";
			}
		}

	private:
		void initiate(){
			domain.assign(nrCellsY+1, vector<vector<double> >(nrCellsX+1, vector<double>(9,0.)));
			domainHelper.assign(nrCellsY+1, vector<vector<double> >(nrCellsX+1, vector<double>(9,0.)));
			

			size_t raza = (ballDiameter / 2) ;
			forr(i,ballCenterY-raza, ballCenterY+raza)
				forr(j,ballCenterX - raza, ballCenterX + raza)
					{
						assert(i>=1 && i<(size_t)nrCellsY && j >=1 && j<(size_t)nrCellsX);
						size_t ds = (i-ballCenterY)*(i-ballCenterY) + (j-ballCenterX)*(j-ballCenterX);
						if(raza*raza >= ds){
							ballCells.insert(mp(i,j));
						}
					}

			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX-1){
					for(int k=0;k<9;++k)						 
						domain[i][j][k] = w[k];
					// 0 1 2
					// 7 8 3
					// 6 5 4
					/*#ifdef DEBUG
						if(ballCells.count(mp(i,j))==0)cerr<<1;
						else cerr<<0;
					#endif*/
				}
				#ifdef DEBUG
					//cerr<<"\n";
				#endif
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
				forn(j,3){

					domain[0][i+neighbours[j].x][2-j] = domain[1][i][6-j];//first

					domain[nrCellsY][i+neighbours[j].x][4+j] = domain[nrCellsY - 1][i][j];//last
				}

			//ballCells
			for(auto p : ballCells){// px - Y, py - X
				forn(k,9)
					if(ballCells.count(mp(p.x+neighbours[k].y,p.y+neighbours[k].x))==0){
						//cerr<<"cell;"<<p.x<<" "<<p.y<<" ballN:"<<(p.x+neighbours[k].y)<<" "<<p.y+neighbours[k].x<<"\n";
						//size_t id = (k%2==1?(k+4)%8:(k+4)%8);
						domain[p.x][p.y][k] = domain[p.x + neighbours[k].y][p.y+neighbours[k].x][(k+4)%8];						
					}
			}

		}
		void streamStep(){
			// pull stream
			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue;
					if(ballCells.count(mp(i,j))!=0)continue;
					forn(q,9){
						//if(i == 1 && j == 1)cerr<<domainHelper[i][j][q]<<" - ";
						domainHelper[i][j][q] = domain[i-neighbours[q].y][j-neighbours[q].x][q];	
						//if(i == 1 && j == 1)cerr<<domainHelper[i][j][q]<<"a\n";
					}

				}

			//copy values
			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue;
					forn(q,9)
						domain[i][j][q] = domainHelper[i][j][q];
				}
			//cerr<<"test:"<<domain[10][10][5]<<"\n";
		}
		void collideStep(){//{{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1},{-1,0},{0,0}}
			double q;
			pair<double , double> u = mp(0.,0.);

			forr(i,1,nrCellsY-1)
				forr(j,1,nrCellsX-1){
					// if ball continue
					if(ballCells.count(mp(i,j))!=0)continue;
					q = 0.;
					u = mp(0.,0.);
					forn(k,9)
						q+=domain[i][j][k];
					forn(k,9){//correct
						u.x += domain[i][j][k]*neighbours[k].x;						
						u.y += domain[i][j][k]*neighbours[k].y;
						feq[k] = 0.;
					}
					u.x /= q;
					u.y /= q;
					forn(k,9){
						double pr = neighbours[k].x*u.x+neighbours[k].y*u.y;
						double u2 = u.x*u.x+u.y*u.y;
						feq[k] = w[k]*q*(1.+ 3.*pr + ((9.*pr*pr)/2.)- (3.*u2)/2.);
						
						//if((i == 1||i==(size_t)nrCellsY) && j == 1)cerr<<i<<" "<<feq[k]<<"f - d:";
						domain[i][j][k] = domain[i][j][k] - W*(domain[i][j][k]-feq[k])+3.*w[k]*q*(neighbours[k].x!=0?acc*neighbours[k].x+acc*0:0);//(neighbours[k].y*0+neighbours[k].x*acc);//??
						//if((i == 1||i==(size_t)nrCellsY) && j == 1)cerr<<i<<" "<<domain[i][j][k]<<"\n";
					}

					//test
					#ifdef DEBUG
						q = 0.;
						//u = mp(0.,0.);

						forn(k,9)
							q+=domain[i][j][k];

						if(q>1.5 || q<0.3){cerr<<q<<" error\n"; 
						cerr<<"v:"<<u.x<<" "<<u.y<<"\n";
						  forn(k,9) cerr<<domain[i][j][k]<<" f:"<<feq[k]<<"\n";
						return;}

						u = mp(0.,0.);
					forn(k,9){//correct
						u.x += domain[i][j][k]*neighbours[k].x;						
						u.y += domain[i][j][k]*neighbours[k].y;
						feq[k] = 0.;
					}
					if(u.x < 0 ){
						cerr<<i<<" "<<j<<" - reverse velosity\n";
						//exit(0);
						getResult();
						forn(k,9) cerr<<"b:"<<domainHelper[i][j][k]<<" a:"<<domain[i][j][k]<<" f:"<<feq[k]<<"\n";
						exit(0);
					}
					#endif
				}
		}

		int nrCellsY, nrCellsX, ballCenterX, ballCenterY, ballDiameter;
		double dx, dt, acc , W;
		vector<pair<int,int> > neighbours;
		vi w;
		vi feq;
		ma domain, domainHelper;
		set<pair<int,int> > ballCells;
};