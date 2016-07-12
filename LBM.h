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
	dt = 0.0005;
	vis = given[0] * (dt/(dx*dx));
	w = 1./(3.*vis+0.5);
	//cerr<<"w0:"<<w<<" vis0:"<<vis<<"\n";
	/*while(w>=1.9800){
		dt *= 10.;		
		vis = given[0] * (dt/(dx*dx));
		w = 1./(3.*vis+0.5);
		//cerr<<"w<"<<w<<"\n";
	}*/

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

			int iline, nrCellsX19 = (nrCellsX-1)*9,nrCellsX9=nrCellsX*9;
			int i9, nrCellsYline = nrCellsY*line;
			double q,feq1,uy,ux, eXu, u2;
			int nn ;

			//nr1 = 100;

			int ll=0;
			double qmax = 0., qmin = 2.;
			for(int j=0;j<9;++j)
				cerr<<"w["<<j<<"]="<<w[j]<<"\n";

			for(int i1=0;i1<nr1;++i1){
				// periodic boundary handler
				for(int i=1;i<nrCellsY;++i)
					for(int j=0;j<9;++j){
						iline = i*line;
						domain[j + (iline)] = domain[j + (iline) + (nrCellsX19)];
						domain[j + (iline) + (nrCellsX9)] = domain[j + 9 + (iline)];
					}
		
				// no slip boundary handler 
				for(int i=1;i<nrCellsX;++i)
					for(int j=0;j<3;++j){
						i9 = i*9;
						domain[(2 - j) + (i9+neighbours[j]*9)] = domain[(6 - j) + line + (i9)];
						domain[j + 4 + (nrCellsYline) + (i9+neighbours[j]*9)] = domain[j + (nrCellsYline-line) + (i9)];
					}

				for(int i=0;i<nk;i+=2)
					domain[rel[i]] = domain[rel[i+1]];

				//if(i1==90)  checkBoundary();

				for(int i=1;i<nrCellsY;++i)
					for(int j=1;j<nrCellsX;++j){
						if(ballCellsB[i*nrCellsX+j])continue;
						q = 0.;ux = 0.; uy = 0.;
					        nn = (i*line)+(j*9);
						for(int k=0;k<9;++k){
							// stream step
							feq1 = domain[k+(line*(i-neighbours[k+9]))+(9*(j-neighbours[k]))];
							domainHelper[k+nn] = feq1;
							q   += feq1;
							ux += (feq1*neighbours[k]);						
							uy += (feq1*neighbours[9+k]);
						}
						ux /= q;
						uy /= q;

						//test
						qmax = max(qmax,q);
						qmin = min(qmin,q);
						// collide step
						for(int k=0;k<9 ;++k){
						        eXu = (((double)neighbours[k])*ux)+(((double)neighbours[9+k])*uy);
							u2 = (ux*ux)+(uy*uy);
							feq1 = w[k]*q*(1.+ 3.*eXu + (4.5*eXu*eXu)- (1.5*u2));							
							domain[k+nn] = ((domainHelper[k+nn]*(1.-W)) + (W*feq1)+ (3.*w[k]*q*acc*(neighbours[k]*1.)));
						}
					}
					
				//if(i1==90)
				//  checkStreamCollide();
				
				++ll;
				if(i1%500==0)
				  cerr<<"iteration "<<i1<<"\n";
			}

			cerr<<"nr iteration:"<<ll<<"\n";
			cerr<<"qmax"<<qmax<< " qmin"<<qmin<<"\n";
		}

		void checkBoundary(){
		        ofstream ou("d1.out");
			int nnt = 0;
			size_t raza = (ballDiameter / 2) ;
			forr(i,ballCenterY-raza-1, ballCenterY+raza){
			  
				ou<<"i:"<<left<<setw(10)<<i;
				forr(j,ballCenterX - raza-1, ballCenterX + raza+1) ou<<left<<setw(12)<<j<<" ";
				ou<<"\n";
				
				forn(k,9){
				  ou<<left<<setw(12)<<k<<" ";
				  forr(j,ballCenterX - raza-1, ballCenterX + raza+1){
				                size_t ds = (0.5 + i-ballCenterY)*(0.5+i-ballCenterY) + ( -0.5 + j-ballCenterX)*(-0.5+j-ballCenterX);
						if(raza*raza >= ds)ou<<"*";
						else ou<<" ";
				                 ou<<left<<setw(12)<<domain[(i*line) + (j*9) + k];
				  }
				  ou<<" "<<k<<"\n";
				}
				ou<<"\n";
				
				forr(j,ballCenterX - raza-1, ballCenterX + raza+1)
					{  //cerr<<" j:"<<j<<" ";
						
					  
						assert(i>=1 && i<(size_t)nrCellsY && j >=1 && j<(size_t)nrCellsX);
						size_t ds = (0.5 + i-ballCenterY)*(0.5+i-ballCenterY) + ( -0.5 + j-ballCenterX)*(-0.5+j-ballCenterX);
						if(raza*raza > ds){
							
							//ballCells.insert(mp(i,j));
							//ballCellsB[j + i*nrCellsX] = true;
							
						}else{
							forn(k,9){

								if(ballCellsB[((i+neighbours[k+9])*(nrCellsX))+((j+neighbours[k]))]){
								  //cerr<<i<<" "<<neighbours[k+9]<<" "<<nrCellsX<<" " <<((i+neighbours[k+9])*(nrCellsX))<<" "<<(j+neighbours[k])<<"\n";
								  ++nnt;
									if(domain[(i*line) + (j*9) + k]!=domain[((i+neighbours[k+9])*line)+((j+neighbours[k])*9)+op[k]]){
										cerr<<"bounder error cylinder:"<< i<<" "<<j<<" "<<k<<"\n";
										cerr<<(i+neighbours[k+9])<<" "<<(j+neighbours[k])<<" "<<op[k]<<"\n";
										exit(0);
									}
								}
							}
						}
					}
					//cerr<<"\n";
				}
			cerr<<"nt:"<<nnt<<"\n";
		}
		void checkStreamCollide(){
		  ofstream ou("d22.out");
		  int nnt = 0;
		    size_t raza = (ballDiameter / 2) ;
			forr(i,ballCenterY-raza-1, ballCenterY+raza){
				ou<<"i:"<<left<<setw(10)<<i;
				forr(j,ballCenterX - raza-1, ballCenterX + raza+1) ou<<left<<setw(24)<<j<<" ";
				ou<<"\n";
				
				forn(k,9){
				  ou<<left<<setw(12)<<k<<" ";
				  forr(j,ballCenterX - raza-1, ballCenterX + raza+1){
				                size_t ds = (0.5 + i-ballCenterY)*(0.5+i-ballCenterY) + ( -0.5 + j-ballCenterX)*(-0.5+j-ballCenterX);
						if(raza*raza >= ds){ou<<"*";
						
				                 ou<<left<<setw(12)<<domainHelper[(i*line) + (j*9) + k]<<left<<setw(12)<<domain[(i*line) + (j*9) + k];}
						else{ ou<<" ";ou<<left<<setw(12)<<domainHelper[(i*line) + (j*9) + k]<<left<<setw(12)<<domain[(i*line) + (j*9) + k];}
				  }
				  ou<<" "<<k<<"\n";
				}
				ou<<"\n";
				
				forr(j,ballCenterX - raza-1, ballCenterX + raza+1)
					{   //cerr<<" j:"<<j<<" ";
						
					  
						assert(i>=1 && i<(size_t)nrCellsY && j >=1 && j<(size_t)nrCellsX);
						//size_t ds = (0.5 + i-ballCenterY)*(0.5+i-ballCenterY) + ( -0.5 + j-ballCenterX)*(-0.5+j-ballCenterX);
						if(ballCellsB[i*(nrCellsX) + j]){
							
							//ballCells.insert(mp(i,j));
							//ballCellsB[j + i*nrCellsX] = true;
							
						}else{
							forn(k,9){

								if(ballCellsB[((i+neighbours[k+9])*(nrCellsX))+((j+neighbours[k]))]){// boundary
									if(domainHelper[(i*line) + (j*9) + op[k]]!=domain[((i+neighbours[k+9])*line)+((j+neighbours[k])*9)+op[k]]){
										cerr<<"bounder collide error cylinder:"<< i<<" "<<j<<" "<<k<<"\n";
										cerr<<(i+neighbours[k+9])<<" "<<(j+neighbours[k])<<" "<<op[k]<<"\n";
										cerr<<domainHelper[(i*line) + (j*9) + op[k]]<<" "<<domain[((i+neighbours[k+9])*line)+((j+neighbours[k])*9)+op[k]]<<"\n";
										exit(0);
									}
									++nnt;
								}
							}
						}
					}
					//cerr<<"\n";
				}
				cerr<<"nnt collide:"<<nnt<<"\n";
				cerr<<"nk:"<<nk<<"\n";
				cerr<<"radius:"<<ballCells.size()<<"\n";
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
						u.x += domain[k+i*line+9*j]*neighbours[k];//??				
						u.y += domain[k+i*line+9*j]*neighbours[9+k];
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
						u.x += domain[k+i*line+9*j]*neighbours[k];//??				
						u.y += domain[k+i*line+9*j]*neighbours[9+k];
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

			img.save((ballCenterY==48?"scenario1.png":"scenario2.png"));
			cerrr.close();					
		}	

	private:
		void initiate(){
			domain       = new double[(nrCellsY+1)*(nrCellsX+1)*9];
			domainHelper = new double[(1+nrCellsY)*(nrCellsX+1)*9];
			ballCellsB   = new bool[(1+nrCellsY)*(nrCellsX+1)];
			cerr<<"total:"<<((nrCellsY+1)*(nrCellsX+1)*9)<<" "<<((1+nrCellsY)*(nrCellsX+1))<<"\n";
			line = (nrCellsX+1)*9;
			cerr<<"line"<<line<<"\n";
			ballCenterX += 3;
			size_t raza = (ballDiameter / 2) ;
			forr(i,ballCenterY-raza-1, ballCenterY+raza){
				cerr<<i<<" ";
				int nn = 0;
				forr(j,ballCenterX - raza-1, ballCenterX + raza+1)
					{
						assert(i>=1 && i<(size_t)nrCellsY && j >=1 && j<(size_t)nrCellsX);
						size_t ds = (0.5 + i-ballCenterY)*(0.5+i-ballCenterY) + ( -0.5 + j-ballCenterX)*(-0.5+j-ballCenterX);
						if(raza*raza >= ds){
							if(nn==0)cerr<<j<<"|";
							++nn;
							ballCells.insert(mp(i,j));
							ballCellsB[j + (i*nrCellsX)] = true;
							cerr<<"0";
						}else cerr<<"1";
					}
					cerr<<" "<<nn<<"\n";
				}

			forr(i,1,nrCellsY-1){
				forr(j,1,nrCellsX-1){
					if(ballCellsB[i*nrCellsX+j])continue;
					for(int k=0;k<9;++k)						 
						domain[k + i*line + j*9] = w[k];
						// 0 1 2
						// 7 8 3
						// 6 5 4
				}
			}
			int kk = 0;
			for(auto p : ballCells){// px - Y, py - X
			  bool ok = 0;
					forn(k,9)
						if(ballCells.count(mp(p.x+neighbours[9+k],p.y+neighbours[k]))==0){
							rel[nk] = k + (p.x*line) + (p.y*9);
							rel[nk+1] = ((k+4)%8) + ((p.x + neighbours[9+k])*line)+ ((p.y+neighbours[k])*9);
							assert(((k+4)%8) == op[k]);
							nk += 2;
							ok = 1;
						}
					if(ok)++kk;
				}
				cerr<<"kk:"<<kk<<"\n";
		}
		

		int nrCellsY, nrCellsX, ballCenterX, ballCenterY, ballDiameter, line, nk=0;
		double dx, dt, acc , W; int rel[5000], rell[5000];
		int neighbours[18] =  {-1,0,1,1,1,0,-1,-1,0,  1,1,1,0,-1,-1,-1,0,0};
		size_t op[9] = {4,5,6,7,0,1,2,3,8};
		double  w[9] = {(1./36.),(1./9.), (1./36.),(1./9.), (1./36.),(1./9.),( 1./36.),(1./9.), (4./9.)};
		
		double * domain, * domainHelper;
		bool * ballCellsB;
		set<pair<int,int> > ballCells;
};