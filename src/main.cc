#include <iostream>
#include "FileReader.hh"
#include <string>
#include "Grid.hh"
#include <fstream>
#include <sys/time.h>

enum cellflag {BOUNDARY,ACCELERATION,FLUID};


static void printmethod(lbm::V_Field &velocity, lbm::D_Field &density, size_t sizex,size_t sizey, size_t k){
			//std::string filename;
			//filename = "test" + std::to_string(k)  +".dat";
			std::ofstream vtkFile( "ldc" + std::to_string(k)  +".vtk");

			if(vtkFile.is_open()){
				vtkFile << "# vtk DataFile Version 4.0 \nSiwiRVisFile \nASCII \nDATASET STRUCTURED_POINTS \nDIMENSIONS " << sizex << " ";
				vtkFile << sizey << " " << 1 << "\n" << "ORIGIN 0 0 0"<< "\n";
				vtkFile << "SPACING 1 1 1\nPOINT_DATA " << sizex*sizey << "\n" << "\n" ;
				vtkFile << "SCALARS flags double 1 \nLOOKUP_TABLE default" << "\n" ;
				
				for(int y=2;y<sizex+2;++y){
					for(int x=2;x<sizey+2;++x){
						vtkFile << 1 << "\n";
					}

				}
				
				vtkFile << "\nSCALARS density double 1 \nLOOKUP_TABLE default \n" ;
				
				for(int y=2;y<sizex+2;++y){
					for(int x=2;x<sizey+2;++x){
						vtkFile << density(x,y,0)<< "\n";
					}

				}
			
				vtkFile << "\nVECTORS velocity double\n" ;
				
				for(int y=2;y<sizex+2;++y){
					for(int x=2;x<sizey+2;++x){
						vtkFile << velocity(x,y,0) << " " << velocity(x,y,1) << " " << 0 <<"\n";
					}

				}
				vtkFile.close();
			}

}

static void initializer(lbm::PDF_Field &field,int sizex, int sizey){
	// Initialize Field
	for(int i=0;i<sizex;++i){
		for(int j=0;j<sizey;++j){
			field(i,j,0)=4.0/9.0;
			field(i,j,1)=1.0/9.0;
			field(i,j,2)=1.0/36.0;
			field(i,j,3)=1.0/9.0;
			field(i,j,4)=1.0/36.0;
			field(i,j,5)=1.0/9.0;
			field(i,j,6)=1.0/36.0;
			field(i,j,7)=1.0/9.0;
			field(i,j,8)=1.0/36.0;

			//std::cout << field(0,0,0) << std::endl;
	

		}
	}
}

static void flags_initializer(lbm::Flags &flags,int sizex, int sizey){
	// Initialize Flag-Field
	for(int y=1;y<sizex+3;++y){
		for(int x=1;x<sizey+3;++x){
			if(y==1 || (x==1 && y!=sizey+2) || (x==sizex+2 && y!=sizey+2)){
				flags(x,y,0)=BOUNDARY;
			} else if (y==sizey+2){
				flags(x,y,0)=ACCELERATION;
			} else {
				flags(x,y,0)=FLUID;
			}
		}
	}
}

int main(int argc, char **args)
{
	if(argc != 2)
	{
		std::cout << "Usage: ./lbm <params.dat>" << std::endl;
		return 0;
	}	
	//read in parameters
	FileReader reader(args[1]);
	std::cout << reader.getParameter<size_t>("vtk_step") << std::endl;
	std::cout << reader.getParameter<std::string>("vtk_file") << std::endl;
	std::cout << reader.getParameter<int>("timesteps") << std::endl;
	std::cout << reader.getParameter<double>("omega") << std::endl;
/*	
	const size_t timesteps(reader.getParameter<size_t>("timesteps"));
	std::cout << timesteps << std::endl;

	const size_t sizex(reader.getParameter<size_t>("sizex"));
	const size_t sizey(reader.getParameter<size_t>("sizey"));
*/	const double omega(reader.getParameter<double>("omega"));
	std::cout << "omega: " <<omega << std::endl;
	//const size_t vtk_step(reader.getParameter<size_t>("vtk_step") );
	const std::string vtk_file(reader.getParameter<std::string>("vtk_file"));
	const size_t sizex(5);
	const size_t sizey(5);
	const double timesteps(4);


	// define weighting factors t_alpha
	const double t_0=4.0/9.0;
	const double t_1=1.0/9.0;
	const double t_2=1.0/36.0;
	lbm::PDF_Field src(sizex+4,sizey+4);
	lbm::PDF_Field dst(sizex+4,sizey+4);
	lbm::V_Field velocity(sizex+4,sizey+4);
	lbm::D_Field density(sizex+4,sizey+4);
	lbm::Flags flags(sizex+4,sizey+4);
	
	// initialize Fluid-Field
	initializer(src,sizex+4,sizey+4);
	// initialize Flag-Field
	flags_initializer(flags,sizex,sizey);
	struct timeval t0, t;

	gettimeofday(&t0, NULL);

	for(int k=0;k<timesteps;++k){
		std::cout << "SRC: " << src(0,sizey+3,4) << "DST: " << dst(0,sizey+3,4) << std::endl;
		//if(k==0){
			

		double rho_sum=0.0;
		for(int y=1;y<sizex+3;++y){
			for(int x=1;x<sizey+3;++x){
			//bottom boundary
			/*
			if(y==0 && x!=0 && x!=sizex+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,5)=src(x,y+1,5);
				dst(x,y,6)=src(x+1,y+1,6);
				dst(x,y,4)=src(x-1,y+1,4);
				dst(x,y,3)=src(x-1,y,3);
				dst(x,y,7)=src(x+1,y,7);
				std::cout << "boundary" << std::endl;

			}
			//top boundary
			else if(y==sizey+1 && x!=0 && x!=sizex+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,2)=src(x-1,y-1,2);
				dst(x,y,8)=src(x+1,y-1,8);
				dst(x,y,3)=src(x-1,y,3);
				dst(x,y,7)=src(x+1,y,7);	
				//dst(x,y,5)=1000;
				std::cout << "top" << std::endl;


			}
			//left boundary
			else if(x==0 && y!=0 && y!=sizey+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,7)=src(x+1,y,7);
				dst(x,y,8)=src(x+1,y-1,8);
				dst(x,y,6)=src(x+1,y+1,6);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,5)=src(x,y+1,5);

			}
			//right boundary
			else if(x==sizex+1 && y!=0 && y!=sizey+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,2)=src(x-1,y-1,2);
				dst(x,y,3)=src(x-1,y,3);
				dst(x,y,4)=src(x-1,y+1,4);
				dst(x,y,5)=src(x,y+1,5);

			}
			// bottom left corner
			else if(x==0 && y==0){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,5)=src(x,y+1,5);
				dst(x,y,6)=src(x+1,y+1,6);
				dst(x,y,7)=src(x+1,y,7);
			}
			// top left corner
			else if(x==0 && y==sizey+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,7)=src(x+1,y,7);
				dst(x,y,8)=src(x+1,y-1,8);
				dst(x,y,1)=src(x,y-1,1);
				//dst(x,y,5)=1000;
			}
			//bottom right corner
			else if(x==sizex+1 && y==0){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,5)=src(x,y+1,5);
				dst(x,y,4)=src(x-1,y+1,4);
				dst(x,y,3)=src(x-1,y,3);
				std::cout << "bottom right" << std::endl;
			}
			//top right corner
			else if(x==sizex+1 && y==sizey+1){
				dst(x,y,0)=src(x,y,0);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,2)=src(x-1,y-1,2);
				dst(x,y,3)=src(x-1,y,3);
				std::cout << "top right" << std::endl;

			}

			
			else{
			*/	// stream-step on fluid-cells and boundaries
				dst(x,y,0)=src(x,y,0);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,2)=src(x-1,y-1,2);
				dst(x,y,3)=src(x-1,y,3);
				dst(x,y,4)=src(x-1,y+1,4);
				dst(x,y,5)=src(x,y+1,5);
				dst(x,y,6)=src(x+1,y+1,6);
				dst(x,y,7)=src(x+1,y,7);
				dst(x,y,8)=src(x+1,y-1,8);
					//	}		
			//collide boundaries
			if(flags(x,y,0)==BOUNDARY){

			//	dst(x,y,0)=src(x,y,0);
				std::swap(dst(x,y,1),dst(x,y,5));
				std::swap(dst(x,y,2),dst(x,y,6));
				std::swap(dst(x,y,3),dst(x,y,7));
				std::swap(dst(x,y,4),dst(x,y,8));
				
				//std::cout << "collide boundary " << x << " " <<y << std::endl;
				continue;

			}

			double rho;
			double u_x;
			double u_y;

			//collide acceleration cell 
			if(flags(x,y,0)==ACCELERATION){
				//std::cout << "acc" << std::endl;
				u_x=0.08;	
				u_y=0.0;
	
				dst(x,y,1)=dst(x,y,1)-6.0*t_1*(u_y);
				dst(x,y,2)=dst(x,y,2)-6.0*t_2*(u_x+u_y);
				dst(x,y,3)=dst(x,y,3)-6.0*t_1*(u_x);
				dst(x,y,4)=dst(x,y,4)-6.0*t_2*(u_x-u_y);
				dst(x,y,5)=dst(x,y,5)-6.0*t_1*(-u_y);
				dst(x,y,6)=dst(x,y,6)-6.0*t_2*(-u_y-u_x);
				dst(x,y,7)=dst(x,y,7)-6.0*t_1*(-u_x);
				dst(x,y,8)=dst(x,y,8)-6.0*t_2*(-u_x+u_y);
				
				std::swap(dst(x,y,1),dst(x,y,5));
				std::swap(dst(x,y,2),dst(x,y,6));
				std::swap(dst(x,y,3),dst(x,y,7));
				std::swap(dst(x,y,4),dst(x,y,8));
			//	std::cout << "collide acc " << x << " " <<y << std::endl;
				//std::cout << "density acc" << dst(x,y,0)+dst(x,y,1)+dst(x,y,2)+dst(x,y,3)+dst(x,y,4)+dst(x,y,5)+dst(x,y,6)+dst(x,y,7)+dst(x,y,8)<<std::endl;
				//std::cout << "density acc single" << dst(x,y,0) << " "<<dst(x,y,1) << " "<<dst(x,y,2) << " "<<dst(x,y,3) << " "<<dst(x,y,4) << " "<<dst(x,y,5) << " "<<dst(x,y,6) << " "<<dst(x,y,7) << " "<<dst(x,y,8)<<std::endl;

				continue;


			}
		
			//collide regular fluid cell
			
			u_x=dst(x,y,2)+dst(x,y,3)+dst(x,y,4);
			u_y=dst(x,y,8)+dst(x,y,1)+dst(x,y,2);
			rho=dst(x,y,5)+dst(x,y,6)+dst(x,y,7)+dst(x,y,0)-dst(x,y,2)+u_x+u_y;
			u_x=(u_x-dst(x,y,6)-dst(x,y,7)-dst(x,y,8));//rho;
			u_y=(u_y-dst(x,y,4)-dst(x,y,5)-dst(x,y,6));//rho;
			
			//std::cout << rho << std::endl;
			std::cout<< "velo1" << x <<" "<< y << " "<< u_x << " " << u_y << std::endl;
		//	velocity(x,y,0)=u_x;
		//	velocity(x,y,1)=u_y;
		//	density(x,y,0)=rho;	
			
			rho_sum+=rho;
		
			double u_sqr=1.5*(u_x*u_x+u_y*u_y);
			
			dst(x,y,0)=(1.0-omega)*dst(x,y,0)+t_0*omega*(rho-u_sqr);
			dst(x,y,1)=(1.0-omega)*dst(x,y,1)+t_1*omega*(rho+3.0*u_y+4.5*(u_y*u_y)-u_sqr);
			dst(x,y,2)=(1.0-omega)*dst(x,y,2)+t_2*omega*(rho+3.0*(u_x+u_y)+4.5*((u_x+u_y)*(u_x+u_y))-u_sqr);
			dst(x,y,3)=(1.0-omega)*dst(x,y,3)+t_1*omega*(rho+3.0*u_x+4.5*(u_x*u_x)-u_sqr);
			dst(x,y,4)=(1.0-omega)*dst(x,y,4)+t_2*omega*(rho+3.0*(u_x-u_y)+4.5*((u_x-u_y)*(u_x-u_y))-u_sqr);
			dst(x,y,5)=(1.0-omega)*dst(x,y,5)+t_1*omega*(rho-3.0*u_y+4.5*(u_y*u_y)-u_sqr);
			dst(x,y,6)=(1.0-omega)*dst(x,y,6)+t_2*omega*(rho+3.0*(-u_x-u_y)+4.5*((-u_x-u_y)*(-u_x-u_y))-u_sqr);
			dst(x,y,7)=(1.0-omega)*dst(x,y,7)+t_1*omega*(rho-3.0*u_x+4.5*(u_x*u_x)-u_sqr);
			dst(x,y,8)=(1.0-omega)*dst(x,y,8)+t_2*omega*(rho+3.0*(u_y-u_x)+4.5*((-u_x+u_y)*(-u_x+u_y))-u_sqr);
			/*
			u_x=dst(x,y,2)+dst(x,y,3)+dst(x,y,4);
			u_y=dst(x,y,8)+dst(x,y,1)+dst(x,y,2);
			rho=dst(x,y,5)+dst(x,y,6)+dst(x,y,7)+dst(x,y,0)-dst(x,y,2)+u_x+u_y;
			u_x=(u_x-dst(x,y,6)-dst(x,y,7)-dst(x,y,8));//rho;
			u_y=(u_y-dst(x,y,4)-dst(x,y,5)-dst(x,y,6));//rho;
			*/	
			//std::cout << rho << std::endl;
			//std::cout<< "velo2" << x <<" "<< y << " "<< u_x << " " << u_y << std::endl;

			
			if(k!=0 && k%100==0){
				// call print-function
				u_x=dst(x,y,2)+dst(x,y,3)+dst(x,y,4);
				u_y=dst(x,y,8)+dst(x,y,1)+dst(x,y,2);
				rho=dst(x,y,5)+dst(x,y,6)+dst(x,y,7)+dst(x,y,0)-dst(x,y,2)+u_x+u_y;
				u_x=(u_x-dst(x,y,6)-dst(x,y,7)-dst(x,y,8))/rho;
				u_y=(u_y-dst(x,y,4)-dst(x,y,5)-dst(x,y,6))/rho;
				velocity(x,y,0)=u_x;
				velocity(x,y,1)=u_y;
				density(x,y,0)=rho;	

				printmethod(velocity,density,sizex,sizey,k);			
			
			}

			
			//std::cout << rho << std::endl;

			}
		}
	//std::cout << rho_sum << std::endl;
	//std::cout << "next" << std::endl;

	src.swap(dst);
	
	

	}
	gettimeofday(&t, NULL);
	
	std::cout << "Wall clock time of execution: " << (((int64_t) (t.tv_sec - t0.tv_sec) * (int64_t)1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3) << " ms" << std::endl;


	return 0;

}
