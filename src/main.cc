#include <iostream>
#include "FileReader.hh"
#include <string>
#include "Grid.hh"
#include <fstream>
#include <sys/time.h>
#include <sys/stat.h>

enum cellflag {BOUNDARY,ACCELERATION,FLUID};


static void printmethod(lbm::V_Field &velocity, lbm::D_Field &density, size_t sizex,size_t sizey, size_t k, std::string vtk_file){

	size_t start=vtk_file.find(".");	
	char* filestring= new char [start];	
	std::size_t length = vtk_file.copy(filestring,start,0);
	filestring[length]='\0';
	char* fileextension = new char[vtk_file.size()-start+1];
	length = vtk_file.copy(fileextension,vtk_file.size()-start+1,start);
	fileextension[length]='\0';

	std::ofstream vtkFile( filestring + std::to_string(k)  + fileextension);

	if(vtkFile.is_open()){
		vtkFile << "# vtk DataFile Version 4.0 \nSiwiRVisFile \nASCII \nDATASET STRUCTURED_POINTS \nDIMENSIONS " << sizex << " ";
		vtkFile << sizey << " " << 1 << "\n" << "ORIGIN 0 0 0"<< "\n";
		vtkFile << "SPACING 1 1 1\nPOINT_DATA " << sizex*sizey << "\n" << "\n" ;
		vtkFile << "SCALARS flags double 1 \nLOOKUP_TABLE default" << "\n" ;

		for(size_t y=2;y<sizey+2;++y){
			for(size_t x=2;x<sizex+2;++x){
				vtkFile << 1 << "\n";
			}

		}

		vtkFile << "\nSCALARS density double 1 \nLOOKUP_TABLE default \n" ;

		for(size_t y=2;y<sizey+2;++y){
			for(size_t x=2;x<sizex+2;++x){
				vtkFile << density(x,y,0)<< "\n";
			}

		}

		vtkFile << "\nVECTORS velocity double\n" ;

		for(size_t y=2;y<sizey+2;++y){
			for(size_t x=2;x<sizex+2;++x){
				vtkFile << velocity(x,y,0) << " " << velocity(x,y,1) << " " << 0 <<"\n";
			}

		}
		vtkFile.close();
	} else {
		std::cout << "Please provide existing Directory for Output in the Parameter-File" << std::endl;
		std::exit(EXIT_FAILURE);
	}

}

static void initializer(lbm::PDF_Field &field,size_t sizex, size_t sizey){
	// Initialize Field
	for(size_t y=0;y<sizey;++y){
		for(size_t x=0;x<sizex;++x){
			field(x,y,0)=4.0/9.0;
			field(x,y,1)=1.0/9.0;
			field(x,y,2)=1.0/36.0;
			field(x,y,3)=1.0/9.0;
			field(x,y,4)=1.0/36.0;
			field(x,y,5)=1.0/9.0;
			field(x,y,6)=1.0/36.0;
			field(x,y,7)=1.0/9.0;
			field(x,y,8)=1.0/36.0;
		}
	}
}

static void flags_initializer(lbm::Flags &flags,size_t sizex, size_t sizey){
	// Initialize Flag-Field
	for(size_t y=1;y<sizey+3;++y){
		for(size_t x=1;x<sizex+3;++x){
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

	const double omega(reader.getParameter<double>("omega"));
	const std::string vtk_file(reader.getParameter<std::string>("vtk_file"));
	const size_t sizex(reader.getParameter<size_t>("sizex"));
	const size_t sizey(reader.getParameter<size_t>("sizey"));
	const size_t timesteps(reader.getParameter<size_t>("timesteps"));
	const size_t vtk_step(reader.getParameter<size_t>("vtk_step"));

	//const size_t sizex(100);
	//const size_t sizey(100);
	//const double timesteps(10000);
	//const size_t vtk_step(100);

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

	for(size_t k=0;k<timesteps;++k){
		for(size_t y=1;y<sizey+3;++y){
			for(size_t x=1;x<sizex+3;++x){

				// stream-step on fluid-cells and boundaries
				dst(x,y,0)=src(x,y,0);
				dst(x,y,1)=src(x,y-1,1);
				dst(x,y,2)=src(x-1,y-1,2);
				dst(x,y,3)=src(x-1,y,3);
				dst(x,y,4)=src(x-1,y+1,4);
				dst(x,y,5)=src(x,y+1,5);
				dst(x,y,6)=src(x+1,y+1,6);
				dst(x,y,7)=src(x+1,y,7);
				dst(x,y,8)=src(x+1,y-1,8);

				//handle boundaries
				if(flags(x,y,0)==BOUNDARY){

					std::swap(dst(x,y,1),dst(x,y,5));
					std::swap(dst(x,y,2),dst(x,y,6));
					std::swap(dst(x,y,3),dst(x,y,7));
					std::swap(dst(x,y,4),dst(x,y,8));

					continue;

				}

				double rho;
				double u_x;
				double u_y;

				//handle acceleration cell 
				if(flags(x,y,0)==ACCELERATION){
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

					continue;

				}

				//collide regular fluid cell

				u_x=dst(x,y,2)+dst(x,y,3)+dst(x,y,4);
				u_y=dst(x,y,8)+dst(x,y,1)+dst(x,y,2);
				rho=dst(x,y,5)+dst(x,y,6)+dst(x,y,7)+dst(x,y,0)-dst(x,y,2)+u_x+u_y;
				u_x=(u_x-dst(x,y,6)-dst(x,y,7)-dst(x,y,8));
				u_y=(u_y-dst(x,y,4)-dst(x,y,5)-dst(x,y,6));

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

				if(k!=0 && k%vtk_step==0 && vtk_step!=0 ){
					// keep data for print-function
					velocity(x,y,0)=u_x;
					velocity(x,y,1)=u_y;
					density(x,y,0)=rho;				

				}
			}
		}
		if(k!=0 && k%vtk_step==0 && vtk_step!=0 ){
			// call print-function
			printmethod(velocity,density,sizex,sizey,k,vtk_file);			
		}
		src.swap(dst);
	}
	gettimeofday(&t, NULL);

	std::cout << "MLUps: " << (sizex*sizey*timesteps)/(((int64_t) (t.tv_sec - t0.tv_sec) * (int64_t)1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3)*1e-3 << std::endl;


	return 0;

}
