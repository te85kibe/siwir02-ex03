#include <iostream>
#include "FileReader.hh"
#include <string>
#include "Grid.hh"
int main(int argc, char **args)
{
	if(argc != 2)
	{
		std::cout << "Usage: ./lbm <params.dat>" << std::endl;
		return 0;
	}	

	FileReader reader(args[1]);
	std::cout << reader.getParameter<size_t>("vtk_step") << std::endl;
	std::cout << reader.getParameter<std::string>("vtk_file") << std::endl;
	std::cout << reader.getParameter<double>("timesteps") << std::endl;
	std::cout << reader.getParameter<double>("omega") << std::endl;
	
	const size_t timesteps(reader.getParameter<size_t>("timesteps"));
	std::cout << timesteps << std::endl;

	Grid<double,3> grid(9,9);
	
	return 0;

}
#if 0
if(argc != 2){
	std::cout << "Usage: ./mgsolve <number_of_levels>" << std::endl;
	return 0;
}

int l;

//time
//	siwir::Timer ti;
//	double time;

std::istringstream iss(args[1]);
if(!(iss >> l)){
	std::cerr << "Could not parse number of level argument: " << args[1] << std::endl;
	return 1;
}
iss.str("");
iss.clear();

Smoother smoother;
MGSolver solver(l, smoother);
solver.initialize_seminar();
solver.saveToFile("init.dat");
std::cout<<"Your Alias: "<<"Thomas+Nils"<<std::endl;

struct timeval t0, t;

gettimeofday(&t0, NULL);

solver.v_cycle(2, 1, 4);

gettimeofday(&t, NULL);
std::cout << "Wall clock time of MG execution: " << ((int64_t) (t.tv_sec - t0.tv_sec) * (int64_t)1000000 + 					(int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3 << " ms" << std::endl;

//	time = ti.elapsed();
//	std::cout << "Time: " << "\t" << time << std::endl;

solver.saveToFile("solution.dat");

return 0;
#endif


