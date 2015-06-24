#include <iostream>
#include <fstream>
#include <sstream>
//template <typename T>
class FileReader{
public:
	FileReader(const std::string& file);
template<typename T>
	T getParameter(const std::string& file);

private:
    std::string parameters;
};

FileReader::FileReader(const std::string& file){
	std::ifstream in;
	in.open(file);
	std::string buffer;

	while (in.good()) {
		getline(in,buffer);
		parameters.append(buffer + " ");
	}	
	in.close();
	//std::cout << parameters << std::endl;
}

template<typename T>
	T FileReader::getParameter(const std::string& key){
	// find parametervalue	
	size_t start=parameters.find(key);
	start=parameters.find(' ',start);
	size_t end=parameters.find(' ',start+1);
	char* parameterString= new char [end-start+1];
	std::size_t length = parameters.copy(parameterString,end-start,start+1);
	parameterString[length]='\0';
	// convert string to type
	std::stringstream ss;
	ss << parameterString;
	T returnvalue;
	ss >> returnvalue;
	return returnvalue;
    }
