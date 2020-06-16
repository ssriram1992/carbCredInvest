#include "carbCredInv.h"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
	if (argc != 4){
		std::cerr << "Usage: readWrite <location_fileName> <solution_fileName> <output_fileName>\n";
		return 0;
	}
	std::string locName(argv[1]);
	std::string solName(argv[2]);
	std::string outName(argv[3]);
	std :: cout << solName << '\n';
	return 0;
}
