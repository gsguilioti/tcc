#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>

std::vector<int> generateTour(int size)
{
    std::vector<int> tour(size);
    for (int i = 0; i < size; ++i) 
        tour[i] = i + 1;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(tour.begin(), tour.end(), rng);

    return tour;
}

int main(int argc, char* argv[])
{
    std::string filename {"../instances/"};
    std::ifstream instance(filename.append(argv[1]));

    std::string node;
    while(std::getline(instance, node))
    {
        if(node.find("NODE_COORD_SECTION"))
            break;
    }

    std::map<int, std::pair<double, double>> cities;
    while(std::getline(instance, node))
    {
         if (node == "EOF")
            break;
        
        std::istringstream  iss(node);

        int key;
        double value1, value2;
        if(iss >> key >> value1 >> value2)
            cities[key] = std::make_pair(value1, value2);
    }

    instance.close();

    std::string outputFilename = "../tours/";
    outputFilename.append(argv[1]);
    outputFilename.append(".txt");
    std::ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not create file " << outputFilename << "\n";
        return 1;
    }

    for(int i = 0; i < 5; ++i)
    {
        std::vector<int> tour = generateTour(cities.size());
        for (size_t j = 0; j < tour.size(); ++j) 
        {
            outputFile << tour[j];
            if (j < tour.size() - 1)
                outputFile << " ";
        }
        outputFile << "\n";
    }

    outputFile.close();

    return 0;
}
