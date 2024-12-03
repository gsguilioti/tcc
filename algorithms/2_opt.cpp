#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

double euclideanDistance(std::pair<double, double> a, std::pair<double, double> b) 
{
    return std::sqrt(std::pow(a.first - b.first, 2) + std::pow(a.second - b.second, 2));
}

double calculateTotalDistance(const std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
        totalDistance += euclideanDistance(cities.at(tour[i]), cities.at(tour[i+1]));

    totalDistance += euclideanDistance(cities.at(tour.back()), cities.at(tour.front()));
    return totalDistance;
}

void swap_edges(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities, int i, int j) 
{
    i += 1;
    while (i < j)
    {
        std::swap(tour[i], tour[j]);
        i++;
        j--;
    }
}

void two_opt(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities)
{
    float currentLength = calculateTotalDistance(tour, cities);
    bool improved {true};
    while(improved)
    {
        improved = false;

        for(int i = 0; i < tour.size() -1; ++i)
        {
            for (int j = i + 2; j < tour.size(); ++j)
            {
                std::vector<int> auxTour(tour);
                swap_edges(auxTour, cities, i, j);

                float auxLength = calculateTotalDistance(auxTour, cities);

                if(auxLength < currentLength)
                {
                    currentLength = auxLength;
                    tour = auxTour;
                    improved = true;
                }
            }
        }
    }
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

    std::string toursFilepath = "../tours/";
    toursFilepath.append(argv[1]);
    toursFilepath.append(".txt");
    std::ifstream toursFile(toursFilepath);
    if (!toursFile.is_open()) {
        std::cerr << "Error: Could not open tours file " << toursFilepath << "\n";
        return 1;
    }

    std::string resultsFilepath = "../results/";
    resultsFilepath.append(argv[1]);
    resultsFilepath.append(".result");
    resultsFilepath.append(".txt");
    std::ofstream resultsFile(resultsFilepath, std::ios::app);
    if (!resultsFile.is_open()) {
        std::cerr << "Error: Could not open results file " << resultsFilepath << "\n";
        return 1;
    }

    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(resultsFile.rdbuf());

    std::cout << "2OPT --------------------------------------------------------------------------------------  \n";

    int tourIndex = 1;
    std::string line;
    while (std::getline(toursFile, line)) {
        std::istringstream iss(line);
        std::vector<int> tour;
        int city;
        while (iss >> city) {
            tour.push_back(city);
        }

        if (tour.empty()) {
            std::cerr << "Error: Empty tour on line " << tourIndex << "\n";
            continue;
        }

        std::cout << "Processing Tour #" << tourIndex << "\n";

        std::vector<int> currentTour = tour;
        double initialDistance = calculateTotalDistance(currentTour, cities);
        std::cout << "Initial Distance: " << initialDistance << "\n";

        auto start = std::chrono::high_resolution_clock::now();
        two_opt(currentTour, cities);
        auto end = std::chrono::high_resolution_clock::now();

        double optimizedDistance = calculateTotalDistance(currentTour, cities);
        std::cout << "Optimized Distance: " << optimizedDistance << "\n";

        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Execution Time: " << elapsed.count() << " seconds\n";

        std::cout << "----------------------------------------------- \n";
        
        ++tourIndex;
    }

    toursFile.close();
     resultsFile.close();
    std::cout.rdbuf(originalCoutBuffer);

    return 0;
}
