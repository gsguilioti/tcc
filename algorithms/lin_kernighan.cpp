#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
#include <limits>
#include <fstream>
#include <sstream>

double euclideanDistance(std::pair<double, double> a, std::pair<double, double> b) 
{
    return std::sqrt(std::pow(a.first - b.first, 2) + std::pow(a.second - b.second, 2));
}

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

double calculateTotalDistance(const std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
        totalDistance += euclideanDistance(cities.at(tour[i]), cities.at(tour[i+1]));

    totalDistance += euclideanDistance(cities.at(tour.back()), cities.at(tour.front()));
    return totalDistance;
}

// Perform a k-opt move based on selected edges (x and y)
void performKOpt(std::vector<int>& tour, const std::vector<int>& selectedXs, const std::vector<int>& selectedYs) 
{
    for (size_t i = 0; i < selectedYs.size(); ++i) 
    {
        int t1 = selectedXs[i];
        int t2 = selectedYs[i];
        std::swap(tour[t1], tour[t2]);
    }
}

double calculateGain(const std::map<int, std::pair<double, double>>& cities, int t1, int t2, int t3) 
{
    auto cityT1 = cities.find(t1);
    auto cityT2 = cities.find(t2);
    auto cityT3 = cities.find(t3);

    if (cityT1 == cities.end() || cityT2 == cities.end() || cityT3 == cities.end())
        return 0;

    return euclideanDistance(cityT1->second, cityT2->second) - euclideanDistance(cityT2->second, cityT3->second);
}

void linKernighan(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    bool improvement = true;
    double bestDistance = calculateTotalDistance(tour, cities);
    std::vector<int> bestTour = tour;

    while (improvement) 
    {
        improvement = false;
        double maxGain = 0.0;

        for (size_t t1 = 0; t1 < tour.size(); ++t1) 
        {
            std::vector<int> selectedXs;
            std::vector<int> selectedYs;
            int i = 1;

            while (true) 
            {
                int t2 = (t1 + i) % tour.size();
                selectedXs.push_back(t1);
                selectedXs.push_back(t2);

                int t3 = (t2 + 1) % tour.size();
                double gain = calculateGain(cities, t1, t2, t3);

                if (gain > maxGain) 
                {
                    maxGain = gain;
                    selectedYs.push_back(t2);
                    selectedYs.push_back(t3);

                    std::vector<int> tempTour = tour;
                    performKOpt(tempTour, selectedXs, selectedYs);
                    double newDistance = calculateTotalDistance(tempTour, cities);

                    if (newDistance < bestDistance) 
                    {
                        bestTour = tempTour;
                        bestDistance = newDistance;
                        improvement = true;
                    }
                }

                if (gain <= 0) 
                    break;

                i++;
            }

            if (improvement) 
            {
                tour = bestTour;
                break;
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

    std::vector<int> tour = generateTour(cities.size());

    double initialDistance = calculateTotalDistance(tour, cities);
    std::cout << "Initial Tour Distance: " << initialDistance << "\n";

    linKernighan(tour, cities);

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
