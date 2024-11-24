#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <thread>
#include <mutex> 

std::mutex tourMutex;

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
    std::vector<int> copyTour(tour);
    double totalDistance = 0.0;
    for (size_t i = 0; i < copyTour.size() - 1; ++i) 
    {
        if (cities.find(copyTour[i]) == cities.end() || cities.find(copyTour[i + 1]) == cities.end()) 
            continue;

        totalDistance += euclideanDistance(cities.at(copyTour[i]), cities.at(copyTour[i + 1]));
    }

    if (cities.find(copyTour.back()) == cities.end() || cities.find(copyTour.front()) != cities.end())
         totalDistance += euclideanDistance(cities.at(copyTour.back()), cities.at(copyTour.front()));

    return totalDistance;
}

void applyReconnection(std::vector<int>& tour, int i, int j, int k, int option)
{
    switch (option)
    {
        case 1: 
            std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
            break;
        case 2:
            std::reverse(tour.begin() + j + 1, tour.begin() + k + 1);
            break;
        case 3:
            std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
            std::reverse(tour.begin() + j + 1, tour.begin() + k + 1);
            break;
        case 4:
            std::reverse(tour.begin() + k + 1, tour.end());
            std::reverse(tour.begin(), tour.begin() + i + 1);
            break;
        case 5:
            std::reverse(tour.begin() + i + 1, tour.end());
            break;
        default:
            break;
    }
}

void best_option_search(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities, 
                           int i, int j, int k, bool& improved)
{
    double bestGain = 0.0;
    int bestOption = -1;

    for (int option = 1; option <= 5; ++option) 
    {
        std::vector<int> candidateTour = tour;
        applyReconnection(candidateTour, i, j, k, option);
        double gain = calculateTotalDistance(tour, cities) - calculateTotalDistance(candidateTour, cities);

        if (gain > bestGain) 
        {
            bestGain = gain;
            bestOption = option;
        }
    }

    if (bestOption != -1) 
    {
        std::lock_guard<std::mutex> lock(tourMutex);
        applyReconnection(tour, i, j, k, bestOption);
        improved = true;
    }
}

void three_opt(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    bool improved = true;

    while (improved) 
    {
        improved = false;
        std::vector<std::thread> threads;

        for (int i = 0; i < tour.size() - 2; ++i) 
        {
            for (int j = i + 1; j < tour.size() - 1; ++j) 
            {
                for (int k = j + 1; k < tour.size(); ++k) 
                    threads.push_back(std::thread(best_option_search, std::ref(tour), std::cref(cities), i, j, k, std::ref(improved)));
            }
        }

        for (auto& t : threads)
        {
            t.join();
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
        std::istringstream iss(node);
        int key;
        double value1, value2;
        if(iss >> key >> value1 >> value2)
            cities[key] = std::make_pair(value1, value2);
    }

    instance.close();

    std::vector<int> tour = generateTour(cities.size());

    double initialDistance = calculateTotalDistance(tour, cities);
    std::cout << "Initial Tour Distance: " << initialDistance << "\n";

    auto start = std::chrono::high_resolution_clock::now();

    three_opt(tour, cities);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution Time: " << elapsed.count() << " seconds\n";

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
