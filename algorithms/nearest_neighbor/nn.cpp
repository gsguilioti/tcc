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
#include <set>
#include <cstdlib>

#define K 3

double euclideanDistance(std::pair<double, double> a, std::pair<double, double> b) {
    return std::sqrt(std::pow(a.first - b.first, 2) + std::pow(a.second - b.second, 2));
}

std::vector<int> knn(const std::map<int, std::pair<double, double>>& cities, const std::set<int>& visited, int index, int k)
{
    auto city = cities.find(index)->second;

    std::vector<std::pair<int, double>> distances;
    for(const auto& [key, value] : cities)
    {
        if(key == index)
            continue;

        distances.emplace_back(key, euclideanDistance(city, value));
    }

    std::sort(distances.begin(), distances.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    std::vector<int> neighbors;
    for (const auto& item : distances) 
    {
        if(visited.find(item.first) == visited.end())
            neighbors.push_back(item.first);
    }

    if (neighbors.size() > K)
        neighbors.resize(K);

    return neighbors;
}

std::vector<int> generateTour(const std::map<int, std::pair<double, double>>& cities, int k)
{
    srand(static_cast<unsigned int>(time(0)));
    int currentCity = rand() % cities.size() + 1;

    std::vector<int> tour;
    std::set<int> visited;
    tour.push_back(currentCity);
    visited.insert(currentCity);

    for (const auto& city : cities)
    {
        std::vector<int> neighbors = knn(cities, visited, city.first, K);

        int nextCity = -1;
        double minDistance = std::numeric_limits<double>::max();
        
        for (int neighbor : neighbors) 
        {
            if (visited.find(neighbor) == visited.end()) 
            {
                double distance = euclideanDistance(cities.find(currentCity)->second, cities.find(neighbor)->second);
                if (distance < minDistance) 
                {
                    minDistance = distance;
                    nextCity = neighbor;
                }
            }
        }

        if (nextCity != -1) 
        {
            tour.push_back(nextCity);
            visited.insert(nextCity);
            currentCity = nextCity;
        } 
        else
            continue;
    }

    return tour;
}

int main(int argc, char* argv[])
{
    std::string filename {"../../instances/"};
    std::ifstream instance(filename.append(argv[1]));

    std::string node;
    while(std::getline(instance, node))
    {
        if(node == "NODE_COORD_SECTION")
            break;
    }

    std::map<int, std::pair<double, double>> cities;
    while(std::getline(instance, node))
    {
        std::istringstream  iss(node);

        int key;
        double value1, value2;
        if(iss >> key >> value1 >> value2)
            cities[key] = std::make_pair(value1, value2);
    }

    instance.close();

    std::vector<int> tour = generateTour(cities, K);

    double totalDistance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) 
    {
        int cityA = tour[i];
        int cityB = tour[i + 1];
        totalDistance += euclideanDistance(cities.at(cityA), cities.at(cityB));
    }
    totalDistance += euclideanDistance(cities.at(tour.back()), cities.at(tour.front()));
    std::cout << totalDistance << "\n";

    return 0;
}
