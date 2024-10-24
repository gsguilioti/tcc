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

#define K 20
#define MAX_LK_ITER 10

double euclideanDistance(std::pair<double, double> a, std::pair<double, double> b) 
{
    return std::sqrt(std::pow(a.first - b.first, 2) + std::pow(a.second - b.second, 2));
}

int getRandomCity(const std::map<int, std::pair<double, double>>& cities) 
{
    std::vector<int> cityKeys;
    for (const auto& [key, _] : cities) 
        cityKeys.push_back(key);

    int randomIndex = rand() % cityKeys.size();
    return cityKeys[randomIndex];
}

std::vector<int> generateTour(const std::map<int, std::pair<double, double>>& cities, const std::map<int, std::vector<int>>& candidateList)
{
    srand(static_cast<unsigned int>(time(0)));
    int currentCity = getRandomCity(cities);

    std::vector<int> tour;
    std::set<int> visited;
    tour.push_back(currentCity);
    visited.insert(currentCity);

    for (const auto& city : cities)
    {
        const auto& neighbors = candidateList.at(currentCity);
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

        if (nextCity == -1) 
        {
            for (const auto& [key, _] : cities) 
            {
                if (visited.find(key) == visited.end()) 
                {
                    nextCity = key;
                    break;
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
            break;
    }

    return tour;
}

std::map<int, std::vector<int>> generateCandidateList(const std::map<int, std::pair<double, double>>& cities, int k) 
{
    std::map<int, std::vector<int>> candidateList;
    for (const auto& [key, value] : cities) 
    {
        std::vector<std::pair<int, double>> distances;
        for (const auto& [otherKey, otherValue] : cities) 
        {
            if (otherKey == key) continue;
            distances.emplace_back(otherKey, euclideanDistance(value, otherValue));
        }

        std::sort(distances.begin(), distances.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) 
        {
            return a.second < b.second;
        });

        std::vector<int> neighbors;
        for (int i = 0; i < std::min(k, static_cast<int>(distances.size())); ++i) 
        {
            neighbors.push_back(distances[i].first);
        }
        candidateList[key] = neighbors;
    }
    return candidateList;
}

double calculateTotalDistance(const std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i)
        totalDistance += euclideanDistance(cities.at(tour[i]), cities.at(tour[i+1]));

    totalDistance += euclideanDistance(cities.at(tour.back()), cities.at(tour.front()));
    return totalDistance;
}

int findPosition(const std::vector<int>& tour, int city) 
{
    auto it = std::find(tour.begin(), tour.end(), city);
    if (it != tour.end()) 
        return std::distance(tour.begin(), it);
        
    return -1;
}

void reverseSegment(std::vector<int>& tour, int start, int end) 
{
    int n = tour.size();
    start = start % n;
    end = end % n;

    if (start < end) 
    {
        std::reverse(tour.begin() + start, tour.begin() + end + 1);
    } 
    else if (start > end) 
    {
        std::reverse(tour.begin() + start, tour.end());
        std::reverse(tour.begin(), tour.begin() + end + 1);
        std::reverse(tour.begin(), tour.end());
    }
}

bool isValidTour(const std::vector<int>& tour, int numCities) 
{
    if (tour.size() != numCities) return false;
    std::set<int> uniqueCities(tour.begin(), tour.end());
    return uniqueCities.size() == numCities;
}

void linKernighan(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities, const std::map<int, std::vector<int>>& candidateList) 
{
    bool improved = true;
    int iterations = 0;

    while (improved && iterations < MAX_LK_ITER) 
    {
        improved = false;
        double bestGain = 0.0;
        int best_i = -1, best_j = -1;

        for (size_t i = 0; i < tour.size(); ++i) 
        {
            int A = tour[i];

            for (int B : candidateList.at(A)) 
            {
                int posB = findPosition(tour, B);
                int posA_next = (i + 1) % tour.size();
                int C = tour[posA_next];

                double gain = euclideanDistance(cities.at(A), cities.at(C)) - euclideanDistance(cities.at(A), cities.at(B));

                if (gain > bestGain)
                {
                    std::vector<int> newTour = tour;
                    reverseSegment(newTour, posA_next, posB);

                    if (isValidTour(newTour, cities.size())) 
                    {
                        double newDistance = calculateTotalDistance(newTour, cities);
                        double currentDistance = calculateTotalDistance(tour, cities);
                        double totalGain = currentDistance - newDistance;

                        if (totalGain > bestGain) {
                            bestGain = totalGain;
                            best_i = i;
                            best_j = posB;
                            improved = true;
                            tour = newTour;
                            break;
                        }
                    }
                }
            }
            if (improved) break;
        }

        if (improved) {
            iterations++;
            std::cout << "Iteration " << iterations << ", Gain: " << bestGain << "\n";
        }
    }

    if (iterations >= MAX_LK_ITER) {
        std::cout << "Reached maximum LK iterations.\n";
    } else {
        std::cout << "LK Optimization completed in " << iterations << " iterations.\n";
    }
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

    std::map<int, std::vector<int>> candidateList = generateCandidateList(cities, K);
    std::vector<int> tour = generateTour(cities, candidateList);

    double initialDistance = calculateTotalDistance(tour, cities);
    std::cout << "Initial Tour Distance: " << initialDistance << "\n";

    linKernighan(tour, cities, candidateList);

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
