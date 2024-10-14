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

bool twoOptSwap(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    bool improved = false;
    double bestGain = 0.0;
    int best_i = 0, best_k = 0;

    for (size_t i = 1; i < tour.size() - 1; ++i)
    {
        for (size_t k = i + 1; k < tour.size(); ++k) 
        {
            int A = tour[i - 1];
            int B = tour[i];
            int C = tour[k];
            int D = tour[(k + 1) % tour.size()];

            double currentDist = euclideanDistance(cities.at(A), cities.at(B)) + euclideanDistance(cities.at(C), cities.at(D));
            double newDist = euclideanDistance(cities.at(A), cities.at(C)) + euclideanDistance(cities.at(B), cities.at(D));
            double gain = currentDist - newDist;

            if (gain > bestGain) 
            {
                bestGain = gain;
                best_i = i;
                best_k = k;
                improved = true;
            }
        }
    }

    if (improved) 
        reverseSegment(tour, best_i, best_k);

    return improved;
}

bool threeOptSwap(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    bool improved = false;
    double bestGain = 0.0;
    int best_i = 0, best_j = 0, best_k = 0;

    for (size_t i = 1; i < tour.size() - 2; ++i) 
    {
        for (size_t j = i + 1; j < tour.size() - 1; ++j) 
        {
            for (size_t k = j + 1; k < tour.size(); ++k) 
            {
                int A = tour[i - 1];
                int B = tour[i];
                int C = tour[j];
                int D = tour[(j + 1) % tour.size()];
                int E = tour[k];
                int F = tour[(k + 1) % tour.size()];

                double currentDist = euclideanDistance(cities.at(A), cities.at(B)) +
                                     euclideanDistance(cities.at(C), cities.at(D)) +
                                     euclideanDistance(cities.at(E), cities.at(F));

                std::vector<double> newDists(4);

                newDists[0] = euclideanDistance(cities.at(A), cities.at(C)) +
                              euclideanDistance(cities.at(B), cities.at(D)) +
                              euclideanDistance(cities.at(E), cities.at(F));

                newDists[1] = euclideanDistance(cities.at(A), cities.at(B)) +
                              euclideanDistance(cities.at(C), cities.at(E)) +
                              euclideanDistance(cities.at(D), cities.at(F));

                newDists[2] = euclideanDistance(cities.at(A), cities.at(C)) +
                              euclideanDistance(cities.at(B), cities.at(E)) +
                              euclideanDistance(cities.at(D), cities.at(F));

                newDists[3] = euclideanDistance(cities.at(A), cities.at(B)) +
                              euclideanDistance(cities.at(C), cities.at(F)) +
                              euclideanDistance(cities.at(D), cities.at(E));

                auto minIter = std::min_element(newDists.begin(), newDists.end());
                double newDist = *minIter;
                double gain = currentDist - newDist;

                if (gain > bestGain) 
                {
                    bestGain = gain;
                    best_i = i;
                    best_j = j;
                    best_k = k;
                    improved = true;

                    switch (std::distance(newDists.begin(), minIter)) 
                    {
                        case 0:
                            reverseSegment(tour, best_i, best_j);
                            break;
                        case 1:
                            reverseSegment(tour, best_j + 1, best_k);
                            break;
                        case 2:
                            reverseSegment(tour, best_i, best_j);
                            reverseSegment(tour, best_j + 1, best_k);
                            break;
                        case 3:
                            reverseSegment(tour, best_i + 1, best_j);
                            reverseSegment(tour, best_j + 1, best_k);
                            break;
                    }
                }
            }
        }
    }

    return improved;
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

    bool improved = true;
    int iteration = 0;
    while (improved && iteration < MAX_LK_ITER) 
    {
        improved = false;

        bool twoOptImprovement = twoOptSwap(tour, cities);
        if (twoOptImprovement) 
        {
            improved = true;
            std::cout << "2-opt improvement found on iteration " << iteration << "\n";
        }

        bool threeOptImprovement = threeOptSwap(tour, cities);
        if (threeOptImprovement) 
        {
            improved = true;
            std::cout << "3-opt improvement found on iteration " << iteration << "\n";
        }

        if (improved)
        {
            iteration++;
            double currentDistance = calculateTotalDistance(tour, cities);
            std::cout << "Iteration " << iteration << ", Tour Distance: " << currentDistance << "\n";
        }
    }

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
