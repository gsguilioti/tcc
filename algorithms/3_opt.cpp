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

void reverseSegment(std::vector<int>& tour, int start, int end) {
    while (start < end) {
        std::swap(tour[start], tour[end]);
        start++;
        end--;
    }
}

void swap_edges_3opt(std::vector<int>& tour, int i, int j, int k) {
    reverseSegment(tour, i + 1, j);
    reverseSegment(tour, j + 1, k);
    reverseSegment(tour, i + 1, k);
}

void three_opt(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) {
    float currentLength = calculateTotalDistance(tour, cities);
    bool improved {true};

    while (improved) {
        improved = false;

        for (int i = 0; i < tour.size() - 2; ++i) {
            for (int j = i + 1; j < tour.size() - 1; ++j) {
                for (int k = j + 1; k < tour.size(); ++k) {
                    std::vector<int> auxTour(tour);

                    swap_edges_3opt(auxTour, i, j, k);

                    float auxLength = calculateTotalDistance(auxTour, cities);

                    if (auxLength < currentLength) {
                        currentLength = auxLength;
                        tour = auxTour;
                        improved = true;
                    }
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

    three_opt(tour, cities);

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
