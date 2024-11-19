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

void reverseSegment(std::vector<int>& tour, int start, int end) 
{
    while (start < end) {
        std::swap(tour[start], tour[end]);
        start++;
        end--;
    }
}

void swap_edges_3opt(std::vector<int>& tour, int i, int j, int k) 
{
    reverseSegment(tour, i + 1, j);
    reverseSegment(tour, j + 1, k);
    reverseSegment(tour, i + 1, k);
}

void three_opt(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities, int start, int end, double& minDistance, std::mutex& mutex) 
{
    double currentLength = calculateTotalDistance(tour, cities);
    bool improved {true};

    while (improved) 
    {
        improved = false;

        for (int i = start; i < end - 2; ++i) 
        {
            for (int j = i + 1; j < end - 1; ++j) 
            {
                for (int k = j + 1; k < end; ++k) 
                {
                    std::vector<int> auxTour(tour);
                    swap_edges_3opt(auxTour, i, j, k);
                    double auxLength = calculateTotalDistance(auxTour, cities);

                    if (auxLength < currentLength) {
                        std::lock_guard<std::mutex> lock(mutex);
                        if (auxLength < minDistance) 
                        {
                            minDistance = auxLength;
                            tour = auxTour;
                        }
                        currentLength = auxLength;
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
        if (node.find("NODE_COORD_SECTION"))
            break;
    }

    std::map<int, std::pair<double, double>> cities;
    while(std::getline(instance, node)) 
    {
        std::istringstream iss(node);
        int key;
        double value1, value2;
        if (iss >> key >> value1 >> value2)
            cities[key] = std::make_pair(value1, value2);
    }

    instance.close();

    std::vector<int> tour = generateTour(cities.size());

    double initialDistance = calculateTotalDistance(tour, cities);
    std::cout << "Initial Tour Distance: " << initialDistance << "\n";

    auto start = std::chrono::high_resolution_clock::now();

    int numThreads = std::thread::hardware_concurrency();
    int segmentSize = tour.size() / numThreads;
    double minDistance = calculateTotalDistance(tour, cities);
    std::mutex mutex;
    std::vector<std::thread> threads;

    for (int t = 0; t < numThreads; ++t) 
    {
        int start = t * segmentSize;
        int end = (t == numThreads - 1) ? tour.size() : (t + 1) * segmentSize;
        threads.emplace_back(three_opt, std::ref(tour), std::ref(cities), start, end, std::ref(minDistance), std::ref(mutex));
    }

    for (auto& thread : threads) 
        thread.join();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution Time: " << elapsed.count() << " seconds\n";

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
