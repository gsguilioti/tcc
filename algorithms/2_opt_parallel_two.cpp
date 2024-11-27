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
#include <thread>
#include <mutex>

#define NUM_THREADS 8

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

void process_j_chunk(int i, int startIdx, int endIdx, std::vector<int>& tour, 
                     const std::map<int, std::pair<double, double>>& cities, 
                     double& currentLength, bool& improved, std::mutex& mtx) 
{
    for (int j = startIdx; j < endIdx; ++j) 
    {
        std::vector<int> auxTour(tour);
        swap_edges(auxTour, cities, i, j);

        float auxLength = calculateTotalDistance(auxTour, cities);

        if (auxLength < currentLength) 
        {
            std::lock_guard<std::mutex> lock(mtx);
            if (auxLength < currentLength) 
            { 
                currentLength = auxLength;
                tour = auxTour;
                improved = true;
            }
        }
    }
}

void parallel_two_opt(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities) 
{
    double currentLength = calculateTotalDistance(tour, cities);
    bool improved = true;
    std::mutex mtx;

    while (improved) 
    {
        improved = false;

        for (int i = 0; i < tour.size() - 1; ++i) 
        {
            std::vector<std::thread> threads;

            int chunkSize = (tour.size() - (i + 2)) / NUM_THREADS;
            for (int t = 0; t < NUM_THREADS; ++t) 
            {
                int startIdx = (i + 2) + t * chunkSize;
                int endIdx = (t == NUM_THREADS - 1) ? (tour.size()) : startIdx + chunkSize;

                threads.push_back(std::thread(process_j_chunk, i, startIdx, endIdx, 
                                              std::ref(tour), std::cref(cities), 
                                              std::ref(currentLength), std::ref(improved), std::ref(mtx)));
            }

            for (auto& th : threads) 
            {
                th.join();
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

    auto start = std::chrono::high_resolution_clock::now();

    parallel_two_opt(tour, cities);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution Time: " << elapsed.count() << " seconds\n";

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
