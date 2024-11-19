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
#include <thread>
#include <mutex>
#include <chrono>

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

void swap_edges(std::vector<int>& tour, int i, int j) 
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

        for(int i = 0; i < tour.size() - 1; ++i)
        {
            for (int j = i + 2; j < tour.size(); ++j)
            {
                std::vector<int> auxTour(tour);
                swap_edges(auxTour, i, j);

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

void optimizeSegment(std::vector<int>& tour, const std::map<int, std::pair<double, double>>& cities, int start, int end)
{
    std::vector<int> segment(tour.begin() + start, tour.begin() + end);
    two_opt(segment, cities);

    std::copy(segment.begin(), segment.end(), tour.begin() + start);
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
    
    int num_threads = std::thread::hardware_concurrency();
    int segment_size = tour.size() / num_threads;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) 
    {
        int start = i * segment_size;
        int end = (i == num_threads - 1) ? tour.size() : (i + 1) * segment_size;
        threads.emplace_back(optimizeSegment, std::ref(tour), std::ref(cities), start, end);
    }

    for (auto& t : threads)
        t.join();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution Time: " << elapsed.count() << " seconds\n";

    double optimizedDistance = calculateTotalDistance(tour, cities);
    std::cout << "Optimized Tour Distance: " << optimizedDistance << "\n";

    return 0;
}
