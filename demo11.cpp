#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>
#include <algorithm>
#include <limits>
#include <ctime>

struct Point {
    double latitude;
    double longitude;
};

double distance(const Point& p1, const Point& p2) {
    double latDiff = p1.latitude - p2.latitude;
    double lonDiff = p1.longitude - p2.longitude;
    return std::sqrt(latDiff * latDiff + lonDiff * lonDiff);
}

int findClosestCentroid(const std::vector<Point>& centroids, const Point& point) {
    int closestIndex = 0;
    double minDistance = distance(centroids[0], point);

    for (int i = 1; i < centroids.size(); ++i) {
        double dist = distance(centroids[i], point);
        if (dist < minDistance) {
            minDistance = dist;
            closestIndex = i;
        }
    }

    return closestIndex;
}

std::vector<Point> calculateCentroids(const std::vector<std::vector<Point>>& clusters) {
    std::vector<Point> centroids;

    for (const auto& cluster : clusters) {
        double sumLat = 0.0;
        double sumLon = 0.0;

        for (const auto& point : cluster) {
            sumLat += point.latitude;
            sumLon += point.longitude;
        }

        double avgLat = sumLat / cluster.size();
        double avgLon = sumLon / cluster.size();

        centroids.push_back({ avgLat, avgLon });
    }

    return centroids;
}

int findLargestCluster(const std::vector<std::vector<Point>>& clusters) {
    int largestClusterIndex = 0;
    int maxClusterSize = clusters[0].size();

    for (int i = 1; i < clusters.size(); ++i) {
        if (clusters[i].size() > maxClusterSize) {
            maxClusterSize = clusters[i].size();
            largestClusterIndex = i;
        }
    }

    return largestClusterIndex;
}

void writeClusterDataToCSV(const std::vector<std::vector<Point>>& clusters, const std::string& filename) {
    std::ofstream outputFile(filename);

    if (!outputFile) {
        std::cout << "Failed to create the output file." << std::endl;
        return;
    }

    for (int i = 0; i < clusters.size(); ++i) {
        for (const auto& point : clusters[i]) {
            outputFile << point.latitude << "," << point.longitude << "," << i << std::endl;
        }
    }

    outputFile.close();
}

int main() {
    std::ifstream inputFile("CSE299_G4_TEXTFILE.txt");
    if (!inputFile) {
        std::cout << "Failed to open file." << std::endl;
        return 1;
    }

    std::vector<Point> dataPoints;
    std::string line;

    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string latStr, lonStr;

        if (std::getline(ss, lonStr, ',') && std::getline(ss, latStr, ',')) {
            try {
                double lat = std::stod(latStr);
                double lon = std::stod(lonStr);

                if (lat >= -90.0 && lat <= 90.0 && lon >= -180.0 && lon <= 180.0) {
                    dataPoints.push_back({ lat, lon });
                }
            } catch (...) {
                continue;
            }
        }
    }

    inputFile.close();

    std::random_device rd;
    std::mt19937 gen(rd());
    gen.seed(std::time(nullptr));
    std::uniform_int_distribution<> dist(30, 200);
    int numCustomers = dist(gen);
    std::vector<Point> selectedDataPoints;

    std::shuffle(dataPoints.begin(), dataPoints.end(), gen);

    for (int i = 0; i < numCustomers; ++i) {
        selectedDataPoints.push_back(dataPoints[i]);
    }

    const int minNumClusters = 2;  // Minimum number of clusters to try
    const int maxNumClusters = 10; // Maximum number of clusters to try
    const int numIterations = 500; // Number of iterations for k-means

    std::uniform_int_distribution<> dis(0, dataPoints.size() - 1);
    std::vector<double> distortions;

    for (int k = minNumClusters; k <= maxNumClusters; ++k) {
        std::vector<Point> centroids;
        std::vector<std::vector<Point>> clusters(k);

        for (int i = 0; i < k; ++i) {
            centroids.push_back(selectedDataPoints[dis(gen)]);
        }

        for (int iter = 0; iter < numIterations; ++iter) {
            for (int i = 0; i < selectedDataPoints.size(); ++i) {
                int closestCentroid = findClosestCentroid(centroids, selectedDataPoints[i]);
                clusters[closestCentroid].push_back(selectedDataPoints[i]);
            }

            centroids = calculateCentroids(clusters);

            for (auto& cluster : clusters) {
                cluster.clear();
            }
        }

        double distortion = 0.0;
        for (int i = 0; i < k; ++i) {
            for (const auto& point : clusters[i]) {
                distortion += distance(point, centroids[i]);
            }
        }

        distortions.push_back(distortion);
    }

    int bestNumClusters = minNumClusters;
    double minDistortion = distortions[minNumClusters - 1];

    for (int i = minNumClusters; i < distortions.size(); ++i) {
        if (distortions[i] < minDistortion) {
            minDistortion = distortions[i];
            bestNumClusters = i + 1;
        }
    }

    std::cout << "Best number of clusters: " << bestNumClusters << std::endl;

    std::vector<Point> centroids;
    std::vector<std::vector<Point>> clusters(bestNumClusters);

    centroids.push_back(selectedDataPoints[dis(gen)]);

    for (int i = 1; i < bestNumClusters; ++i) {
        double maxProbability = 0.0;
        int nextCentroidIndex = 0;

        for (int j = 0; j < selectedDataPoints.size(); ++j) {
            double minDistance = distance(centroids[0], selectedDataPoints[j]);

            for (int k = 1; k < centroids.size(); ++k) {
                double dist = distance(centroids[k], selectedDataPoints[j]);
                if (dist < minDistance) {
                    minDistance = dist;
                }
            }

            double probability = minDistance * minDistance;
            if (probability > maxProbability) {
                maxProbability = probability;
                nextCentroidIndex = j;
            }
        }

        centroids.push_back(selectedDataPoints[nextCentroidIndex]);
    }

    for (const auto& point : selectedDataPoints) {
        int closestCentroid = findClosestCentroid(centroids, point);
        clusters[closestCentroid].push_back(point);
    }

    int largestClusterIndex = findLargestCluster(clusters);

    std::cout << "Number of data points: " << selectedDataPoints.size() << std::endl;
    std::cout << "Largest cluster size: " << clusters[largestClusterIndex].size() << std::endl;

    std::string filename = "clusters.csv";
    writeClusterDataToCSV(clusters, filename);
    std::cout << "Cluster data written to: " << filename << std::endl;

    return 0;
}
