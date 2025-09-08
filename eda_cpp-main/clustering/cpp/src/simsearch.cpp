#include "simsearch.hpp"

SimSearch::SimSearch(const Matrix &_mat_data, const Cluster &_mat_clusters): mat_data(_mat_data), mat_clusters(_mat_clusters) {
    // Constructor implementation (if needed)
}

std::vector<size_t> SimSearch::search_with_clusters(const float *query, size_t top_k) {
    // Implementation of search using clusters
    std::vector<size_t> results;
    // ... (search logic here)

    return results;
}

std::vector<size_t> SimSearch::search_without(const float *query, size_t top_k) {
    // Implementation of search without using clusters
    std::vector<size_t> results;
    // ... (search logic here)

    return results;
}