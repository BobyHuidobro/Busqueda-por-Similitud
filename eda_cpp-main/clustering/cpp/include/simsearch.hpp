#ifndef _SIM_SEARCH_HPP
#define _SIM_SEARCH_HPP
#include "cluster.hpp"

class SimSearch{
private:
    const Matrix &mat_data; // data to be processed
    const Cluster &mat_clusters; // centroids to be used

public:
    SimSearch(const Matrix &_mat_data, const Cluster &_mat_clusters);
    std::vector<size_t> search_with_clusters(const float *query, size_t top_k);
    std::vector<size_t> search(const float *query, size_t top_k);
};

#endif 