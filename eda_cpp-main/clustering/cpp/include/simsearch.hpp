#ifndef _SIM_SEARCH_HPP
#define _SIM_SEARCH_HPP

class SimSearch{
private:
    const Matrix &mat_data; // data to be processed
    const Matrix &mat_clusters; // centroids to be used

public:
    SimSearch(Matrix, clusters);
    std::vector<size_t> search_with_clusters(const float *query,  top_k);
    std::vector<size_t> search_without(const float *query, top_k);
}

#endif 