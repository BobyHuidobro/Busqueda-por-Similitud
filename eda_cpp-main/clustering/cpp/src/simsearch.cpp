#include "simsearch.hpp"
#include "utils.hpp"

SimSearch::SimSearch(const Matrix &_mat_data, const Cluster &_mat_clusters): mat_data(_mat_data), mat_clusters(_mat_clusters) {
    // Constructor implementation (if needed)
}

std::vector<size_t> SimSearch::search_with_clusters(const float *query, size_t top_k) {
    size_t N = mat_data.getN();
    if (top_k == 0 || N == 0) return {};
    if (top_k > N) top_k = N;

    size_t K = mat_clusters.getK();
    if(K == 0) {
        return search_without(query, top_k);
    }

    size_t D = mat_data.getDim();

    std::vector<float> dist_centroids(K);
    std::vector<size_t> order(K);
    for (size_t c = 0; c < K; c++) {
        order[c] = c;
        const float *centroid = mat_clusters.getCentroid(c);
        dist_centroids[c] = vec_compute_distance(query, centroid, D);
    }

    std::sort(order.begin(), order.end(), [&](size_t a, size_t b){
        return dist_centroids[a] < dist_centroids[b];
    });

    using Node = std::pair<float, size_t>;
    auto cmp = [](const Node &A, const Node &B)

}

std::vector<size_t> SimSearch::search_without(const float *query, size_t top_k) {
    // Implementation of search without using clusters
    std::vector<size_t> results;
    // ... (search logic here)

    return results;
}