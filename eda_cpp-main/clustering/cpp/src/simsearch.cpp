#include "simsearch.hpp"
#include "utils.hpp"
#include <queue>

SimSearch::SimSearch(const Matrix &_mat_data, const Cluster &_mat_clusters): mat_data(_mat_data), mat_clusters(_mat_clusters) {
    // Constructor implementation (if needed)
}

std::vector<size_t> SimSearch::search_with_clusters(const float *query, size_t top_k) {
    size_t N = mat_data.getN();
    if(top_k == 0 || N == 0) return {};
    if(top_k > N) top_k = N;

    size_t k = mat_clusters.getK();
    if (k == 0) return search_without(query, top_k);

    size_t d = mat_data.getDim();

    std::vector<float> dist_centorids;
    dist_centorids.resize(k);
    size_t c;
    for (c = 0; c < k; c++) {
        const float *centroid = mat_clusters.getCentroid(c);
        dist_centorids[c] = vec_compute_distance(query, centroid, d);
    }

    std::vector<size_t> centroids_order = argsort(dist_centorids);

    std::vector<std::pair<float, size_t>> candidates;
    candidates.reserve(top_k * 2);

    size_t pos;
    for(pos = 0; pos < centroids_order.size(); pos++) {
        size_t cid = centroids_order[pos];
        const std::vector<size_t> &members = mat_clusters.getInds(cid);

        size_t m;
        for(m = 0; m < members.size(); m++) {
            size_t idxP = members[m];
            const float *row = mat_data.getRow(idxP);
            float d = vec_compute_distance(query, row, d);
            candidates.push_back(std::pair<float, size_t>(d, idxP));
        }

        if (candidates.size() >= top_k) {
            break;
        }
    }

    if(candidates.empty()) return {};

    std::vector<float> distances;
    distances.reserve(candidates.size());
    size_t i;
    for (i = 0; i < candidates.size(); i++) {
        distances[i] = candidates[i].first;
    }

    std::vector<size_t> order = argsort(distances);

    size_t limit = top_k;
    if(candidates.size() < limit) {
        limit = candidates.size();
    }

    std::vector<size_t> result;
    result.reserve(limit);
    for( i = 0; i < limit; i++) {
        size_t pos_in_candidates = order[i];
        result.push_back(candidates[pos_in_candidates].second);
    }
    return result;
}

std::vector<size_t> SimSearch::search_without(const float *query, size_t top_k) {
    size_t N = mat_data.getN();
    if(top_k == 0 || N == 0) return {};
    if(top_k > N) top_k = N;

    size_t d = mat_data.getDim();

    std::vector<std::pair<float, size_t>> all;
    all.reserve(N);

    size_t i;
    for(i = 0; i < N; i++) {
        const float *row = mat_data.getRow(i);
        float dist = vec_compute_distance(query, row, d);
        all.push_back(std::pair<float, size_t>(dist, i));
    }

    std::vector<float> distances;
    distances.resize(N);
    for (i = 0; i < N; i++) {
        distances[i] = all[i].first;
    }

    std::vector<size_t> order = argsort(distances);

    std::vector<size_t> result;
    result.reserve(top_k);
    for( i = 0; i < top_k; i++) {
        size_t pos_in_all = order[i];
        result.push_back(all[pos_in_all].second);
    }
    return result;
}
