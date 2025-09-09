#include "utils.hpp"
#include "matrix.hpp"
#include "cluster.hpp"
#include "simsearch.hpp"

#include <iostream>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <algorithm>

using Clock = std::chrono::steady_clock;
using Ms = std::chrono::duration<double, std::milli>;

struct Metrics {
    double comps_sum = 0.0;          
    double t_no_sort_ms_sum = 0.0;   
    double t_with_sort_ms_sum = 0.0; 
    double error_sum = 0.0;          
    size_t count = 0;                

    void add(size_t comps, double t_no_sort, double t_with_sort, double err) {
        comps_sum += static_cast<double>(comps);
        t_no_sort_ms_sum += t_no_sort;
        t_with_sort_ms_sum += t_with_sort;
        error_sum += err;
        count++;
    }

    double comps_avg() const { return (count ? comps_sum / count : 0.0); }
    double t_no_sort_ms_avg() const { return (count ? t_no_sort_ms_sum / count : 0.0); }
    double t_with_sort_ms_avg() const { return (count ? t_with_sort_ms_sum / count : 0.0); }
    double error_avg() const { return (count ? error_sum / count : 0.0); }
};

static std::vector<size_t> topk_linear_measure(
    const Matrix& data, const float* q, size_t top_k,
    size_t& out_comps, double& out_t_no_sort_ms, double& out_t_with_sort_ms)
{
    const size_t N = data.getN();
    const size_t d = data.getDim();
    if (top_k == 0 || N == 0) {
        out_comps = 0;
        out_t_no_sort_ms = 0.0;
        out_t_with_sort_ms = 0.0;
        return {};
    }
    size_t k = std::min(top_k, N);

    std::vector<float> distances(N);
    auto t0 = Clock::now();
    for (size_t i = 0; i < N; ++i) {
        const float* row = data.getRow(i);
        distances[i] = vec_compute_distance(q, row, d);
    }
    auto t1 = Clock::now();

    auto order = argsort(distances);
    auto t2 = Clock::now();

    std::vector<size_t> result;
    result.reserve(k);
    for (size_t i = 0; i < k; ++i) {
        result.push_back(order[i]);
    }
    auto t3 = Clock::now();

    out_comps = N; // comparaciones a data
    out_t_no_sort_ms = Ms(t1 - t0).count();
    out_t_with_sort_ms = Ms(t3 - t0).count();
    return result;
}


static std::vector<size_t> topk_clustered_measure(
    const Matrix& data, const Cluster& cluster, const float* q, size_t top_k,
    size_t& out_comps, double& out_t_no_sort_ms, double& out_t_with_sort_ms)
{
    const size_t N = data.getN();
    const size_t d = data.getDim();
    if (top_k == 0 || N == 0) {
        out_comps = 0;
        out_t_no_sort_ms = 0.0;
        out_t_with_sort_ms = 0.0;
        return {};
    }
    size_t k_top = std::min(top_k, N);

    const size_t K = cluster.getK();
    if (K == 0) {

        return topk_linear_measure(data, q, k_top, out_comps, out_t_no_sort_ms, out_t_with_sort_ms);
    }


    auto t0 = Clock::now();
    std::vector<float> dist_centroids(K);
    for (size_t c = 0; c < K; ++c) {
        const float* centroid = cluster.getCentroid(c);
        dist_centroids[c] = vec_compute_distance(q, centroid, d);
    }
    auto t1 = Clock::now();


    auto order_centroids = argsort(dist_centroids);
    auto t2 = Clock::now();


    std::vector<size_t> candidates;
    candidates.reserve(k_top * 2);

    for (size_t pos = 0; pos < order_centroids.size() && candidates.size() < k_top; ++pos) {
        size_t cid = order_centroids[pos];
        std::vector<size_t> members = cluster.getInds(cid);
        for (size_t m = 0; m < members.size(); ++m) {
            candidates.push_back(members[m]);
            if (candidates.size() >= k_top) break;
        }
    }


    std::vector<float> cand_dist(candidates.size());
    for (size_t i = 0; i < candidates.size(); ++i) {
        const float* row = data.getRow(candidates[i]);
        cand_dist[i] = vec_compute_distance(q, row, d);
    }
    auto t3 = Clock::now();


    auto order_cand = argsort(cand_dist);
    std::vector<size_t> result;
    result.reserve(k_top);
    for (size_t i = 0; i < k_top && i < order_cand.size(); ++i) {
        result.push_back(candidates[order_cand[i]]);
    }
    auto t4 = Clock::now();


    size_t comps = candidates.size();


    double t_dist_only_ms = Ms((t1 - t0) + (t3 - t2)).count(); 
    double t_total_ms = Ms(t4 - t0).count(); 

    out_comps = comps;
    out_t_no_sort_ms = t_dist_only_ms;
    out_t_with_sort_ms = t_total_ms;
    return result;
}

static double error_one(const std::vector<size_t>& truth, const std::vector<size_t>& approx, size_t m) {

    std::unordered_set<size_t> T(truth.begin(), truth.end());
    size_t inter = 0;
    for (size_t x : approx) {
        if (T.find(x) != T.end()) inter++;
    }
    double denom = static_cast<double>(m == 0 ? 1 : m);
    double recall = static_cast<double>(inter) / denom;
    return 1.0 - recall;
}

int main() {
    std::ios::sync_with_stdio(false);


    const std::string dataPath = "../../data_eda.npy";
    const std::string queriesPath = "../../queries_eda.npy";

    Matrix data(dataPath);
    Matrix queries(queriesPath);

    if (data.getDim() == 0 || queries.getDim() == 0) {
        std::cerr << "Error cargando .npy. Verifica rutas:\n"
                  << "  data: " << dataPath << "\n"
                  << "  queries: " << queriesPath << "\n";
        return 1;
    }
    if (data.getDim() != queries.getDim()) {
        std::cerr << "Error: dim(data) != dim(queries)\n";
        return 1;
    }

    std::cout << "data: rows=" << data.getN() << " dim=" << data.getDim() << "\n";
    std::cout << "queries: rows=" << queries.getN() << " dim=" << queries.getDim() << "\n";


    const std::vector<size_t> Ks = {0, 8, 16, 32, 64, 128};

    const std::vector<size_t> Ms_topk = {8, 16, 32, 64, 128};


    std::vector<Cluster> clusters;
    clusters.reserve(Ks.size());
    for (size_t K : Ks) {
        if (K == 0) {

            clusters.emplace_back(data, 1);
        } else {
            clusters.emplace_back(data, K);         
            clusters.back().applyClustering();      
        }
    }


    const size_t max_m = Ms_topk.back();
    SimSearch sim_linear(data, clusters[0]); 
    std::vector<std::vector<size_t>> truth_all; 
    truth_all.reserve(queries.getN());
    for (size_t qi = 0; qi < queries.getN(); ++qi) {
        const float* q = queries.getRow(qi);
        auto truth = sim_linear.search(q, std::min(max_m, data.getN()));
        truth_all.push_back(std::move(truth));
    }

    std::cout << std::fixed << std::setprecision(4);

    for (size_t m : Ms_topk) {
        std::cout << "\nTabla de resultados para m = " << m << "\n";
        std::cout << "k-clusters\t#comparaciones\t\ttiempo (ms) sin ordenar\t\ttiempo (ms) con ordenar\t\terror\n";

        for (size_t i = 0; i < Ks.size(); ++i) {
            size_t K = Ks[i];
            Metrics agg;

            for (size_t qi = 0; qi < queries.getN(); ++qi) {
                const float* q = queries.getRow(qi);
                size_t comps = 0;
                double t_no_sort = 0.0, t_with_sort = 0.0;
                std::vector<size_t> res;

                if (K == 0) {
                    res = topk_linear_measure(data, q, m, comps, t_no_sort, t_with_sort);
                    // error contra ground truth debe ser 0
                    agg.add(comps, t_no_sort, t_with_sort, 0.0);
                } else {
                    res = topk_clustered_measure(data, clusters[i], q, m, comps, t_no_sort, t_with_sort);
                    // verdad de esta query para m
                    const auto& truth_full = truth_all[qi];
                    std::vector<size_t> truth_m;
                    truth_m.reserve(std::min(m, truth_full.size()));
                    for (size_t j = 0; j < std::min(m, truth_full.size()); ++j) truth_m.push_back(truth_full[j]);

                    double err = error_one(truth_m, res, std::min(m, data.getN()));
                    agg.add(comps, t_no_sort, t_with_sort, err);
                }
            }

            std::cout << K
                      << "\t\t" << static_cast<size_t>(agg.comps_avg())
                      << "\t\t\t" << agg.t_no_sort_ms_avg()
                      << "\t\t\t\t" << agg.t_with_sort_ms_avg()
                      << "\t\t\t\t" << agg.error_avg()
                      << "\n";
        }
    }


    {
        std::vector<float> v = {5, 2, 7, 1};
        auto idx = argsort(v); 
        std::cout << "\nargsort([5,2,7,1]) => [";
        for (size_t i = 0; i < idx.size(); ++i) {
            std::cout << idx[i] << (i + 1 < idx.size() ? ", " : "");
        }
        std::cout << "]\n";
    }

    return 0;
}