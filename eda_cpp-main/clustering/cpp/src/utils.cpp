#include "utils.hpp"
#include <cmath>
#include <iostream>
#include <vector>

float vec_compute_distance(const float* u, const float* v, size_t dim){
    float d = 0;
    float s = 0;
    float diff = 0;
    for (size_t i = 0;  i < dim; i++){
        diff = (u[i] - v[i]);
        s = s + diff * diff;
    }
    d = std::sqrt(s);
    return d;
}

void vec_add(float* s, const float* u,  size_t dim){
    for (size_t i = 0;  i < dim; i++){
        s[i] = u[i] + s[i];
    }
}

void vec_divide(float* u, float scalar, size_t dim){
    for (size_t i = 0;  i < dim; i++){
        u[i] = u[i] / scalar;
    }
}

void set_array(float *array, size_t dim, float val){
    for (size_t i = 0;  i < dim; i++){
        array[i] = val;
    }
}

float vec_compute_avg_dif(const float *u, const float* v,  size_t dim){

    float dif = 0;
    for (size_t i = 0;  i < dim; i++){        
        dif = dif + std::abs(u[i] - v[i]);
        //std::cout << u[i] << " - " << v[i] << "  " << std::abs(u[i] - v[i]) << std::endl;
    }   
    //std::cout << dif <<std::endl;
    return dif / dim;
}

void print_array(const float *array, size_t d){
    for (size_t i = 0; i< d; i ++){
        std::cout << array[i] << " ";
    }
    std::cout<< std::endl;
}

static int partition(const std::vector<float> &values, std::vector<size_t> &indices, int low, int high) {
    float pivotVal= values[indices[high]];
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if(values[indices[j]] <= pivotVal) {
            i++;
            std::swap(indices[i], indices[j]);
        }
    }
    std::swap(indices[i + 1], indices[high]);
    return i + 1;
}

static void quicksort_indices(const std::vector<float> &values, std::vector<size_t> &indices, int low, int high) {
    if(low < high) {
        int p = partition(values, indices, low, high);
        quicksort_indices(values, indices, low, p - 1);
        quicksort_indices(values, indices, p + 1, high);
    }
}

std::vector<size_t> argsort(const std::vector<float> &values) {
    std::vector<size_t> indices(values.size());
    for (size_t i = 0; i < values.size(); i++) {
        indices[i] = i;
    }
    quicksort_indices(values, indices, 0, static_cast<int>(values.size()) - 1);
    return indices;
}