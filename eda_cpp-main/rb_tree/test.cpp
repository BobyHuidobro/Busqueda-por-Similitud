#include "trees/rb.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

int* readKeysFromFile(std::string filename, int* n_keys){
    std::ifstream  fin(filename, std::ios::binary); 
    char* val = new char[4];
    int n = 0;
    fin.read(val, 4);
    while (!fin.eof()){
        n = n + 1;        
        fin.read(val, 4);
    }
    fin.close();
    fin.open(filename, std::ios::binary); 
    int* keys = new int[n];
    for(int i=0; i < n; i++){
        fin.read(val, 4);
        keys[i] = *reinterpret_cast<int*>(val);
    }
    fin.close();
    *n_keys = n;
    delete[] val;
    return keys;
}

int main(int nargs, char** vargs){
    int n_data = 0;
    std::vector<std::string> vec(5);
    vec[0] = "/home/jmsaavedrar/Research/git/otros/keys_1024.bin";
    vec[1] = "/home/jmsaavedrar/Research/git/otros/keys_32768.bin";
    vec[2] = "/home/jmsaavedrar/Research/git/otros/keys_1048576.bin"; 
    vec[3] = "/home/jmsaavedrar/Research/git/otros/keys_33554432.bin";
    vec[4] = "/home/jmsaavedrar/Research/git/otros/queries_1000.bin";
    for (auto f:vec){
        int* data = readKeysFromFile(f, &n_data);
        std::cout<<n_data <<std::endl;
        delete[] data;
    }

    // trees::RB rbtree;
    // for(int i=0; i<n_data;i++){
    //     std::cout<< "inserting " << data[i] << std::endl;
    //     rbtree.insert(data[i]);
    // } 
    // rbtree.traverse();

    return 0;
}