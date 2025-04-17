#include "../tools/biGraph.hpp"
#ifndef _MBCS_H
#define _MBCS_H

class linearSet {
private:
    uint32_t * vSet = nullptr;
    uint32_t * pSet = nullptr;
    uint32_t sz = 0;

public:
    uint32_t x, c, r;

    linearSet() {}
    linearSet(uint32_t size):sz(size), c(0), x(0), r(0) {
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }
    void resize(uint32_t size) {
        sz = size;
        c = 0;
        x = 0;
        r = 0;
        vSet = new uint32_t[size];
        pSet = new uint32_t[size];

        for(uint32_t i = 0; i < size; i++) {
            vSet[i] = pSet[i] = i;
        }
    }

    void init() {
        for(uint32_t i = 0; i < c; i++) {
            pSet[vSet[i]] = i;
        }
    }

    void reset() {
        c = 0;
        x = 0;
        r = 0;
    }

    ~linearSet() {
        if (vSet != nullptr) delete [] vSet;
        if (pSet != nullptr) delete [] pSet;
    }

    void swapByPos(uint32_t i, uint32_t j) {
        std::swap(vSet[i], vSet[j]);
        std::swap(pSet[vSet[i]], pSet[vSet[j]]);
    }

    uint32_t operator [] (uint32_t i) {
        if(i >= sz) {
            printf("%u %u\n", i, sz);fflush(stdout);
            exit(-1);
        }
        return  vSet[i];
    }

    uint32_t pos(uint32_t v) {
        return pSet[v];
    }

    bool CIsEmpty() { return c == r; }
    bool XIsEmpty() { return x == c; }

    uint32_t cSize() {
        return c - r;
    }

    uint32_t * begin() { return vSet; }

    // void print() {
    //     printf("X:");
    //     for(uint32_t i = x; i < c; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    //     printf("C:");
    //     for(uint32_t i = c; i < r; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    //     printf("R:");
    //     for(uint32_t i = r; i < sz; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    // }

    void print() {
        printf("R:");
        for(uint32_t i = 0; i < r; i++) {
            printf("%u ", vSet[i]);
        }
        printf("\nC:");
        for(uint32_t i = r; i < c; i++) {
            printf("%u ", vSet[i]);
        }
        printf("\n");
    }

    // void print(uint32_t x, uint32_t c, uint32_t r) {
    //     printf("X:");
    //     for(uint32_t i = x; i < c; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    //     printf("C:");
    //     for(uint32_t i = c; i < r; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    //     printf("R:");
    //     for(uint32_t i = r; i < sz; i++) {
    //         printf("%u ", vSet[i]);
    //     }
    // }

    void print(uint32_t r, uint32_t c) {
        printf("R:");
        for(uint32_t i = 0; i < r; i++) {
            printf("%u ", vSet[i]);
        }
        printf("\nC:");
        for(uint32_t i = r; i < c; i++) {
            printf("%u ", vSet[i]);
        }
        printf("\n");
    }

    uint32_t capacity(){ return sz;}
};

class MBCS
{
private:
    biGraph *g = NULL, *gs = NULL;
    uint32_t mSize = 0, realsz[2] = {0}, sz[2] = {1}, csz[2] = {1}, gt = 0, mxdeep=0;
    unsigned long numbranches = 0, maxsubbranches = 0, subbcnt = 0;
    double stp = 1.6;
    bool BIBRANCH = false;
    int alg;
    uint32_t cvlb = 0;


    std::vector<uint32_t> sdeg[2], bin[2], _deg[2];
    std::vector<uint32_t> realRes[2];
    std::vector<uint32_t> vqueue[2];

    linearSet S[2], Sc[2];
    std::vector<std::vector<uint32_t>> ws;

    uint32_t lrs[2] = {0};

    void printRes();
    uint32_t upperCan(uint32_t deep);
    void gemSize(uint32_t t, uint32_t tr, uint32_t zr, uint32_t ne);

    void branch(uint32_t deep, uint32_t ne);
    void branchNew(uint32_t deep, uint32_t ne);
    void branchSubG(uint32_t deep, uint32_t ne);
    void branch_minCover(uint32_t u, uint32_t deep, uint32_t ne);
    void initSize();

    void reSetdeg(uint32_t *c);
    void reSetdeg(uint32_t t, uint32_t *c);

    void graph_convert();
    void subGraph_generate();
    bool is_nbr_sub(uint32_t t, uint32_t u, uint32_t v);
    bool is_nbr_cvert(uint32_t t, uint32_t u, uint32_t v);
    void sort_deg(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t);
public:
    MBCS(/* args */);
    MBCS(const std::string & filePath, uint32_t ls, uint32_t rs);
    ~MBCS();

    void run(int alg);
    void sort_in(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t);
    void sort_de(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t);

    void setStp(double s){ stp = s; assert(stp > 1);}

    uint32_t find_max (int n1, int n2, int x) {
        if (n1 > n2) {
            if (n1-x >= n2) return (n1-x)*n2;
            else {
                x = x+n2-n1;
                return (n2-x/2)*(n2-x/2+x%2);
            }
        }
        else {
            if (n2-x >= n1) return (n2-x)*n1;
            else {
                x = x+n1-n2;
                return (n1-x/2)*(n1-x/2+x%2);
            }
        }
    }
};

#endif