#pragma once
#ifndef _BIGRAPH_HPP
#define _BIGRAPH_HPP

#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>
#include <chrono>
#include<cstring>
#include <cassert>
#include <emmintrin.h>

#define unfilled -1
class CuckooHash
{
private:
    /* data */
    uint32_t capacity;
    uint32_t mask;
    uint32_t size;
    const uint32_t buff_size = sizeof(uint32_t);
    uint32_t *hashtable;

    void rehash(uint32_t **_table) {
        uint32_t oldcapacity = capacity;
        mask = mask == 0 ? 1 : ((mask << 1) | 1);
        capacity = (mask + 1) * buff_size;
        uint32_t *newhash = new uint32_t[capacity];
        memset((newhash), unfilled, sizeof(uint32_t) * capacity);
        for (uint32_t i = 0; i < oldcapacity; ++i){
            if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
        }
        std::swap((*_table), newhash);
        delete[] newhash;
    }
    
    void insert(const uint32_t &_u, uint32_t **_table) {
        
        uint32_t hs = hash1(_u);
        for (uint32_t i = 0; i < buff_size; ++i) {
            if ((*_table)[hs * buff_size + i] == unfilled){
                (*_table)[hs * buff_size + i] = _u;
                return;
            }
        }
        hs = hash2(_u);
        for (uint32_t i = 0; i < buff_size; ++i) {
            if ((*_table)[hs * buff_size + i] == unfilled){
                (*_table)[hs * buff_size + i] = _u;
                return;
            }
        }

        bool use_hash1 = true;
        uint32_t u = _u;
        for (uint32_t i = 0; i < mask; ++i) {
            uint32_t replaced;
            if (use_hash1) hs = hash1(u);
            else hs = hash2(u);
            uint32_t j = 0;
            for (; j < buff_size; ++j) {
                if ((*_table)[hs * buff_size + j] == unfilled) break;
            }
            if (buff_size == j) {
                replaced = std::move((*_table)[hs * buff_size]);
                j = 1;
                for (; j < buff_size; j++) {
                    (*_table)[hs * buff_size + j - 1] =
                        std::move((*_table)[hs * buff_size + j]);
                }
                (*_table)[hs * buff_size + j - 1] = u;
            }
            else {
                replaced = std::move((*_table)[hs * buff_size + j]);
                (*_table)[hs * buff_size + j] = u;
            }
            use_hash1 = hs == hash2(replaced);
            u = std::move(replaced);
            if (u == unfilled) return;
        }
        rehash(_table);
        insert(u, _table);
    }

    uint32_t hash1(const uint32_t &x) { return x & mask;}
    uint32_t hash2(const uint32_t &x) { return ~x & mask;}

public:
    CuckooHash(/* args */) {
        capacity = 0;
        hashtable = NULL;
        mask = 0;
        size = 0;
    }
    ~CuckooHash() {
        if (hashtable) delete[] hashtable;
    }

    void reserve(uint32_t _size) {
        if (capacity >= _size) return;
        mask = mask == 0 ? 1 : ((mask << 1) | 1);
        while (_size >= mask * buff_size) mask = (mask << 1) | 1;
        capacity = (mask + 1) * buff_size;
        if (hashtable) delete[] hashtable;
        hashtable = new uint32_t[capacity];
        memset(hashtable, unfilled, sizeof(uint32_t) * capacity);
    }

    void insert(const uint32_t &_u) {
        if (find(_u)) return;
        insert(_u, &hashtable);
        size++;
    }

    bool find(const uint32_t &_u) {
        uint32_t hs1 = hash1(_u);
        uint32_t hs2 = hash2(_u);
// if(!(buff_size == 4 && sizeof(uint32_t) == 4))
// printf("%u %u\n", buff_size, sizeof(uint32_t));fflush(stdout);
        // assert(buff_size == 4 && sizeof(uint32_t) == 4);
        __m128i cmp = _mm_set1_epi32(_u);
        __m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
        __m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
        __m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

        return _mm_movemask_epi8(flag) != 0;
    }
    uint32_t getcapacity() {return capacity;}
    uint32_t getsize() {return size;}
    uint32_t getmask() {return mask;}
    uint32_t *gethashtable() {return hashtable;}
};

struct biGraph
{
    uint32_t n1, n2, m, maxDu = 0, maxDv = 0;
    uint32_t n[2] = {0}, sz[2] = {1}, mxd[2] = {0};

    struct Edge{
        uint32_t u, v;
    };
    std::vector<Edge> edges;
    
    std::vector<uint32_t> pU, e1, pV, e2;
    std::vector<uint32_t> p[2];
    std::vector<uint32_t> e[2];
    std::vector<uint32_t> cores[2];

    std::vector<CuckooHash> nbrIndex[2];

    std::vector<uint32_t> realLable[2];
    std::vector<uint32_t> colors[2], colorOrdering[2];

    uint32_t lrs[2];
    uint32_t k = 0;
    uint32_t maxcore = 0, mxcl = 0;
    uint32_t maxcols[2] = {0};

    biGraph(/* args */) {

    }

    biGraph(const std::string & filePath, uint32_t ls, uint32_t rs) {
        sz[0] = ls; sz[1] = rs;
        read(filePath);
        reduction_and_relabel();
        // degreeOrdering();
        nbrIndex[0].resize(n[0]);
        nbrIndex[1].resize(n[1]);
        printf("n1=%d,n2=%d\n", n1,n2);
        for(uint32_t t = 0; t <= 1; t++) {
            for(uint32_t u = 0; u < n[t]; u++) {
                uint32_t d = p[t][u + 1] - p[t][u];
                nbrIndex[t][u].reserve(d+1);
                for (uint32_t j = p[t][u]; j < p[t][u + 1];++j)
                    nbrIndex[t][u].insert(e[t][j]);
            }
        }
    }

    ~biGraph() {

    }

    void read(const std::string & filePath) {
        // printf("file: %s\n", filePath.c_str());
    
        bool is_bin = false;
        auto t1 = std::chrono::steady_clock::now();
        if (strstr(filePath.c_str(),".bin")) is_bin = true;
        if (is_bin) {
            FILE *in = fopen(filePath.c_str(), "rb");
            if (in == NULL) {
                printf("No such file: %s\n", filePath.c_str());
                exit(1);
            }
            if(fread(&n1, sizeof(uint32_t), 1, in)!=1) printf("err: read n1!\n");
            if(fread(&n2, sizeof(uint32_t), 1, in)!=1) printf("err: read n2!\n");
            if(fread(&m, sizeof(uint32_t), 1, in)!=1) printf("err: read m!\n");

            std::vector<uint32_t> deg[2];
            deg[0].resize(n1); 
            deg[1].resize(n2);
            edges.resize(m+5);
            e1.resize(m+5);
            e2.resize(m+5);
            pU.resize(n1 + 5);
            pV.resize(n2 + 5);

            n[0] = n1; n[1] = n2;
            for (uint32_t t = 0; t <= 1; t++)
            for (uint32_t i = 0; i < n[t]; ++i)
                if (fread(&deg[t][i], sizeof(uint32_t), 1, in)!=1)
                    printf("err: read deg[%d][%d]!\n",t, i);

            for (uint32_t t = 0; t <= 1; t++){
                uint32_t cntm = 0;
                for (uint32_t i = 0; i < n[t]; ++i)
                    cntm += deg[t][m];
                assert(cntm == m);
            }

            int v, d, offs = 0;
            for (uint32_t i = 0; i < n[0]; ++i) {
                d = deg[0][i];
                for (uint32_t j = 0; j < d; j++) {
                    if (fread(&v, sizeof(uint32_t), 1, in)!=1)
                        printf("err: read adj[%d][%d]!\n",i, j);
                    edges[offs].u = i;
                    edges[offs].v = v;
                    offs++;
                }
            }
            fclose(in);

            pU[0] = 0;
            for(uint32_t u = 0; u < n1; u++) {
                pU[u + 1] = deg[0][u] + pU[u];
                maxDu = std::max(deg[0][u], maxDu);
            }
            for(uint32_t i = 0; i < m; i++) {
                e1[pU[edges[i].u]++] = edges[i].v; 
            }
            pU[0] = 0;
            for(uint32_t u = 0; u < n1; u++) {
                pU[u + 1] = deg[0][u] + pU[u];
            }

            pV[0] = 0;
            for(uint32_t v = 0; v < n2; v++) {
                pV[v + 1] = deg[1][v] + pV[v];
                maxDv = std::max(deg[0][v], maxDv);
            }
            for(uint32_t i = 0; i < m; i++) {
                e2[pV[edges[i].v]++] = edges[i].u; 
            }
            pV[0] = 0;
            for(uint32_t v = 0; v < n2; v++) {
                pV[v + 1] = deg[1][v] + pV[v];
            } 
        }
        else {
            FILE *in = fopen(filePath.c_str(), "r");
            if (in == NULL) {
                printf("No such file: %s\n", filePath.c_str());
                exit(1);
            }
            char line[128];
            fgets(line, 128, in);
            if (sscanf(line, "%d %d %d", &n1, &n2, &m) != 3) exit(1);
            assert(n1 > 0); assert(n2 > 0); assert(m > 0);

            std::vector<uint32_t> deg[2];
            deg[0].resize(n1); 
            deg[1].resize(n2);
            edges.resize(m+5);
            e1.resize(m+5);
            e2.resize(m+5);
            pU.resize(n1 + 5);
            pV.resize(n2 + 5);

            uint32_t u, v, cnt = 0;
            uint32_t maxid = 0;
            for (long i = 0; i < m; ++i) {
                char *r = fgets(line, 128, in);
                if (feof(in)) break;
                sscanf(line, "%d %d", &u, &v);
                // if (u >= v) continue;
                assert(u < n1 && u >= 0);
                assert(v < n2 && v >= 0);
                edges[cnt].u = u;
                edges[cnt].v = v; cnt++;

                deg[0][u]++; deg[1][v]++;
            }
            fclose(in);

            // printf("cnt=%d, m=%d\n", cnt, m);

            pU[0] = 0;
            for(uint32_t u = 0; u < n1; u++) {
                pU[u + 1] = deg[0][u] + pU[u];
                maxDu = std::max(deg[0][u], maxDu);
            }
            for(uint32_t i = 0; i < m; i++) {
                e1[pU[edges[i].u]++] = edges[i].v; 
            }
            pU[0] = 0;
            for(uint32_t u = 0; u < n1; u++) {
                pU[u + 1] = deg[0][u] + pU[u];
            }

            pV[0] = 0;
            for(uint32_t v = 0; v < n2; v++) {
                pV[v + 1] = deg[1][v] + pV[v];
                maxDv = std::max(deg[1][v], maxDv);
            }
            for(uint32_t i = 0; i < m; i++) {
                e2[pV[edges[i].v]++] = edges[i].u; 
            }
            pV[0] = 0;
            for(uint32_t v = 0; v < n2; v++) {
                pV[v + 1] = deg[1][v] + pV[v];
            } 
        }

        printf("Graph info: n1=%d,n2=%d,m=%d,maxd1=%d,maxd2=%d\n", n1,n2,m,maxDu, maxDv);

        n[0] = n1;
        n[1] = n2;
        p[0] = std::move(pU);
        p[1] = std::move(pV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
        this->mxd[0] = maxDu;
        this->mxd[1] = maxDv;
    }

    void reduction_and_relabel() {
        uint32_t qsz[2] = {0};
        std::vector<uint32_t> deg[2], vqueue[2];
        deg[0].resize(n[0]); deg[1].resize(n[1]);
        vqueue[0].resize(n[0]); vqueue[1].resize(n[1]);

        for (uint32_t t = 0; t <= 1; ++t)  {
            uint32_t z = t ^ 1;
            uint32_t d;
            for (uint32_t i = 0; i < n[t]; ++i) {
                d = p[t][i+1] - p[t][i];
                if (d < sz[z]) vqueue[t][qsz[t]++] = i;
                deg[t][i] = d;
            }
        }

        uint32_t r[2] = {0};
        while (r[0] < qsz[0] || r[1] < qsz[1]) {

            for (uint32_t t = 0; t <= 1; ++t) {
                uint32_t z = t^1;
                for (uint32_t i = r[t]; i < qsz[t]; ++i) {
                    uint32_t u = vqueue[t][i];
                    for (uint32_t j = p[t][u]; j < p[t][u+1]; ++j) {
                        uint32_t v = e[t][j];
                        if (deg[z][v] >= sz[t]) {
                            deg[z][v]--;
                            if (deg[z][v] < sz[t]) vqueue[z][qsz[z]++] = v;
                        }
                    }
                }
                r[t] = qsz[t];
            }
        }
        
        r[0] = 0; r[1] = 0;
        qsz[0] = n[0]; qsz[1] = n[1];

        realLable[0].resize(n[0]);
        realLable[1].resize(n[1]);
        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t z = t^1;
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t v = deg[t][u] >= sz[z] ? r[t]++ : --qsz[t]; 

                realLable[t][v] = u;
                vqueue[t][u] = v; 
            }
        }
    

        r[0] = 0; r[1] = 0;
        for (long i = 0; i < m; ++i) {
            edges[i].u = vqueue[0][edges[i].u];
            edges[i].v = vqueue[1][edges[i].v];

            if (edges[i].u < qsz[0] && edges[i].v < qsz[1]) r[0]++;
            if (edges[i].u >= qsz[0] || edges[i].v >= qsz[1]) { 
                m--;
                std::swap(edges[i], edges[m]);
                i--;
            }
        }

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t z = t^1;
            p[t][0] = 0; mxd[t] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t oldu = realLable[t][u];
                p[t][u + 1] = deg[t][oldu] + p[t][u];
                mxd[t] = std::max(deg[t][oldu], mxd[t]);
            }
        }
        
        // long mnew = 0;
        for(long i = 0; i < m; i++) {
            // if (edges[i].u < qsz[0] && edges[i].v < qsz[1]) {
                e[0][p[0][edges[i].u]++] = edges[i].v; 
                e[1][p[1][edges[i].v]++] = edges[i].u;
                // edges[mnew].u = edges[i].u;
                // edges[mnew].v = edges[i].v;
                // mnew++;
            // }
        }
        // m = mnew;

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t z = t^1;
            p[t][0] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t oldu = realLable[t][u];
                p[t][u + 1] = deg[t][oldu] + p[t][u];
                // std::sort(e[t].begin()+p[t][u], e[t].begin()+p[t][u+1]);
            }
        }
        maxDu = mxd[0]; maxDv = mxd[1];

        // for (uint32_t t = 0; t <= 1; ++t) {
        //     uint32_t z = t^1;
        //     p[t][0] = 0;
        //     for (uint32_t u = 0; u < n[t]; ++u) 
        //         deg[t][u] = 0;
        // }

        n1 = qsz[0]; n2 = qsz[1];
        this->n[0] = n1; this->n[1] = n2;
        this->mxd[0] = maxDu;
        this->mxd[1] = maxDv;

        printf("New n1=%d, n2=%d, m=%d, maxDu=%d, maxDv=%d\n", n1, n2, m, maxDu, maxDv);

        relabelCore();
        // relableCore_old();
    }

    void relabelCore(){

        std::vector<uint32_t> deg[2], bin[2], vpos[2], _lable[2], vqueue[2];

        deg[0].resize(n[0]); deg[1].resize(n[1]);
        bin[0].resize(mxd[0]+2, 0); bin[1].resize(mxd[1]+2, 0);
        cores[0].resize(n[0]);  cores[1].resize(n[1]);
        vpos[0].resize(n[0]); vpos[1].resize(n[1]);
        _lable[0].resize(n[0]); _lable[1].resize(n[1]);
        vqueue[0].resize(n[0]); vqueue[1].resize(n[1]);
        colorOrdering[0].resize(n[0]); colorOrdering[1].resize(n[1]);

        if (realLable[0].capacity() > 0) {
            memcpy(_lable[0].data(), realLable[0].data(), sizeof(uint32_t)* n[0]);
            memcpy(_lable[1].data(), realLable[1].data(), sizeof(uint32_t)* n[1]);
        }
        else {
            for (uint32_t t=0; t<=1; ++t) 
            for (uint32_t u = 0; u < n[t]; ++u)
                _lable[t][u] = u;

            realLable[0].resize(n[0]);
            realLable[1].resize(n[1]);
        }

        // printf("n1=%d, n2=%d, m=%d, maxDu=%d, maxDv=%d\n", n[0], n[1], m, maxDu, maxDv);
        
        for (uint32_t t = 0; t <= 1; ++t) {
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t d = p[t][u+1] - p[t][u];
                deg[t][u] = d;
                bin[t][d] ++;
            }
        }

        for (uint32_t t = 0; t <= 1; ++t){
            uint32_t ucnt = 0;
            for (uint32_t i = 0; i <= mxd[t]; ++i) 
                ucnt += bin[t][i];
        }

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t cnt = 0, dct = 0;
            for (uint32_t i = 0; i <= mxd[t]; ++i) {
                dct = bin[t][i];
                bin[t][i] = cnt;
                cnt += dct;
            }
            bin[t][mxd[t]+1] = n[t];
        }

        for (uint32_t t = 0; t <= 1; ++t)
        for (uint32_t u = 0; u < n[t]; ++u) {
            uint32_t d = deg[t][u];
            uint32_t i = bin[t][d]++;
            vqueue[t][i] = u;
            vpos[t][u] = i; 
        }


        uint32_t kk = 0;
        uint32_t ss[2] = {0};
        uint32_t it = 0, cnts = 0;
        while (ss[0] < n[0] || ss[1] < n[1]) {
            uint32_t t = 0, z = 1;
            while (ss[t] < bin[t][kk] || ss[z] < bin[z][kk]){
                for (t = 0; t <= 1; t++) {
                    z = t^1;
                    for (uint32_t i = ss[t]; i < bin[t][kk]; ++i) {
                        uint32_t u = vqueue[t][i];
                        cores[t][u] = kk;
                        colorOrdering[t][i] = cnts++;
                        assert(u < n[t]);
                        for (uint32_t j = p[t][u]; j < p[t][u+1]; ++j) {
                            uint32_t v = e[t][j];
                            if (deg[z][v] > kk) {
                                uint32_t dz = --deg[z][v];
                                uint32_t pz = bin[z][dz]++;
                                uint32_t vp = vpos[z][v];
                                vqueue[z][vp] = vqueue[z][pz];
                                vqueue[z][pz] = v;
                                vpos[z][v] = pz;
                                vpos[z][vqueue[z][vp]] = vp;
                            }
                        }
                    }
                    ss[t] = bin[t][kk];
                    // if (it < 10)
                    // printf("it=%d, kk=%d, ss[0]=%d, bin[0]=%d, ss[1]=%d, bin[1]=%d\n", it, kk, ss[0], bin[0][kk], ss[1], bin[1][kk]);
                }
                maxcore = std::max(kk,maxcore);
            }   
            kk++;
            it ++;
        }
        // maxcore = kk - 1;
        printf("ss[0]=%d, bin[0]=%d, ss[1]=%d, bin[1]=%d\n", ss[0], bin[0][kk], ss[1], bin[1][kk]);
        assert(ss[0] == bin[0][kk]);
        assert(ss[1] == bin[1][kk]);
        printf("maxcore=%d\n", maxcore);

        for (uint32_t t = 0; t <= 1; ++t) {
            bin[t].resize(n[t], 0);
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t up = vpos[t][u];
                realLable[t][up] = _lable[t][u];
                cores[t][up] = deg[t][u];
                bin[t][up] = p[t][u+1]-p[t][u];
            }
            p[t][0] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                p[t][u+1] = p[t][u] + bin[t][u];
            }
        }
     
        for (uint32_t i = 0; i < m; ++i) {
            uint32_t u = edges[i].u;
            uint32_t v = edges[i].v;
            edges[i].u = vpos[0][u];
            edges[i].v = vpos[1][v];

            e[0][p[0][edges[i].u]++] = edges[i].v;
            e[1][p[1][edges[i].v]++] = edges[i].u;
        }

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t z = t^1;
            p[t][0] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                p[t][u + 1] = bin[t][u] + p[t][u];
                std::sort(e[t].begin()+p[t][u], e[t].begin()+p[t][u+1]);
            }
        }
        // exit(1);
    }

    void degreeOrdering() {
        std::vector<uint32_t> deg[2], bin[2], vpos[2], _lable[2], vqueue[2];

        deg[0].resize(n[0]); deg[1].resize(n[1]);
        bin[0].resize(mxd[0]+2, 0); bin[1].resize(mxd[1]+2, 0);
        cores[0].resize(n[0]);  cores[1].resize(n[1]);
        vpos[0].resize(n[0]); vpos[1].resize(n[1]);
        _lable[0].resize(n[0]); _lable[1].resize(n[1]);
        vqueue[0].resize(n[0]); vqueue[1].resize(n[1]);
        colorOrdering[0].resize(n[0]); colorOrdering[1].resize(n[1]);

        if (realLable[0].capacity() > 0) {
            memcpy(_lable[0].data(), realLable[0].data(), sizeof(uint32_t)* n[0]);
            memcpy(_lable[1].data(), realLable[1].data(), sizeof(uint32_t)* n[1]);
        }
        else {
            for (uint32_t t=0; t<=1; ++t) 
            for (uint32_t u = 0; u < n[t]; ++u)
                _lable[t][u] = u;

            realLable[0].resize(n[0]);
            realLable[1].resize(n[1]);
        }

        // printf("n1=%d, n2=%d, m=%d, maxDu=%d, maxDv=%d\n", n[0], n[1], m, maxDu, maxDv);
        
        for (uint32_t t = 0; t <= 1; ++t) {
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t d = p[t][u+1] - p[t][u];
                deg[t][u] = d;
                bin[t][d] ++;
            }
        }

        for (uint32_t t = 0; t <= 1; ++t){
            uint32_t ucnt = 0;
            for (uint32_t i = 0; i <= mxd[t]; ++i) 
                ucnt += bin[t][i];
        }

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t cnt = 0, dct = 0;
            for (uint32_t i = 0; i <= mxd[t]; ++i) {
                dct = bin[t][i];
                bin[t][i] = cnt;
                cnt += dct;
            }
            bin[t][mxd[t]+1] = n[t];
        }

        for (uint32_t t = 0; t <= 1; ++t)
        for (uint32_t u = 0; u < n[t]; ++u) {
            uint32_t d = deg[t][u];
            uint32_t i = bin[t][d]++;
            vqueue[t][i] = u;
            vpos[t][u] = i; 
        }


        // uint32_t kk = 0;
        // uint32_t ss[2] = {0};
        // uint32_t it = 0, cnts = 0;
        // while (ss[0] < n[0] || ss[1] < n[1]) {
        //     uint32_t t = 0, z = 1;
        //     while (ss[t] < bin[t][kk] || ss[z] < bin[z][kk]){
        //         for (t = 0; t <= 1; t++) {
        //             z = t^1;
        //             for (uint32_t i = ss[t]; i < bin[t][kk]; ++i) {
        //                 uint32_t u = vqueue[t][i];
        //                 cores[t][u] = kk;
        //                 colorOrdering[t][i] = cnts++;
        //                 assert(u < n[t]);
        //                 for (uint32_t j = p[t][u]; j < p[t][u+1]; ++j) {
        //                     uint32_t v = e[t][j];
        //                     if (deg[z][v] > kk) {
        //                         uint32_t dz = --deg[z][v];
        //                         uint32_t pz = bin[z][dz]++;
        //                         uint32_t vp = vpos[z][v];
        //                         vqueue[z][vp] = vqueue[z][pz];
        //                         vqueue[z][pz] = v;
        //                         vpos[z][v] = pz;
        //                         vpos[z][vqueue[z][vp]] = vp;
        //                     }
        //                 }
        //             }
        //             ss[t] = bin[t][kk];
        //             // if (it < 10)
        //             // printf("it=%d, kk=%d, ss[0]=%d, bin[0]=%d, ss[1]=%d, bin[1]=%d\n", it, kk, ss[0], bin[0][kk], ss[1], bin[1][kk]);
        //         }
        //         maxcore = std::max(kk,maxcore);
        //     }   
        //     kk++;
        //     it ++;
        // }
        // // maxcore = kk - 1;
        // printf("ss[0]=%d, bin[0]=%d, ss[1]=%d, bin[1]=%d\n", ss[0], bin[0][kk], ss[1], bin[1][kk]);
        // assert(ss[0] == bin[0][kk]);
        // assert(ss[1] == bin[1][kk]);
        // printf("maxcore=%d\n", maxcore);

        for (uint32_t t = 0; t <= 1; ++t) {
            bin[t].resize(n[t], 0);
            for (uint32_t u = 0; u < n[t]; ++u) {
                uint32_t up = vpos[t][u];
                realLable[t][up] = _lable[t][u];
                cores[t][up] = deg[t][u];
                bin[t][up] = p[t][u+1]-p[t][u];
            }
            p[t][0] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                p[t][u+1] = p[t][u] + bin[t][u];
            }
        }
     
        for (uint32_t i = 0; i < m; ++i) {
            uint32_t u = edges[i].u;
            uint32_t v = edges[i].v;
            edges[i].u = vpos[0][u];
            edges[i].v = vpos[1][v];

            e[0][p[0][edges[i].u]++] = edges[i].v;
            e[1][p[1][edges[i].v]++] = edges[i].u;
        }

        for (uint32_t t = 0; t <= 1; ++t) {
            uint32_t z = t^1;
            p[t][0] = 0;
            for (uint32_t u = 0; u < n[t]; ++u) {
                p[t][u + 1] = bin[t][u] + p[t][u];
                std::sort(e[t].begin()+p[t][u], e[t].begin()+p[t][u+1]);
            }
        }
        // exit(1);
    }

    void relableCore_old() {
        uint32_t n = n1+n2;
        uint32_t maxd = std::max(maxDu, maxDv);
        std::vector<uint32_t> deg, vpos, bin, vqueue;
        deg.resize(n);
        vpos.resize(n);
        bin.resize(n, 0);
        vqueue.resize(n);
        cores[0].resize(n1);
        cores[1].resize(n2);

        uint32_t t = 0;
        for (uint32_t u = 0; u < n1; ++u) {
            uint32_t d = p[t][u+1] - p[t][u];
            deg[u] = d;
            bin[d]++;
        }
        t = t^1;
        for (uint32_t u = 0; u < n2; ++u) {
            uint32_t d = p[t][u+1] - p[t][u];
            deg[u+n1] = d;
            bin[d]++;
        }
        uint32_t s = 0;
        for (uint32_t i = 0; i <= maxd; ++i) {
            uint32_t d = bin[i];
            bin[i] = s;
            s += d;
        }

        for (uint32_t u = 0; u < n; ++u) {
            uint32_t d = deg[u];
            uint32_t p = bin[d]++;
            vqueue[p] = u;
            vpos[u] = p;
        }
        bin[maxd+1] = n;

        for (uint32_t i = maxd; i >= 1; --i) {
            bin[i] = bin[i-1];
        }
        bin[0] = 0;

        uint32_t kk = 0;
        s = 0;
        while (s < n) {
            if (bin[kk] >= n) break;
            uint32_t u, v, t, d, vp;
            for (uint32_t i = bin[kk]; i < bin[kk+1]; ++i) {
                u = vqueue[i];
                t = 0; s++;
                if (u >= n1) { u-= n1; t = 1;}

                cores[t][u] = kk;

                for (uint32_t j = p[t][u]; j < p[t][u+1]; ++j) {
                    v = e[t][j];
                    if (t == 0) v += n1;
                    d = deg[v] --;
                    if (d > kk) {
                        vp = vpos[v];
                        std::swap(vqueue[vp], vqueue[bin[d]]);
                        vpos[vqueue[vp]] = vp;
                        vpos[v] = bin[d]++;

                    }
                }
                maxcore = std::max(kk,maxcore);
            }

            kk++;
        }
        
        printf("relabel done, maxcore=%d\n", maxcore);

        s = 0;
        printf("check\n");
        std::vector<uint32_t> check(n, 0);
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t u = vqueue[i];
            if (check[u] <= 0) {
                s++; check[u]++;
            }
            else assert(check[u] == 1);

            assert(i == vqueue[vpos[i]]);
        }
        assert(s == n);

        uint32_t nn1 = 0, nn2 = 0;
        for (uint32_t i = 0; i < n; i++) {
            uint32_t x = vqueue[i];
            if (x < n1) {
                // if (cores[0][x] < sz[1]) {
                //     vpos[x] = n+1; vqueue[i] = n+1;
                //     continue;
                // }
                // if (x == 61340) {printf("realu=61340, u=%d\n", nn1); exit(1);}
                vpos[x] = nn1;
                vqueue[i] = nn1;
                deg[nn1] = cores[0][x];
                realLable[0][nn1++] = x;
            }
        }
        for (uint32_t i = 0; i < n; i++) {
            uint32_t x = vqueue[i];
            if (x >= n1 && x < n) {
                // if (cores[1][x-n1] < sz[0]) {
                //     vpos[x] = n+1; vqueue[i] = n+1;
                //     continue;
                // }
                vpos[x] = nn2;
                vqueue[i] = nn2+nn1;
                deg[nn2+n1] = cores[1][x-n1];
                realLable[1][nn2++] = x-n1;
            }
        }
        for (uint32_t u = 0; u < n; u++) {
            if (u < n1) cores[0][u] = deg[u];
            else cores[1][u-n1] = deg[u];
            deg[u] = 0;
        }
        printf("nn1=%d, nn2=%d, deg[nn1+nn2]=%d\n", nn1,nn2, deg[nn1+nn2]);

        for (uint32_t i = 0; i < n; ++i) assert(deg[i] == 0);
        for (long i = 0; i < m; ++i) {
            if (vpos[edges[i].u] > n || vpos[edges[i].v+n1] > n) {
                edges[i].u = n1; edges[i].v = n2;
                continue;
            }
            edges[i].u = vpos[edges[i].u];
            edges[i].v = vpos[edges[i].v+n1];
            assert(edges[i].u < nn1);
            assert(edges[i].v < nn2);
            deg[edges[i].u]++;
            deg[edges[i].v+n1]++;
        }
        p[0][0] = 0; maxDu = 0;
        for(uint32_t u = 0; u < n1; u++) {
            p[0][u + 1] = deg[u] + p[0][u];
            maxDu = std::max(deg[u], maxDu);
        }
        for(uint32_t i = 0; i < m; i++) {
            if (edges[i].u >= n1 || edges[i].v >= n2) continue;
            assert(edges[i].u < nn1);
            assert(edges[i].v < nn2);
            e[0][p[0][edges[i].u]++] = edges[i].v; 
        }
        p[0][0] = 0;
        for(uint32_t u = 0; u < nn1; u++) {
            p[0][u + 1] = deg[u] + p[0][u];
            std::sort(e[0].begin()+p[0][u], e[0].begin()+p[0][u+1]);
        }

        p[1][0] = 0; maxDv = 0;
        for(uint32_t v = 0; v < n2; v++) {
            p[1][v + 1] = deg[v+n1] + p[1][v];
            maxDv = std::max(deg[v+n1], maxDv);
        }
        for(uint32_t i = 0; i < m; i++) {
            if (edges[i].u >= n1 || edges[i].v >= n2) continue;
            e[1][p[1][edges[i].v]++] = edges[i].u; 
        }
        p[1][0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            p[1][v + 1] = deg[v+n1] + p[1][v];
            std::sort(e[1].begin()+p[1][v], e[1].begin()+p[1][v+1]);
        }

        // exit(1);

        n1 = nn1; n2 = nn2;
        this->n[0] = n1; this->n[1] = n2;
        this->mxd[0] = maxDu;
        this->mxd[1] = maxDv;
    }


    void coloring2hop(){

        std::vector<uint32_t> colorcnts, vqueue;
        colors[0].resize(n[0]);
        colors[1].resize(n[1]);
        colorcnts.resize(n[0]+n[1]+1);
        vqueue.resize(n[0]+n[1]+1);


        uint32_t cnrmxcl = 0, t = 0, z = 1;

        if (mxd[0] < mxd[1]) { t = 1; z = 0;}

        for (uint32_t u = n[t]-1; u >= 0; u--) {

            cnrmxcl = 0;
            for (uint32_t i = p[t][u]; i < p[t][u+1]; i++) {
                uint32_t v = e[t][i];
                for (uint32_t j = p[z][v+1]-1; j >= p[z][v]; j--) {
                    uint32_t w = e[z][j];
                    uint32_t c = colors[t][w];
                    cnrmxcl = std::max(c, cnrmxcl);
                    colorcnts[c] = 1;
                    if (w < u || j == 0) break;
                }
            }

            uint32_t setc = 0;
            for (uint32_t i = 0; i <= cnrmxcl; ++i) {
                if (i > 0 && setc == 0 && colorcnts[i] == 0) setc = i;
                colorcnts[i] = 0;
            }
            if (setc == 0) {
                colors[t][u] = cnrmxcl+1;
                mxcl = std::max(cnrmxcl+1, mxcl);
            }
            else colors[t][u] = setc;

            if (u == 0) break;
        }
        for (uint32_t u = 0; u < n[t]; u++) {
            assert(colors[t][u] > 0);
            assert(colors[t][u] <= mxcl);
        }
        printf("n[t]=%d, mxcl=%d\n", n[t], mxcl);
        maxcols[t] = mxcl; 
        maxcols[z] = 0;
        // for (uint32_t u = n[z]-1; u >= 0; u--) {

        //     cnrmxcl = 0;
        //     for (uint32_t i = p[z][u]; i < p[z][u+1]; i++) {
        //         uint32_t v = e[z][i];
        //         uint32_t c = colors[t][v];
        //         cnrmxcl = std::max(c, cnrmxcl);
        //         colorcnts[c] = 1;
        //         for (uint32_t j = p[t][v+1]-1; j >= p[t][v]; j--) {
        //             uint32_t w = e[t][j];
        //             c = colors[z][w];
        //             cnrmxcl = std::max(c, cnrmxcl);
        //             colorcnts[c] = 1;
        //             if (w < u || j == 0) break;
        //         }
        //     }

        //     uint32_t setc = 0;
        //     for (uint32_t i = 0; i <= cnrmxcl; ++i) {
        //         if (i > 0 && setc == 0 && colorcnts[i] == 0) setc = i;
        //         colorcnts[i] = 0;
        //     }
        //     if (setc == 0) {
        //         colors[z][u] = cnrmxcl+1;
        //         mxcl = std::max(cnrmxcl+1, mxcl);
        //     }
        //     else colors[z][u] = setc;

        //     if (u == 0) break;
        // }
        // for (uint32_t u = 0; u < n[z]; u++) {
        //     assert(colors[z][u] > 0);
        //     assert(colors[z][u] <= mxcl);
        // }
        // printf("n[z]=%d, mxcl=%d\n", n[z], mxcl);
        // maxcols[z] = mxcl;
    }

    uint32_t deg(uint32_t t, uint32_t u) {
        return p[t][u+1] - p[t][u];
    }

    bool is_nbr(uint32_t t, uint32_t u, uint32_t v) {
        return nbrIndex[t][u].find(v);
    }


    //     void coloring2hop(){

    //     std::vector<uint32_t> colorcnts, vqueue;
    //     colors[0].resize(n[0]);
    //     colors[1].resize(n[1]);
    //     colorcnts.resize(n[0]+n[1]+1);
    //     vqueue.resize(n[0]+n[1]+1);
    //     uint32_t cnrmxcl = 0, t = 0, z = 1;

    //     if (mxd[t] < mxd[z]) { z = t; t = z^1;}

    //     for (uint32_t u = n[t]-1; u >= 0; u--) {

    //         cnrmxcl = 0;
    //         for (uint32_t i = p[t][u]; i < p[t][u+1]; i++) {
    //             uint32_t v = e[t][i];
    //             for (uint32_t j = p[z][v+1]-1; j >= p[z][v]; j--) {
    //                 uint32_t w = e[z][j];
    //                 uint32_t c = colors[t][w];
    //                 cnrmxcl = std::max(c, cnrmxcl);
    //                 colorcnts[c] = 1;
    //                 if (w < u || j == 0) break;
    //             }
    //         }

    //         uint32_t setc = 0;
    //         for (uint32_t i = 0; i <= cnrmxcl; ++i) {
    //             if (i > 0 && setc == 0 && colorcnts[i] == 0) setc = i;
    //             colorcnts[i] = 0;
    //         }
    //         if (setc == 0) {
    //             colors[t][u] = cnrmxcl+1;
    //             mxcl = std::max(cnrmxcl+1, mxcl);
    //         }
    //         else colors[t][u] = setc;

    //         if (u == 0) break;
    //     }
    //     for (uint32_t u = 0; u < n[t]; u++) {
    //         assert(colors[t][u] > 0);
    //         assert(colors[t][u] <= mxcl);
    //     }
    //     printf("n[t]=%d, mxcl=%d\n", n[t], mxcl);
    //     maxcols[t] = mxcl;
    //     for (uint32_t u = n[z]-1; u >= 0; u--) {

    //         cnrmxcl = 0;
    //         for (uint32_t i = p[z][u]; i < p[z][u+1]; i++) {
    //             uint32_t v = e[z][i];
    //             uint32_t c = colors[t][v];
    //             cnrmxcl = std::max(c, cnrmxcl);
    //             colorcnts[c] = 1;
    //             for (uint32_t j = p[t][v+1]-1; j >= p[t][v]; j--) {
    //                 uint32_t w = e[t][j];
    //                 c = colors[z][w];
    //                 cnrmxcl = std::max(c, cnrmxcl);
    //                 colorcnts[c] = 1;
    //                 if (w < u || j == 0) break;
    //             }
    //         }

    //         uint32_t setc = 0;
    //         for (uint32_t i = 0; i <= cnrmxcl; ++i) {
    //             if (i > 0 && setc == 0 && colorcnts[i] == 0) setc = i;
    //             colorcnts[i] = 0;
    //         }
    //         if (setc == 0) {
    //             colors[z][u] = cnrmxcl+1;
    //             mxcl = std::max(cnrmxcl+1, mxcl);
    //         }
    //         else colors[z][u] = setc;

    //         if (u == 0) break;
    //     }
    //     for (uint32_t u = 0; u < n[z]; u++) {
    //         assert(colors[z][u] > 0);
    //         assert(colors[z][u] <= mxcl);
    //     }
    //     printf("n[z]=%d, mxcl=%d\n", n[z], mxcl);
    //     maxcols[z] = mxcl;
    // }
};

#endif