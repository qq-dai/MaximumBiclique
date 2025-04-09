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

    void branch_basic(uint32_t deep);

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
    void run_basic(int alg);
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


// #ifdef UB1_
//         uint32_t qsz[2] = {0};
//         for(uint32_t j = S[z].r; j < c[z]; j++) _deg[z][S[z][j]] = 0;
//         for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
//         // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
//             uint32_t u = S[t][i];
//             uint32_t e = 0;
//             if(g->deg(t, u) > 2*S[z].cSize()) {
//                 for(uint32_t j = S[z].r; j < S[z].c; j++) {
//                     uint32_t v = S[z][j];
//                     if(g->is_nbr(t, u, v)) {
//                         e++; _deg[z][v]++;
//                     }
//                 }
//             }
//             else {
//                 for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
//                     uint32_t v = g->e[t][j];
//                     uint32_t pv = S[z].pos(v);
//                     if(S[z].r <= pv && pv < S[z].c) {//in C
//                         e++; _deg[z][v]++;
//                     }
//                 }
//             }
            
//             _deg[t][u] = e;
//             if (S[z].r+e < sz[z] || S[t].c * (S[z].r+e) <= mSize) {
//                 S[t].swapByPos(--S[t].c, i--);
//                 vqueue[t][qsz[t]++] = u;
//                 if (S[t].c < sz[t] || S[t].c * S[z].c <= mSize)  {
//                     S[t].c = c[t];
//                     S[z].c = c[z];
//                     return;
//                 }
//             }
//         }

//         for(int j = S[z].r; j < S[z].c; j++) {
//             uint32_t u = S[z][j];
//             uint32_t e = _deg[z][u];
//             if (r[t]+e < sz[t] || S[z].c *(r[t]+e) < mSize) {
//                 S[z].swapByPos(--S[z].c, j--);
//                 vqueue[z][qsz[z]++] = u;
//             }
//         }

//         // if (deep == 1 && S[z][0] == 26285) {
//         //     printf("sz[t]=%d, sz[z]=%d\n", sz[t], sz[z]);
//         //     for (uint32_t i = S[t].r; i < S[t].c; ++i) {
//         //         printf("t=%d, v=%d, deg[v]=%d, _deg[v]=%d\n", t, S[t][i], _deg[t][S[t][i]]+S[z].r, S[z].c - S[z].r - _deg[t][S[t][i]]);
//         //     }
//         // }

//         uint32_t qr[2] = {0};
//         while(qr[0] < qsz[0] || qr[1] < qsz[1]) {
//             for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
//                 uint32_t zzz = ttt^1;
//                 for (uint32_t i = qr[ttt]; i < qsz[ttt]; ++i) {
//                     uint32_t u = vqueue[ttt][i];
//                     _deg[ttt][u] = 0;

//                     if(g->deg(ttt, u) > 10*S[zzz].cSize()) {
//                         for(int j = S[zzz].r; j < S[zzz].c; j++) {
//                             uint32_t v = S[zzz][j];
//                             if(g->is_nbr(ttt, u, v)) {
//                                 uint32_t d = --_deg[zzz][v];
//                                 d += S[ttt].r;
//                                 if (d < sz[ttt] || d * S[zzz].c <= mSize) {
//                                     S[zzz].swapByPos(--S[zzz].c, j--);
//                                     vqueue[zzz][qsz[zzz]++] = v;
//                                     if (S[t].c < sz[t] && S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  break;
//                                 }
//                             }
//                         }
//                     }
//                     else 
//                     for (uint32_t j = g->p[ttt][u]; j < g->p[ttt][u+1]; ++j) {
//                         uint32_t v = g->e[ttt][j];
//                         uint32_t pv = S[zzz].pos(v);
//                         if (pv < S[zzz].c && S[zzz].r <= pv) {
//                             uint32_t d = --_deg[zzz][v];
//                             d += S[ttt].r;
//                             if (d < sz[ttt] || d * S[zzz].c <= mSize) {
//                                 S[zzz].swapByPos(--S[zzz].c, pv);
//                                 vqueue[zzz][qsz[zzz]++] = v;
//                                 if (S[t].c < sz[t] && S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  break;
//                             }
//                         }
//                     }
//                     if (S[t].c < sz[t] && S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  {
//                         i = qsz[ttt]; qr[zzz] = qsz[zzz];
//                         break;
//                     }
//                 }
//                 qr[ttt] = qsz[ttt];

//             }
//         }

//         // if (sz[gt] != uint32_t(mSize/sz[gt^1]/2)) { printf("sz[gt]=%d, sz[gt^1]=%d, mSize=%d\n", sz[gt], sz[gt^1], mSize); }
//         // assert(sz[gt] == uint32_t(mSize/sz[gt^1]/2) );

//         maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
//         uint32_t minI = 0, minU = 0, mnume = 0, minE = std::max(S[t].c, S[z].c)+1;

//         mnume = S[t].r*S[z].c+S[z].r*(S[t].c-S[t].r);
//         if (S[t].cSize() > S[z].cSize()) {
//             t = z; z = t^1;
//         }

//         qsz[t] = 0; qsz[z] = 0;
//         for(uint32_t j = S[t].r; j < S[t].c; j++) {
//             uint32_t u = S[t][j];
//             uint32_t e = _deg[t][u];

//             if (e == S[z].c - r[z])  {
//                 vqueue[t][j] = vqueue[t][S[t].r];
//                 S[t].swapByPos(S[t].r++, j);
//                 continue;
//             }
//             else vqueue[t][j] = S[z].c - r[z] - e;


//             if(e > maxE[t]) {
//                 maxE[t] = e;
//                 maxI[t] = j;
//                 maxUV[t] = u;
//             }
//             if (e < minE) {
//                 minE = e;
//                 minI = j;
//                 minU = u;
//             }
//             mnume += e;
//         }

//         for(int j = S[z].r; j < S[z].c; j++) {
//             uint32_t u = S[z][j];
//             uint32_t e = _deg[z][u]; 

//             if (e == S[t].c - r[t]) {
//                 vqueue[z][j] = vqueue[z][S[z].r];
//                 S[z].swapByPos(S[z].r++, j);
//                 continue;
//             }else vqueue[z][j] = S[t].c - r[t] - e;

//             if(e > maxE[z]) {
//                 maxE[z] = e;
//                 maxI[z] = j;
//                 maxUV[z] = u;
//             }
//         }

//         if (S[t].c < sz[t] && S[z].c < sz[z]|| S[0].c * S[1].c <= mSize)  {
//             S[t].c = c[t]; S[t].r = r[t]; 
//             S[z].c = c[z]; S[z].r = r[z];
//             return;
//         }

//         if (S[t].cSize() > 0 || S[z].cSize() > 0){

//             // // test
//             // uint32_t cntd[2] = {0};
//             // for (uint32_t i = S[t].r; i < S[t].c; ++i)
//             //     cntd[t] += vqueue[t][i];
//             // for (uint32_t i = S[z].r; i < S[z].c; ++i)
//             //     cntd[z] += vqueue[z][i];
//             // uint32_t cnt = 0;
//             // for (uint32_t i = S[t].r; i < S[t].c; ++i) {
//             //     uint32_t v = S[t][i];
//             //     uint32_t d = 0;
//             //     for (uint32_t j = S[z].r; j < S[z].c; ++j)
//             //         if (g->is_nbr(t,v,S[z][j])) d++;
//             //     cnt += (S[z].cSize() - d);
//             // }
//             // if (cntd[t] != cnt || cntd[z] != cnt) {
//             //     S[t].print();
//             //     S[z].print();
//             //     printf("deep=%d, cntd[t]=%d, cntd[z]=%d\n", deep, cntd[t], cntd[z]);

//             //     printf("correct cntd=%d\n", cnt); 
//             //     for (uint32_t i = S[t].r; i < S[t].c; ++i) {
//             //         uint32_t v = S[t][i];
//             //         int d = 0;
//             //         for (uint32_t j = S[z].r; j < S[z].c; ++j)
//             //             if (g->is_nbr(t,v,S[z][j])) d++;
//             //         printf("v =%d, d=%d, _d=%d, vqueue[v]=%d\n", v, d+S[z].r, S[z].cSize()-d, vqueue[t][i]);
//             //     }
//             //     exit(1);
//             // }
//             // assert(cntd[t] == cntd[z]);
//             // //end test

//             // std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c, [](int a, int b){return a > b;});
//             // std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c, [](int a, int b){return a > b;});

//             // std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c);
//             // std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c);


//             uint32_t numubs = 0, ub[2] = {0};
//             uint32_t qsz[2] = {0};
//             for (int i = S[t].r; i < S[t].c; ++i) {
//                 uint32_t v = S[t][i];
//                 uint32_t c = g->colors[t][v];
//                 if (bin[t][c] == 0) ub[t]++;
//                 if (bin[t][c] == 1) vqueue[t][qsz[t]++] = c;
//                 bin[t][c]++;
//             }

//             if (S[t].r + ub[t] < sz[t]) {

//                 for (int i = S[t].r; i < S[t].c; ++i) bin[t][g->colors[t][S[t][i]]] = 0;

//                 S[t].c = c[t]; S[t].r = r[t]; 
//                 S[z].c = c[z]; S[z].r = r[z];
//                 return;
//             }

//             for (int i = S[z].r; i < S[z].c; ++i) {
//                 uint32_t v = S[z][i];
//                 uint32_t c = g->colors[z][v];
//                 if (bin[z][c] == 0) ub[z]++;
//                 if (bin[z][c] == 1) vqueue[z][qsz[z]++] = c;
//                 bin[z][c]++;
//             }

//             if (S[z].r + ub[z] < sz[z]) {

//                 for (int i = S[z].r; i < S[z].c; ++i) bin[z][g->colors[z][S[z][i]]] = 0;

//                 S[t].c = c[t]; S[t].r = r[t]; 
//                 S[z].c = c[z]; S[z].r = r[z];
//                 return;
//             }

//             // for (int i = S[t].r; i < S[t].c; ++i) {
//             //     uint32_t _d = vqueue[t][i];
//             //     uint32_t d  =  S[z].c - _d;
//             //     if ((i+1) > sz[t])
//             //         numubs = std::max(numubs, (i+1) * d);
//             // }
//             // for (uint32_t i = S[z].r; i < S[z].c; ++i) {
//             //     uint32_t _d = vqueue[z][i];
//             //     uint32_t d  =  S[t].c - _d;
//             //     if ((i+1) > sz[z])
//             //         numubs = std::max(numubs, (i+1) * d);
//             // }
//             // if (numubs <= mSize) {
//             //     S[t].c = c[t]; S[t].r = r[t]; 
//             //     S[z].c = c[z]; S[z].r = r[z];
//             //     return;
//             // } 



//             // uint32_t head[2] = {0}, tail[2] = {0},
//             // head[t] = S[t].r; head[z] = S[z].r;
//             // tail[t] = S[t].c; tail[z] = S[z].c;
//             // ub[t] = S[t].c; ub[z] = S[z].c;

//             // if (S[z][0] == 10622 && deep == 0) {
//             //     printf("deep=%d, head[t]=%d, tail[t]=%d, cntd[t]=%d\n", deep, head[t], tail[t], cntd[t]);
//             //     printf("deep=%d, head[z]=%d, tail[z]=%d, cntd[z]=%d\n", deep, head[z], tail[z], cntd[z]);

//             //     for (uint32_t i = S[t].r; i < S[t].c; ++i) {
//             //         printf("v=%d, _d=%d\n", S[t][i], vqueue[t][i]);
//             //     }
//             //     printf("Z\n");
//             //     for (uint32_t i = S[z].r; i < S[z].c; ++i) {
//             //         printf("i=%d, v=%d, _d=%d\n", i, S[z][i], vqueue[z][i]);
//             //     }
//             // }

//             // while(true) {
//             //     if (head[t] >= tail[t] || head[z] >= tail[z] || vqueue[t][head[t]] <= 0 || vqueue[z][head[z]] <= 0) break;

//             //     uint32_t d[2] = {0};
//             //     d[t] = S[z].c - head[z] - vqueue[t][head[t]] + S[z].r;
//             //     d[z] = S[t].c - head[t] - vqueue[z][head[z]] + S[t].r;

//             //     if (d[z] < 0) d[z] = 0;

//             //     if (S[z][0] == 10622 && deep == 0) {
//             //         printf("d[t]=%d, head[t]=%d, tail[t]=%d, vqueue[t][head[t]]=%d\n", d[t], head[t], tail[t], vqueue[t][head[t]]);
//             //         printf("d[z]=%d, head[z]=%d, tail[z]=%d, vqueue[z][head[z]]=%d\n", d[z], head[z], tail[z], vqueue[z][head[z]]);
//             //         printf("ub[t]=%d, ub[z]=%d, ubSize=%d\n", ub[t], ub[z], ub[t] *ub[z]);
//             //     }

//             //     if (d[t] == S[z].r) {

//             //         if (S[z].r * (S[t].c - head[t] + S[t].r) >= mSize) {break;}
//             //         d[z] = d[t]+1;
//             //     }
//             //     else if (d[z] == S[t].r) {

//             //         if (S[t].r * (S[z].c - head[z] + S[z].r) >= mSize) {break;}
//             //         d[t] = d[z]+1;
//             //     }

//             //     if (d[t] <= d[z]) {
//             //         uint32_t len = tail[z], _d = vqueue[t][head[t]];

//             //         if (S[z][0] == 10622 && deep == 0) {
//             //             printf("t,head[t]=%d, tail[t]=%d, _d=%d, vqueue[t][head[t]]=%d\n", head[t], tail[t], _d, vqueue[t][head[t]]);
//             //             printf("t,tail[z]=%d, tail[z] - _d=%d\n", tail[z], tail[z]-_d);

//             //         }

//             //         for (int j = len-1; j >= len - _d; --j) {
//             //             assert(vqueue[z][j] > 0);
//             //             if (j < head[z]) {
//             //                 S[t].print();
//             //                 S[z].print();
//             //                 printf("deep=%d\n", deep);
//             //             }
//             //             assert(j >= head[z]);
//             //             uint32_t dv = vqueue[z][j]--;
//             //             if (dv == 1) tail[z] = j;
//             //         }
//             //         head[t]++; ub[t]--;
//             //     }
//             //     else {
//             //         uint32_t len = tail[t], _d = vqueue[z][head[z]];
//             //         if (S[z][0] == 10622 && deep == 0) {
//             //             printf("head[z]=%d, tail[z]=%d, _d=%d, vqueue[z][head[z]]=%d\n", head[z], tail[z], _d, vqueue[z][head[z]]);
//             //             printf("tail[t]=%d, tail[t] - _d=%d\n", tail[t], tail[t]-_d);

//             //         }
//             //         for (uint32_t j = len-1; j >=  std::max((uint32_t) 0, (len - _d)); --j) {
//             //             if (j < 0) {
//             //                 S[t].print();
//             //                 S[z].print();
//             //                 printf("deep=%d\n", deep);
//             //                 printf("tail[t]=%d, tail[t] - _d=%d, j = %d\n", len, len -_d, j);
//             //             }
//             //             assert(vqueue[t][j] > 0);
//             //             if (j < head[t]) {
//             //                 S[t].print();
//             //                 S[z].print();
//             //                 printf("deep=%d\n", deep);
//             //             }
//             //             assert(j >= head[t]);
//             //             uint32_t dv = vqueue[t][j]--;
//             //             if (S[z][0] == 10622 && deep == 0) {
//             //                 printf("dv=%d, j=%d\n", dv, j);
//             //             }
//             //             if (dv == 1) tail[t] = j;
//             //             if (j == 0) break;
//             //         }
//             //         head[z]++; ub[z]--;
//             //     }

//             //     if (S[z][0] == 10622 && deep == 0)
//             //     for (uint32_t i = S[t].r; i < S[t].c; ++i) {
//             //         printf("v=%d, _d=%d\n", S[t][i], vqueue[t][i]);
//             //     }

//             //     if (ub[t] * ub[z] <= mSize) break;
//             // }
//             // if (S[z][0] == 23948 && deep == 0) exit(1);

//             // if (ub[t] * ub[z] <= mSize) {
//             //     S[t].c = c[t]; S[t].r = r[t]; 
//             //     S[z].c = c[z]; S[z].r = r[z];
//             //     return;
//             // }
//         }

//         // if (mSize > 0 && S[gt^1].c > 0 && sz[gt] == uint32_t(mSize/sz[gt^1]/stp) && uint32_t(mSize / S[gt^1].c) > sz[gt]) {
//         //     // subbcnt++;
//         //     uint32_t oldsz[2] = {0}, Szc = S[gt^1].c;
//         //     oldsz[0] = sz[0]; oldsz[1] = sz[1];
//         //     // printf("in, S[gt][0]=%d, deep=%d, gt=%d, t=%d, sz[t]=%d, sz[z]=%d, S[gt].c=%d, S[gt^1].c=%d, mSize=%d\n", S[gt][0], deep, gt, t, sz[t], sz[z], S[gt].c, S[gt^1].c, mSize);
//         //     sz[gt] = oldsz[gt]; sz[gt^1] = Szc;
//         //     branchNew(deep, ne);
//         //     // printf("in1, deep=%d, gt=%d, t=%d, sz[t]=%d, sz[z]=%d, S[gt].c=%d, S[gt^1].c=%d\n", deep, gt, t, sz[t], sz[z], S[gt].c, S[gt^1].c);
//         //     S[t].c = c[t]; S[z].c = c[z];
//         //     sz[gt] = uint32_t(mSize/Szc); sz[gt^1] = oldsz[gt^1];
//         //     branchNew(deep+1, ne);
//         //     // printf("out, deep=%d, gt=%d, t=%d, sz[t]=%d, sz[z]=%d, S[gt].c=%d, S[gt^1].c=%d\n", deep, gt, t, sz[t], sz[z], S[gt].c, S[gt^1].c);
//         //     // exit(1);
//         //     sz[0] = oldsz[0]; sz[1] = oldsz[1];
//         //     S[t].c = c[t]; S[z].c = c[z];
//         //     return;
//         // }

// #else
//     // r[0] = S[0].r; r[1] = S[1].r;
//     // c[0] = S[0].c; c[1] = S[1].c;

//     // sube = 0;

//     maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
//     uint32_t minI = 0, minU = 0, minE = S[z].c+1;
//     for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
//     // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
//         uint32_t u = S[t][i];
//         uint32_t e = 0;
//         if(g->deg(t, u) > 2*S[z].cSize()) {
//             for(uint32_t j = S[z].r; j < S[z].c; j++) {
//                 uint32_t v = S[z][j];
//                 if(g->is_nbr(t, u, v)) {
//                     e++; _deg[z][v]++;
//                 }
//             }
//         }
//         else {
//             for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
//                 uint32_t v = g->e[t][j];
//                 uint32_t pv = S[z].pos(v);
//                 if(S[z].r <= pv && pv < S[z].c) {//in C
//                     e++; _deg[z][v]++;
//                 }
//             }
//         }

//         assert(e == _deg[t][u]);

//         // if (e == S[z].cSize()) {
//         //     S[t].swapByPos(S[t].r++, i);
//         //     continue;
//         // }

//         if(e > maxE[t]) {
//             maxE[t] = e;
//             maxI[t] = i;
//             maxUV[t] = u;
//             // if (e >= S[z].cSize()) break;
//         }
//         if (e < minE) {
//             minE = e;
//             minI = i;
//             minU = u;
//         }
//         // _deg[t][u] = e;
//         // sube += e;
//     }

//     for(int i = S[z].r; i < S[z].c; i++) {
//         uint32_t u = S[z][i];
//         uint32_t e = _deg[z][u];
//         _deg[z][u] = 0;
//         // if (e == S[t].c - r[t]) {
//         //     // uint32_t d = 0;
//         //     // for (uint32_t j = S[t].r; j  < S[t].c; ++j)
//         //     //     if (g->is_nbr(z,u, S[t][j])) d++;
//         //     // assert(d == e);
//         //     // assert(i == S[z].pos(u));
//         //     S[z].swapByPos(S[z].r++, i);
//         // }
//         // else 
//         if(e > maxE[z]) {
//             maxE[z] = e;
//             maxI[z] = i;
//             maxUV[z] = u;
//         }
//     }
// #endif