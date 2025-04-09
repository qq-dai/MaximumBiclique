#include "MBCS.h"

#define UB1_
#define RUN_1

#define _BASIC_

#ifdef RUN_1

void MBCS::run() {

    printf("run\n");
    // for(uint32_t i = 0; i < g->n[1]; ++i) {
    //     uint32_t d = g->p[1][i+1] - g->p[1][i];
    //     if (d >= 7591)
    //         printf("V v=%d, d=%d, is_nbr(12,53)=%d\n", i, d, g->is_nbr(0,12,53));
    // }
    S[0].c = g->n[0]; S[1].c = g->n[1];

    // mSize=15255;
    // mSize= 205180;
    // mSize = 880;

    initSize();

    // mSize = 30000;

    // if (mSize >= std::min(g->maxDu, g->maxDu) * 100){
    //     #define BIBRANCH_
    // }

    uint32_t qsz[2] = {0};

    stp = 1.5;

    // return;

    // sz[0] = 7; sz[1] = 13;

    // printf("g->n[0]=%d, g->n[1]=%d, mSize=%d\n", g->n[0], g->n[1], mSize);
    // return;
    uint32_t t = 0, z = 1;
    if (g->maxDu < g->maxDv) {t = 1; z = 0;} 
    printf("g->n[0]=%d, g->n[1]=%d, mSize=%d, t=%d, z=%d, stp=%.3f\n", g->n[0], g->n[1], mSize, t, z, stp);
    // printf("t=%d, u=21878, realu=%d\n", t, g->realLable[t][21878]); return;
 
    unsigned long maxBranchduration = 0;
    uint32_t crsz[2] = {0}, tt = 0, zz = 1;

    std::vector<uint32_t> colorcnts;
    colorcnts.resize(g->n[0]+g->n[1]+1);
    eside = t;

    // tt = z; zz = tt^1;
    tt = t; zz = tt^1; gt = tt;

    crsz[0] = sz[0]; crsz[1] = sz[1];

    // sz[tt] = g->mxd[tt^1];
    // sz[tt] = crsz[tt]; sz[zz] = g->mxd[tt] / 4;
    assert(crsz[0] > 0);
    assert(crsz[1] > 0);
    // assert(sz[t] <= g->n[t]);

    // if (sz[tt] * sz[zz] > mSize) {sz[zz] = std::max(crsz[zz], mSize / sz[tt]);}

    uint32_t it = 0;

    auto t1 = std::chrono::steady_clock::now();

    for (uint32_t u = std::max((uint32_t)0, g->n[t]-crsz[t]); u >= 0; u--) {

    // for (uint32_t u = 61289; u < 61289+1 && u < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); u++) {
    // for (uint32_t u = 0; u < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); u++) {

        S[t].r = S[t].c = 0;
        S[z].r = S[z].c = 0;
        // if (S[t].capacity() <= u) {
            // printf("S[t].capacity()=%d, u=%d, std::max((uint32_t)0, g->n[t]-crsz[t])=%d\n", S[t].capacity(), u, std::max((uint32_t)0, g->n[t]-crsz[t]));
        // }

        assert(S[t].capacity() > u);
        assert(S[t][S[t].pos(u)] == u);
        S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

        uint32_t Szs = g->deg(t,u);
        uint32_t tc = 0, zc = 0;

        if (Szs * g->mxd[z] <= mSize) continue;

        uint32_t cnrmxcl = 0;
        for (uint32_t i = g->p[t][u]; i < g->p[t][u+1]; ++i) {
            uint32_t v = g->e[t][i];
            assert(S[z][S[z].pos(v)] == v);

            int nbrv = 0;

            for (int j = g->p[z][v+1]-1; j >= g->p[z][v]; j--) {
                uint32_t w = g->e[z][j];
                if (w < u) break;
                nbrv ++;
                if (j == 0) break;
            }

            if (Szs * nbrv <= mSize) {
                Szs -- ;
                continue;
            }

            S[z].swapByPos(S[z].c++, S[z].pos(v));
            sdeg[z][v] = nbrv;

            for (int j = g->p[z][v+1]-1; j >= g->p[z][v]; j--) {
                uint32_t w = g->e[z][j];
                uint32_t wp = S[t].pos(w);
                // assert(wp >= 0);
                // assert(wp < g->n[t]);
                // assert(S[t][wp] == w);
                if (w < u) break;
                if (wp >= S[t].c)
                    S[t].swapByPos(S[t].c++, wp);
                sdeg[t][w]++;
                // uint32_t c = g->colors[t][w];
                // cnrmxcl = std::max(c, cnrmxcl);
                // colorcnts[c] ++;
                if (j == 0) break;
            }
            // uint32_t setc = 0;
            // for (uint32_t i = 0; i <= cnrmxcl; ++i) {
            //     if (i > 0 && setc == 0 && colorcnts[i] == 0) {setc = i; break;}
            // }
            // for (int j = g->p[z][v+1]-1; j >= g->p[z][v]; j--) {
            //     uint32_t w = g->e[z][j];
            //     uint32_t wp = S[t].pos(w);
            //     if (w < u) break;
            //     // uint32_t c = g->colors[t][w];
            //     // colorcnts[c] --;
            //     if (j == 0) break;
            // }
            // if (setc == 0) {
            //     g->colors[z][v] = cnrmxcl+1;
            //     cnrmxcl += 1; setc = cnrmxcl;
            // }
            // else g->colors[z][v] = setc;
            // colorcnts[setc] ++;

            // if (u == 8329) {
            //     printf("v=%d, colors[v]=%d\n", v, g->colors[z][v]);
            // }
        }

        // for (uint32_t i = g->p[t][u]; i < g->p[t][u+1]; ++i) colorcnts[g->colors[z][g->e[t][i]]] = 0;
        // g->maxcols[z] = cnrmxcl; 
        // // if (cnrmxcl <= S[z].c) {
        // //     printf("u=%d, cnrmxcl=%d, S[z].c=%d\n", u, cnrmxcl, S[z].c);
        // // }
        // assert (cnrmxcl >= S[z].c);

        tc = S[t].c; zc = S[z].c;

        // reduction

        // sz[tt] = g->mxd[tt^1];
        sz[tt] = crsz[tt]; sz[zz] = g->mxd[tt] / 4;

        // if (sz[tt] * sz[zz] > mSize) {sz[zz] = std::max(crsz[zz], mSize / sz[tt] );}

        sz[zz] = std::max(crsz[zz], mSize / sz[tt]);

    
        it = 0;
        uint32_t old_numbranches = numbranches;
        while (sz[tt] <= g->mxd[zz] && sz[zz] >= crsz[zz]) {

            if (it > 0 && it < 2) {
                sz[tt] = std::min(g->mxd[zz], sz[tt] +1);
                sz[zz] = std::max(crsz[zz], mSize / sz[tt]);
            }
            else if  (it == 2) {
                sz[tt] = std::min(g->mxd[zz], sz[tt] +1);
                sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt] / stp));
            }
            else if (it > 2) {
                sz[tt] = std::min(g->mxd[zz], uint32_t(sz[tt] * stp));
                sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt] / stp));
            }

            old_numbranches = numbranches;

            S[t].c = tc; S[z].c = zc;
            for (uint32_t j = 0; j < tc; ++j){
                uint32_t v = S[t][j];
                _deg[t][v] = sdeg[t][v];
            }
            for (uint32_t j = 0; j < zc; ++j){
                uint32_t v = S[z][j];
                _deg[z][v] = sdeg[z][v];
            }
        
        qsz[0] = 0; qsz[1] = 0;
        for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
            uint32_t zzz = ttt^1;
            uint32_t v, dv;
            for (uint32_t j = 0 ; j < S[ttt].c; ++j) {
                v = S[ttt][j];
                dv = _deg[ttt][v];
                if (dv < sz[zzz] || dv * S[ttt].c <= mSize) {
                    vqueue[ttt][qsz[ttt]++] = v;
                    if ((v == u && ttt == t) || S[ttt].c < sz[ttt] || S[t].c * S[z].c <= mSize) {
                        qsz[ttt] = 0; qsz[zzz] = 0;
                        S[ttt].c = S[ttt].r = 0;
                        S[zzz].c = S[zzz].r = 0;
                        break;
                    }
                    S[ttt].swapByPos(--S[ttt].c, S[ttt].pos(v));
                }
            }
        }

        if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
            uint32_t qr[2] = {0};
            while ((qr[0] < qsz[0] || qr[1] < qsz[1]) && S[t].c >= sz[t]  && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
                for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                    uint32_t zzz = ttt^1;
                    for (uint32_t j = qr[ttt]; j < qsz[ttt]; ++j) {
                        uint32_t v = vqueue[ttt][j];
                        if (S[ttt].pos(v) < S[ttt].c) printf("ttt=%d, v=%d, S[ttt].pos(v)=%d, S[ttt].c=%d, j=%d, u=%d\n", ttt, v, S[ttt].pos(v), S[ttt].c, j, u);
                        assert(S[ttt].pos(v) >= S[ttt].c);

                        _deg[ttt][v] = 0;

                        for (uint32_t jj = g->p[ttt][v]; jj < g->p[ttt][v+1]; ++jj) {
                            uint32_t w = g->e[ttt][jj];
                            uint32_t wp = S[zzz].pos(w);
                            if (wp < S[zzz].c) {
                                _deg[zzz][w]--;
                                if (_deg[zzz][w] < sz[ttt] || _deg[zzz][w] * S[zzz].c <= mSize) {
                                    vqueue[zzz][qsz[zzz]++] = w;

                                    if ((w == u && zzz == t) || S[zzz].c < sz[zzz] || S[t].c * S[z].c <= mSize ) {
                                        qsz[ttt] = 0; qsz[zzz] = 0;
                                        S[ttt].c = S[ttt].r = 0;
                                        S[zzz].c = S[zzz].r = 0;
                                        break;
                                    }
                                    S[zzz].swapByPos(--S[zzz].c, wp);
                                }
                            }
                        }
                    }
                    qr[ttt] = qsz[ttt];
                }
            }
              

                // if (sz[tt] > sz[zz]) branchNew(0, 0);
                // else branch(0, 0);
                uint32_t dens = 0;
                for (uint32_t j = 0; j < S[t].c; ++j){
                    dens += _deg[t][S[t][j]]; 
                    _deg[t][S[t][j]] = 0;
                }
                for (uint32_t j = 0; j < S[z].c; ++j)
                    _deg[z][S[z][j]] = 0;

                assert(gt == tt);
                if (S[t].c / 100 >= S[z].c || S[z].c / 100 >= S[t].c) BIBRANCH = true;
                else BIBRANCH = false;

                // BIBRANCH = true;

                branchNew(0, 0);

                if (numbranches - old_numbranches > maxBranchduration || numbranches - old_numbranches >= 10000000) {
                    printf("it=%d, u=%d, sz[t]=%d, sz[z]=%d, ct=%d, cz=%d, maxBranchduration=%d, mSize=%d, density=%.3f\n", it, u, sz[t], sz[z], S[t].c, S[z].c, numbranches - old_numbranches, mSize, double(dens)/double(S[t].c*S[z].c));
                }
                maxBranchduration = std::max(maxBranchduration, numbranches - old_numbranches);

                // if (sz[tt] * sz[zz] < mSize && it < 2) {sz[zz] = std::max(sz[zz], mSize / sz[tt]);}
                // else if (sz[tt] * sz[zz] * 2 < mSize  && it >= 2) {sz[zz] = std::max(sz[zz], mSize / sz[tt] / 2);}

                // for (uint32_t j = 0; j < S[t].c; ++j)
                //     _deg[t][S[t][j]] = 0;
                // for (uint32_t j = 0; j < S[z].c; ++j)
                //     _deg[z][S[z][j]] = 0;
            }

            if (sz[zz] == crsz[zz] || sz[tt] == g->mxd[zz]) break;
            it++;
            // if (it > 20) { printf("error\n"); exit(1); }
        }

        for (uint32_t j = 0; j < tc; ++j){
            uint32_t v = S[t][j];
            sdeg[t][v] = 0;
        }
        for (uint32_t j = 0; j < zc; ++j){
            uint32_t v = S[z][j];
            sdeg[z][v] = 0;
        }

        if (u == 0) break;
        // if (u == 61205) {
        //     auto t2 = std::chrono::steady_clock::now();
        //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        //     std::cout << "iteration time:" << duration.count() << "ms" << std::endl;
        //     break;
        // }
    }

        // printf("it=%d, sz[0]=%d, sz[1]=%d, mSize=%d\n", it, sz[0], sz[1], mSize);
        // printf("it=%d, crsz[0]=%d, crsz[1]=%d, tt=%d\n", it, crsz[0], crsz[1], tt);
        // it++;
        // assert(it < 25);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        printf("mSize=%d, lsz=%d,rsz=%d, mxdeep=%d, numbranches=%ld, maxsubbranches=%ld, maxBranchduration=%d, subbcnt=%ld\n", mSize, realsz[0], realsz[1], mxdeep, numbranches, maxsubbranches, maxBranchduration, subbcnt);
        std::cout << "iteration time:" << duration.count() << "ms" << std::endl;
        // if (sz[zz] == crsz[zz] || sz[tt] == g->mxd[zz]) break;
    // }

    // printf("mSize=%d, lsz=%d,rsz=%d, mxdeep=%d, numbranches=%ld, maxsubbranches=%ld\n", mSize, realsz[0], realsz[1], mxdeep, numbranches, maxsubbranches);
    printRes();
}

#else 

void MBCS::run() {

    printf("run\n");
    // for(uint32_t i = 0; i < g->n[1]; ++i) {
    //     uint32_t d = g->p[1][i+1] - g->p[1][i];
    //     if (d >= 7591)
    //         printf("V v=%d, d=%d, is_nbr(12,53)=%d\n", i, d, g->is_nbr(0,12,53));
    // }
    S[0].c = g->n[0]; S[1].c = g->n[1];

    // mSize=15255;
    // mSize= 205180;
    // mSize = 880;

    initSize();

    uint32_t qsz[2] = {0};
    std::vector<uint32_t> vqueue[2];
    vqueue[0].resize(g->n[0]);
    vqueue[1].resize(g->n[1]);

    // return;

    // sz[0] = 7; sz[1] = 13;

    // printf("g->n[0]=%d, g->n[1]=%d, mSize=%d\n", g->n[0], g->n[1], mSize);
    // return;
    uint32_t t = 0, z = 1;
    if (g->maxDu < g->maxDv) {t = 1; z = 0;} 
    printf("g->n[0]=%d, g->n[1]=%d, mSize=%d, t=%d, z=%d\n", g->n[0], g->n[1], mSize, t, z);
    // printf("t=%d, u=21878, realu=%d\n", t, g->realLable[t][21878]); return;

    uint32_t crsz[2] = {0}, tt = 0, zz = 1;

    // tt = z; zz = tt^1;
    tt = t; zz = tt^1;

    crsz[0] = sz[0]; crsz[1] = sz[1];
    // sz[tt] = g->mxd[tt^1];
    sz[tt] = crsz[tt]; sz[zz] = g->mxd[tt] / 4;
    assert(crsz[0] > 0);
    assert(crsz[1] > 0);
    assert(sz[t] <= g->n[t]);

    if (sz[tt] * sz[zz] > mSize) { sz[zz] = std::max(crsz[zz], mSize / sz[tt] );}

    uint32_t it = 0;
    while (sz[tt] <= g->mxd[zz] && sz[zz] >= crsz[zz]) {

      
        if (it > 0 && it < 2) {
            sz[tt] = std::min(g->mxd[zz], sz[tt] +1);
            sz[zz] = std::max(crsz[zz], mSize / sz[tt]);
        }
        else if  (it == 2) {
            sz[tt] = std::min(g->mxd[zz], sz[tt] +1);
            sz[zz] = std::max(crsz[zz], mSize / sz[tt] / 2);
        }
        else if (it > 0) {
            sz[tt] = std::min(g->mxd[zz], sz[tt] * 2);
            sz[zz] = std::max(crsz[zz], mSize / sz[tt] / 2);
        }

        auto t1 = std::chrono::steady_clock::now();
    
    for (uint32_t u = std::max(g->n[t]-1, std::max((uint32_t)0,g->n[t]-sz[t])); u >= 0; u--) {

        S[t].r = S[t].c = 0;
        S[z].r = S[z].c = 0;
        assert(u < g->n[t]);
        assert(u >= 0);
        assert(S[t][S[t].pos(u)] == u);
        S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

        uint32_t Szs = g->deg(t,u);
        uint32_t tc = 0, zc = 0;

        if (Szs * g->mxd[z] <= mSize) {
            if (u == 0) break;
            continue;
        }

        for (uint32_t i = g->p[t][u]; i < g->p[t][u+1]; ++i) {
            uint32_t v = g->e[t][i];
            assert(S[z][S[z].pos(v)] == v);
            assert(v < g->n[z]);

            int nbrv = 0;
            for (uint32_t j = g->p[z][v+1]-1; j >= g->p[z][v]; j--) {
                uint32_t w = g->e[z][j];
                if (j == 0 || w < u) break;
                nbrv ++;
            }
            if (Szs * nbrv <= mSize) {
                Szs -- ;
                continue;
            }

            S[z].swapByPos(S[z].c++, S[z].pos(v));
            sdeg[z][v] = nbrv;

            for (uint32_t j = g->p[z][v+1]-1; j >= g->p[z][v]; j--) {
                uint32_t w = g->e[z][j];
                uint32_t wp = S[t].pos(w);
                // assert(wp >= 0);
                // assert(wp < g->n[t]);
                // assert(S[t][wp] == w);
                if (j == 0 || w < u) break;
                if (wp >= S[t].c)
                    S[t].swapByPos(S[t].c++, wp);
                sdeg[t][w]++;
            }
        }

        tc = S[t].c; zc = S[z].c;

        // reduction
        
        qsz[0] = 0; qsz[1] = 0;
        for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
            uint32_t zzz = ttt^1;
            uint32_t v, dv;
            for (uint32_t j = 0 ; j < S[ttt].c; ++j) {
                v = S[ttt][j];
                dv = sdeg[ttt][v];
                if (dv < sz[zzz] || dv * S[ttt].c <= mSize) {
                    vqueue[ttt][qsz[ttt]++] = v;
                    S[ttt].swapByPos(--S[ttt].c, S[ttt].pos(v));
                    if ((v == u && ttt == t) || S[ttt].c < sz[ttt] || S[t].c * S[z].c <= mSize) {
                        qsz[ttt] = 0; qsz[zzz] = 0;
                        S[ttt].c = S[ttt].r = 0;
                        S[zzz].c = S[zzz].r = 0;
                        break;
                    }
                }
            }
        }

        if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
            uint32_t qr[2] = {0};
            while ((qr[0] < qsz[0] || qr[1] < qsz[1]) && S[t].c >= sz[t]  && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
                for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                    uint32_t zzz = ttt^1;
                    for (uint32_t j = qr[ttt]; j < qsz[ttt]; ++j) {
                        uint32_t v = vqueue[ttt][j];
                        if (S[ttt].pos(v) < S[ttt].c) printf("ttt=%d, v=%d, S[ttt].pos(v)=%d, S[ttt].c=%d, j=%d, u=%d\n", ttt, v, S[ttt].pos(v), S[ttt].c, j, u);
                        assert(S[ttt].pos(v) >= S[ttt].c);

                        sdeg[ttt][v] = 0;

                        for (uint32_t jj = g->p[ttt][v]; jj < g->p[ttt][v+1]; ++jj) {
                            uint32_t w = g->e[ttt][jj];
                            uint32_t wp = S[zzz].pos(w);
                            if (wp < S[zzz].c) {
                                sdeg[zzz][w]--;
                                if (sdeg[zzz][w] < sz[ttt] || sdeg[zzz][w] * S[zzz].c <= mSize) {
                                    vqueue[zzz][qsz[zzz]++] = w;
                                    S[zzz].swapByPos(--S[zzz].c, wp);

                                    if ((w == u && zzz == t) || S[zzz].c < sz[zzz] || S[t].c * S[z].c <= mSize ) {
                                        qsz[ttt] = 0; qsz[zzz] = 0;
                                        S[ttt].c = S[ttt].r = 0;
                                        S[zzz].c = S[zzz].r = 0;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    qr[ttt] = qsz[ttt];
                }
            }

            for (uint32_t j = 0; j < tc; ++j)
                sdeg[t][S[t][j]] = 0;
            for (uint32_t j = 0; j < zc; ++j)
                sdeg[z][S[z][j]] = 0;

            // if (sz[tt] > sz[zz]) branchNew(0, 0);
            // else branch(0, 0);

            branchNew(0, 0);

            if (sz[tt] * sz[zz] < mSize && it < 2) {sz[zz] = std::max(sz[zz], mSize / sz[tt]);}
            else if (sz[tt] * sz[zz] * 2 < mSize && it >= 2) {sz[zz] = std::max(sz[zz], mSize / sz[tt] / 2);}

            // for (uint32_t j = 0; j < tc; ++j)
            //     sdeg[t][S[t][j]] = 0;
            // for (uint32_t j = 0; j < zc; ++j)
            //     sdeg[z][S[z][j]] = 0;
        }

        if (u == 0) break;
    }

        printf("it=%d, sz[0]=%d, sz[1]=%d, mSize=%d\n", it, sz[0], sz[1], mSize);
        // printf("it=%d, crsz[0]=%d, crsz[1]=%d, tt=%d\n", it, crsz[0], crsz[1], tt);
        it++;
        // assert(it < 25);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        printf("mSize=%d, lsz=%d,rsz=%d, mxdeep=%d, numbranches=%ld, maxsubbranches=%ld\n", mSize, realsz[0], realsz[1], mxdeep, numbranches, maxsubbranches);
        std::cout << "iteration time:" << duration.count() << "ms" << std::endl;
        if (sz[zz] == crsz[zz] || sz[tt] == g->mxd[zz]) break;
    }

    printf("maxDu=%d, maxDv=%d\n", g->maxDu, g->maxDv);
    // printf("mSize=%d, lsz=%d,rsz=%d, mxdeep=%d, numbranches=%ld, maxsubbranches=%ld\n", mSize, realsz[0], realsz[1], mxdeep, numbranches, maxsubbranches);
    printRes();
}
#endif

void MBCS::branch1(uint32_t deep, uint32_t ne) {
    if (S[0].c < sz[0] || S[1].c < sz[1]) return;
    if (mSize >= S[0].c * S[1].c) return;
    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (S[0].cSize() == 0) {
        if (mSize < S[0].r * S[1].c) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // realRes[0].clear(); realRes[1].clear();
        }
        return;
    }
    else if (S[1].cSize() == 0) {
        if (mSize < S[0].c * S[1].r) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
        }
    
        return;
    }

    uint32_t r[2], c[2];
    r[0] = S[0].r; r[1] = S[1].r;
    c[0] = S[0].c; c[1] = S[1].c;

    for (uint32_t t = 0; t <= 1; ++t) 
    if (S[t].cSize() == 1) {
        uint32_t d = 0;
        uint32_t z = t^1;
        uint32_t u = S[t][S[t].r];
        for (uint32_t i = S[z].r; i < S[z].c; ++i) {
            uint32_t v = S[z][i];
            if (g->is_nbr(t, u, v)) d++;
        }
        if ((S[z].r+d) * (S[t].r+1) > mSize  &&  S[t].r+1 >= sz[t] && S[z].r+d >= sz[z]) {
            // mSize = (S[z].r+d) * (S[t].r+1);
            mSize = (S[z].r+d) * (S[t].r+1); realsz[t] = (S[t].r+1); realsz[z] = (S[z].r+d);

            // gemSize(t,S[t].r+1, S[z].r+d, ne);
        }
        else if (S[z].c * S[t].r > mSize &&  S[t].r >= sz[t] && S[z].c >= sz[z]) {
            // mSize = S[z].c * S[t].r;
            mSize = S[z].c * S[t].r; realsz[t] = S[t].r; realsz[z] = S[z].c;
            // assert(S[z].c * S[t].r != 880);
            // gemSize(t, S[t].r, S[z].c, ne);
        }
        return;
    }

    for (uint32_t t = 0; t <= 1; ++t) {
        if (S[t].c == sz[t]) {
            uint32_t z = t^1;
            for (int i = S[z].r; i < S[z].c; ++i) {
                uint32_t u = S[z][i];
                for (uint32_t j = S[t].r; j < S[t].c; ++j) {
                    uint32_t v = S[t][j];
                    if (!g->is_nbr(z,u,v)) {
                        S[z].swapByPos(i--, --S[z].c);
                        break;
                    }

                }
                if (S[t].c * S[z].c <= mSize) break;
            }

            S[t].r = S[t].c;
            branch(deep+1, ne);

            S[t].r = r[t];
            S[t].c = c[t];
            S[z].c = c[z];
            return;
        }
    }

    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};

    uint32_t t = 0, z = 1;
    uint32_t sube = 0;
    if (S[1].cSize() < S[0].cSize()) {t = 1; z = 0;} 

#ifdef UB1_
    bool fg = true;
    while (fg) {
        fg = false; 

        for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
        // for(uint32_t i = S[t].r; i < S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = sdeg[t][u];

            // if(g->deg(t, u) > S[z].cSize()) {
            //     for(uint32_t j = S[z].r; j < S[z].c; j++) {
            //         uint32_t v = S[z][j];
            //         if(g->is_nbr(t, u, v)) {
            //             e++; sdeg[z][v]++;
            //         }
            //     }
            // }
            // else {
            //     for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            //         uint32_t v = g->e[t][j];
            //         uint32_t pv = S[z].pos(v);
            //         if(S[z].r <= pv && pv < S[z].c) {//in C
            //             e++; sdeg[z][v]++;
            //         }
            //     }
            // }
    
            if (e < sz[z] || S[t].c * e < mSize) {
                S[t].swapByPos(--S[t].c, i--);

                if(g->deg(t, u) > c[z]) {
                    for(uint32_t j = 0; j < c[z]; j++) {
                        uint32_t v = S[z][j];
                        if(g->is_nbr(t, u, v))
                            sdeg[z][v]--;
                    }
                }
                else {
                    for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                        uint32_t v = g->e[t][j];
                        uint32_t pv = S[z].pos(v);
                        if(pv < c[z]) {//in C and R
                            sdeg[z][v]--;
                        }
                    }
                }

                if (S[t].c < sz[t] || S[0].c * S[1].c <= mSize)  {
                    // for(uint32_t j = S[z].r; j < S[z].c; j++) 
                    //     sdeg[z][S[z][j]] = 0;

                    reSetdeg(c);
                    
                    S[t].c = c[t];
                    S[z].c = c[z];
                    return;
                }

            }
 

        }

        for(int i = S[z].r; i < S[z].c; i++) {
            uint32_t u = S[z][i];
            uint32_t e = sdeg[z][u];

            if (e < sz[t] || S[z].c * e < mSize) {
                S[z].swapByPos(--S[z].c, i--);
                fg = true;

                if(g->deg(z, u) > c[t]) {
                    for(uint32_t j = 0; j < c[t]; j++) {
                        uint32_t v = S[t][j];
                        if(g->is_nbr(z, u, v))
                            sdeg[t][v]--;
                    }
                }
                else {
                    for(uint32_t j = g->p[z][u]; j < g->p[z][u + 1]; j++) {
                        uint32_t v = g->e[z][j];
                        uint32_t pv = S[t].pos(v);
                        if(pv < c[t]) {//in C and R
                            sdeg[t][v]--;
                        }
                    }
                }


            }
        }

        if (S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  {

            reSetdeg(c);

            S[t].c = c[t];
            S[z].c = c[z];
            return;
        }
    }
#endif
    // sube = 0;

    maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].r; i < S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = sdeg[t][u];

        if(e > maxE[t]) {
            maxE[t] = e;
            maxI[t] = i;
            maxUV[t] = u;
            // if (e >= S[z].cSize()) break;
        }
        // sdeg[t][u] = e;
        // sube += e;
    }
    for(int i = S[z].r; i < S[z].c; i++) {
        uint32_t u = S[z][i];
        uint32_t e = sdeg[z][u];

        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[z] = i;
            maxUV[z] = u;
        }
    }
    // if (S[z].r > r[z]) {

    //     branch(deep+1, ne);
    //     S[z].r = r[z];
    //     S[t].c = c[t];
    //     S[z].c = c[z];
    //     return;
    // }

    // if(S[0].r * S[1].c + S[1].r * S[0].c - S[0].r*S[1].r +sube  <= mSize) {
    //     S[t].c = c[t];
    //     S[z].c = c[z]; 
    //     return;
    // }

    if (S[0].cSize() == 0) {

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }

        reSetdeg(c);

        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    if (S[1].cSize() == 0) {

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
        
        reSetdeg(c);

        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    
    // printf("maxE[0]=%d, maxE[1]=%d\n", maxE[0],maxE[1]);
    // exit(1);
    if (maxE[t] == 0) {

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }

        reSetdeg(c);
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }

    assert(S[0].cSize() > 0);
    assert(S[1].cSize() > 0);


    // t = 0, z = 1;
    if (S[z].cSize() - maxE[t] > S[t].cSize()-maxE[z]) {
        t = z; z = t^1;
    }
    // printf("t=%d, z=%d\n", t,z);
    // exit(1);

    uint32_t u = maxUV[t], wsSize = 0;

    if(g->deg(t, u) > S[z].cSize()) {
        if (S[z].cSize() > ws[deep].capacity()) ws[deep].resize(S[z].cSize()+1);
        for(uint32_t i = S[z].r; i < S[z].c; i++) {
            if(!g->is_nbr(t, u, S[z][i])) {
                ws[deep][wsSize++] = S[z][i];
            }
        }
    }
    else {
        uint32_t cc = S[z].c;
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].r <= pv && pv < S[z].c) {
                S[z].swapByPos(--cc, pv);
            }
        }
        wsSize = cc - S[z].r;
        if (wsSize > ws[deep].capacity()) ws[deep].resize(wsSize+1);
        memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);
    }

    maxsubbranches = wsSize+1 > maxsubbranches ? wsSize+1 : maxsubbranches;
    fflush(stdout); 
    uint32_t ct = S[t].c, cz = S[z].c;
    for (uint32_t i = 0; i < wsSize; ++i) {

        uint32_t v = ws[deep][i];
      
        S[z].swapByPos(S[z].r++, S[z].pos(v));
        if(g->deg(z, v) > S[t].cSize()) {
            for(int j = S[t].r; j < S[t].c; j++) {
                if(!g->is_nbr(z, v, S[t][j])) {
                    S[t].swapByPos(--S[t].c, j); j--;
                }
            }
        }
        else {
            uint32_t cc = S[t].r;
            for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                uint32_t w = g->e[z][j];
                uint32_t pw = S[t].pos(w);
                if(S[t].r <= pw && pw < S[t].c) {
                    S[t].swapByPos(cc++, pw);
                }
            }
            S[t].c = cc;
        }


        for (uint32_t ii = S[t].c; ii < ct; ++ii) {
            uint32_t  w = S[t][ii];
            if (g->deg(t,w) > cz) {
                for(int j = 0; j < cz; j++) {
                    uint32_t w1 = S[z][j];
                    if(g->is_nbr(t, w, w1)) {
                        sdeg[z][w1]--;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][w]; j < g->p[t][w + 1]; j++) {
                    uint32_t w1 = g->e[t][j];
                    uint32_t pw1 = S[z].pos(w1);
                    if(pw1 < cz) {
                        sdeg[z][w1]--;
                    }
                }
            }

        }
        
        branch(deep+1, ne);

        for (uint32_t ii = S[t].c; ii < ct; ++ii) {
            uint32_t  w = S[t][ii];
            if (g->deg(t,w) > cz) {
                for(int j = 0; j < cz; j++) {
                    uint32_t w1 = S[z][j];
                    if(g->is_nbr(t, w, w1)) {
                        sdeg[z][w1]++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][w]; j < g->p[t][w + 1]; j++) {
                    uint32_t w1 = g->e[t][j];
                    uint32_t pw1 = S[z].pos(w1);
                    if(pw1 < cz) {
                        sdeg[z][w1]++;
                    }
                }
            }

        }


        S[t].c = ct;
        S[z].swapByPos(--S[z].r,--S[z].c);


        if(g->deg(z, v) > ct) {
            for(int j = 0; j < ct; j++) {
                uint32_t w = S[t][j];
                if(g->is_nbr(z, v, w)) {
                    sdeg[t][w]--;
                }
            }
        }
        else {
            for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                uint32_t w = g->e[z][j];
                uint32_t pw = S[t].pos(w);
                if(pw < ct) 
                    sdeg[t][w]--;
            }
        }
    }

    
    // pivot vertex
    S[t].swapByPos(S[t].r++, S[t].pos(u));

    branch(deep+1, ne);

    reSetdeg(c);

    S[t].r = r[t];
    S[t].c = c[t];
    S[z].r = r[z];
    S[z].c = c[z];

}


void MBCS::branch(uint32_t deep, uint32_t ne) {
    if (S[0].c < sz[0] || S[1].c < sz[1]) return;
// if (S[0][0] == 12)
// printf("deep=%d, S[0],r=%d, S[0].c=%d, S[1],r=%d, S[1].c=%d, S[0].size=%d, S[1].size=%d, mSize=%d\n", deep, S[0].r, S[0].c, S[1].r, S[1].c, S[0].cSize(), S[1].cSize(), mSize);
    if (mSize >= S[0].c * S[1].c) return;

    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (S[0].cSize() == 0) {
        if (mSize < S[0].r * S[1].c) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
        return;
    }
    else if (S[1].cSize() == 0) {
        if (mSize < S[0].c * S[1].r) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
    
        return;
    }

    uint32_t r[2], c[2];
    r[0] = S[0].r; r[1] = S[1].r;
    c[0] = S[0].c; c[1] = S[1].c;

    for (uint32_t t = 0; t <= 1; ++t) 
    if (S[t].cSize() == 1) {
        uint32_t d = 0;
        uint32_t z = t^1;
        uint32_t u = S[t][S[t].r];
        for (uint32_t i = S[z].r; i < S[z].c; ++i) {
            uint32_t v = S[z][i];
            if (g->is_nbr(t, u, v)) d++;
        }
        if ((S[z].r+d) * (S[t].r+1) > mSize  &&  S[t].r+1 >= sz[t] && S[z].r+d >= sz[z]) {
            // mSize = (S[z].r+d) * (S[t].r+1);
            mSize = (S[z].r+d) * (S[t].r+1); realsz[t] = (S[t].r+1); realsz[z] = (S[z].r+d);
            // assert((S[z].r+d) * (S[t].r+1) != 880);
            
            // gemSize(t,S[t].r+1, S[z].r+d, ne);
        }
        else if (S[z].c * S[t].r > mSize &&  S[t].r >= sz[t] && S[z].c >= sz[z]) {
            // mSize = S[z].c * S[t].r;
            mSize = S[z].c * S[t].r; realsz[t] = S[t].r; realsz[z] = S[z].c;
            // assert(S[z].c * S[t].r != 880);
            // gemSize(t, S[t].r, S[z].c, ne);
        }
        return;
    }

    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};

    uint32_t t = 0, z = 1;
    uint32_t sube = 0;
    if (S[1].cSize() < S[0].cSize()) {t = 1; z = 0;} 

    if (S[t].c == sz[t]) {

        for (int i = S[z].r; i < S[z].c; ++i) {
            uint32_t u = S[z][i];
            for (uint32_t j = S[t].r; j < S[t].c; ++j) {
                uint32_t v = S[t][j];
                if (!g->is_nbr(z,u,v)) {
                    S[z].swapByPos(i--, --S[z].c);
                    break;
                }

            }
            if (S[t].c * S[z].c <= mSize) break;
        }

        S[t].r = S[t].c;
        branch(deep+1, ne);
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    if (S[z].c == sz[z]) {

        for (int i = S[t].r; i < S[t].c; ++i) {
            uint32_t u = S[t][i];
            for (uint32_t j = S[z].r; j < S[z].c; ++j) {
                uint32_t v = S[z][j];
                if (!g->is_nbr(t,u,v)) {
                    S[t].swapByPos(i--, --S[t].c);
                    break;
                }

            }
            if (S[t].c * S[z].c <= mSize) break;
        }

        S[z].r = S[z].c;
        branch(deep+1, ne);
        S[z].r = r[z];
        S[z].c = c[z];
        S[t].c = c[t];
        return;
    }


#ifdef UB1_
    bool fg = true;
    while (fg) {
        fg = false; 

        for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;
            if(g->deg(t, u) > S[z].cSize()) {
                for(uint32_t j = S[z].r; j < S[z].c; j++) {
                    uint32_t v = S[z][j];
                    if(g->is_nbr(t, u, v)) {
                        e++; sdeg[z][v]++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].r <= pv && pv < S[z].c) {//in C
                        e++; sdeg[z][v]++;
                    }
                }
            }
    
            if (S[z].r+e < sz[z] || S[t].c *(S[z].r+e) < mSize) {
                S[t].swapByPos(--S[t].c, i--);

                if (S[t].c < sz[t])  {
                    for(uint32_t j = S[z].r; j < S[z].c; j++) 
                        sdeg[z][S[z][j]] = 0;
                    
                    S[t].c = c[t];
                    S[z].c = c[z];
                    return;
                }

                // if(g->deg(t, u) > S[z].cSize()) {
                //     for(uint32_t j = S[z].r; j < S[z].c; j++) {
                //         uint32_t v = S[z][j];
                //         if(g->is_nbr(t, u, v)) {
                //             sdeg[z][v]--;
                //         }
                //     }
                // }
                // else {
                //     for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                //         uint32_t v = g->e[t][j];
                //         uint32_t pv = S[z].pos(v);
                //         if(S[z].r <= pv && pv < S[z].c) {//in C
                //             sdeg[z][v]--;
                //         }
                //     }
                // }

            }
            // else if (e == S[z].cSize()) {
            //     S[t].swapByPos(S[t].r++, i);
            //     if(g->deg(t, u) > S[z].cSize()) {
            //     for(uint32_t j = S[z].r; j < S[z].c; j++) {
            //         uint32_t v = S[z][j];
            //         if(g->is_nbr(t, u, v)) {
            //                 sdeg[z][v]--;
            //             }
            //         }
            //     }
            //     else {
            //         for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            //             uint32_t v = g->e[t][j];
            //             uint32_t pv = S[z].pos(v);
            //             if(S[z].r <= pv && pv < S[z].c) {//in C
            //                 sdeg[z][v]--;
            //             }
            //         }
            //     }

            // }
            // else sube -= ();
            // if(e > maxE[t]) {
            //     maxE[t] = e;
            //     maxI[t] = i;
            //     maxUV[t] = u;
            //     // if (e >= S[z].cSize()) break;
            // }

        }

        // if (S[t].r > r[t]) {
        //     for(uint32_t j = S[z].r; j < S[z].c; j++) 
        //         sdeg[z][S[z][j]] = 0;
        //     branch(deep+1, ne);
        //     S[t].r = r[t];
        //     S[t].c = c[t];
        //     S[z].c = c[z];
        //     return;
        // }

    // }

    // printf("t=%d, z=%d, maxE[0]=%d, maxE[1]=%d\n", t,z,maxE[0],maxE[1]);
        // for(int i = S[z].r; fg1 == false && i < S[z].c; i++) {
        //     // uint32_t u = S[z][i];
        //     uint32_t e = sdeg[z][S[z][i]];
        //     if (e == S[t].cSize()) {
        //         S[z].swapByPos(S[z].r++, i);
        //         maxE[t] --;
        //     }
        // }

        for(int i = S[z].r; i < S[z].c; i++) {
            uint32_t u = S[z][i];
            uint32_t e = sdeg[z][u];
            sdeg[z][u] = 0;

            if (S[t].r+e < sz[t] || S[z].c *(S[t].r+e) < mSize) {
                sdeg[z][u] = 0;
                S[z].swapByPos(--S[z].c, i--);
                fg = true;

                // if(g->deg(z, u) > S[t].cSize()) {
                //     for(uint32_t j = S[t].r; j < S[t].c; j++) {
                //         uint32_t v = S[t][j];
                //         if(g->is_nbr(z, u, v)) {
                //             sdeg[t][v]--;
                //         }
                //     }
                // }
                // else {
                //     for(uint32_t j = g->p[z][u]; j < g->p[z][u + 1]; j++) {
                //         uint32_t v = g->e[z][j];
                //         uint32_t pv = S[t].pos(v);
                //         if(S[t].r <= pv && pv < S[t].c) {//in C
                //             sdeg[t][v]--;
                //         }
                //     }
                // }

            }
            // else if(e > maxE[z]) {
            //     maxE[z] = e;
            //     maxI[z] = i;
            //     maxUV[z] = u;
            // }
        }
        if (S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  {
            S[t].c = c[t];
            S[z].c = c[z];
            return;
        }
    }
#endif
    // sube = 0;

    maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;
        if(g->deg(t, u) > S[z].cSize()) {
            for(uint32_t j = S[z].r; j < S[z].c; j++) {
                uint32_t v = S[z][j];
                if(g->is_nbr(t, u, v)) {
                    e++; sdeg[z][v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].r <= pv && pv < S[z].c) {//in C
                    e++; sdeg[z][v]++;
                }
            }
        }

        // if (e == S[z].cSize()) {
        //     S[t].swapByPos(S[t].r++, i);
        //     continue;
        // }

        if(e > maxE[t]) {
            maxE[t] = e;
            maxI[t] = i;
            maxUV[t] = u;
            // if (e >= S[z].cSize()) break;
        }
        // sdeg[t][u] = e;
        // sube += e;
    }
    for(int i = S[z].r; i < S[z].c; i++) {
        uint32_t u = S[z][i];
        uint32_t e = sdeg[z][u];
        sdeg[z][u] = 0;
        // if (e == S[t].c - r[t]) {
        //     // uint32_t d = 0;
        //     // for (uint32_t j = S[t].r; j  < S[t].c; ++j)
        //     //     if (g->is_nbr(z,u, S[t][j])) d++;
        //     // assert(d == e);
        //     // assert(i == S[z].pos(u));
        //     S[z].swapByPos(S[z].r++, i);
        // }
        // else 
        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[z] = i;
            maxUV[z] = u;
        }
    }

    // if (mxdeep == 2000) {
    //     S[0].print();
    //     S[1].print();

    //     printf("maxE[0]=%d, maxE[1]=%d\nmaxI[0]=%d, maxI[1]=%d\nmaxUV[0]=%d, maxUV[1]=%d\n", maxE[0], maxE[1], maxI[0], maxI[1], maxUV[0], maxUV[1]);
    //     printf("S[0].c=%d, S[1].c=%d\n", S[0].c, S[1].c);
    // }
    // if (S[z].r > r[z]) {

    //     branch(deep+1, ne);
    //     S[z].r = r[z];
    //     S[t].c = c[t];
    //     S[z].c = c[z];
    //     return;
    // }

    // if(S[0].r * S[1].c + S[1].r * S[0].c - S[0].r*S[1].r +sube  <= mSize) {
    //     S[t].c = c[t];
    //     S[z].c = c[z]; 
    //     return;
    // }

    if (S[0].cSize() == 0) {

        for(uint32_t i = S[t].r; i < S[t].c; i++) 
            sdeg[t][S[t][i]] = 0;
        for(uint32_t i = S[z].r; i < S[z].c; i++) 
            sdeg[z][S[z][i]] = 0;

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    if (S[1].cSize() == 0) {
        for(uint32_t i = S[t].r; i < S[t].c; i++) 
            sdeg[t][S[t][i]] = 0;
        for(uint32_t i = S[z].r; i < S[z].c; i++) 
            sdeg[z][S[z][i]] = 0;

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
  
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    
    // printf("maxE[0]=%d, maxE[1]=%d\n", maxE[0],maxE[1]);
    // exit(1);
    if (maxE[t] == 0) {
        for(uint32_t i = S[t].r; i < S[t].c; i++) 
            sdeg[t][S[t][i]] = 0;
        for(uint32_t i = S[z].r; i < S[z].c; i++) 
            sdeg[z][S[z][i]] = 0;
        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }

    assert(S[0].cSize() > 0);
    assert(S[1].cSize() > 0);

    // if(upperCan(deep) < mSize) {
    //     S[t].c = c[t];
    //     S[z].c = c[z]; 
    //     return;
    // }

    // if (maxE[t]+1 >= S[z].cSize()) {
    //     if (maxE[t] == S[z].cSize()) {
    //         S[t].swapByPos(S[t].r++, S[t].pos(maxUV[t]));
    //         for(uint32_t i = S[z].r; i < S[z].c; i++) 
    //             subdeg[S[z][i]] = 0;
    //         branch(deep+1, ne);
    //         S[t].r--;
    //     }
    //     else {
            
    //         for(uint32_t i = S[z].r; i < S[z].c; i++) 
    //             subdeg[S[z][i]] = 0;
    //     }
    //     return;
    // }
    
    // for(uint32_t i = S[z].r; i < S[z].c; i++) {
    //     uint32_t u = S[z][i];
    //     uint32_t e = subdeg[u];
    //     subdeg[u] = 0;

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[z] = i;
    //         maxUV[z] = u;
    //     }
    // }

    // t = 0, z = 1;
    if (S[z].cSize() - maxE[t] > S[t].cSize()-maxE[z]) {
        t = z; z = t^1;
    }
    // printf("t=%d, z=%d\n", t,z);
    // exit(1);

    uint32_t u = maxUV[t], wsSize = 0;

    if(g->deg(t, u) > S[z].cSize()) {
        if (S[z].cSize() > ws[deep].capacity()) ws[deep].resize(S[z].cSize()+1);
        for(uint32_t i = S[z].r; i < S[z].c; i++) {
            if(!g->is_nbr(t, u, S[z][i])) {
                ws[deep][wsSize++] = S[z][i];
            }
        }
    }
    else {
        uint32_t cc = S[z].c;
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].r <= pv && pv < S[z].c) {
                S[z].swapByPos(--cc, pv);
            }
        }
        wsSize = cc - S[z].r;
        if (wsSize > ws[deep].capacity()) ws[deep].resize(wsSize+1);
        memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);
    }

    // if (S[0][0] == 12){
    //     printf("deep=%d, S[0],r=%d, S[0].c=%d, S[1],r=%d, S[1].c=%d, wsSize=%d, t=%d, z=%d, S[0].size=%d, S[1].size=%d\n", deep, S[t].r, S[t].c, S[z].r, S[z].c, wsSize, t, z, S[0].cSize(), S[1].cSize());
    //     printf("maxE[0]=%d, maxI[0]=%d, maxUV[0]=%d\n", maxE[0], maxI[0], maxUV[0]);
    //     printf("maxE[1]=%d, maxI[1]=%d, maxUV[1]=%d\n", maxE[1], maxI[1], maxUV[1]);
    //     printf("ws:");
    //     for (auto i = 0; i < wsSize; ++i) 
    //         printf("%d ", ws[deep][i]);
    //     printf("\n");
    // }
    maxsubbranches = wsSize+1 > maxsubbranches ? wsSize+1 : maxsubbranches;
    fflush(stdout); 
    uint32_t ct = S[t].c;
    for (uint32_t i = 0; i < wsSize; ++i) {

        uint32_t v = ws[deep][i];
        // if (deep == 0) {
        //     printf("z=%d, v=%d, i=%d\n", z, v, i);fflush(stdout); 
        //     printf("ws:");
        //     for ( auto i = 0; i < wsSize; ++i) 
        //         printf("%d ", ws[deep][i]);
        //     printf("\n");

        //     if (i >=12 ) exit(1);
        // }
        S[z].swapByPos(S[z].r++, S[z].pos(v));
        // if (g->deg(z, v) < 0) printf("z=%d, v=%d, g->n[z]=%d, g->deg(z, v)=%d\n", z, v, g->n[z], g->deg(z, v)); fflush(stdout);  
        if(g->deg(z, v) > S[t].cSize()) {
            for(int j = S[t].r; j < S[t].c; j++) {
                // assert(v >= 0);
                // if (v >= g->n[z]) {
                //     printf("z=%d, v=%d, g->n[z]=%d, g->deg(z, v)=%d\n", z, v, g->n[z], g->deg(z, v)); fflush(stdout);  
                //     printf("deep=%d, r=%d, wsSize=%d, c=%d, t=%d, z=%d, i=%d\n", deep, S[z].r, wsSize, S[z].c, t, z, i);

                //     printf("ws:");
                //     for ( auto i = 0; i < wsSize; ++i) 
                //         printf("%d ", ws[deep][i]);
                //     printf("\n");
                //     fflush(stdout);  
                // }
                // assert(v < g->n[z]);
                if(!g->is_nbr(z, v, S[t][j])) {
                    S[t].swapByPos(--S[t].c, j); j--;
                }
            }
        }
        else {
            uint32_t cc = S[t].r;
            for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                uint32_t w = g->e[z][j];
                uint32_t pw = S[t].pos(w);
                if(S[t].r <= pw && pw < S[t].c) {
                    S[t].swapByPos(cc++, pw);
                }
            }
            S[t].c = cc;
        }
        // if (S[0][0] == 12) {
        //     printf("z=%d, v=%d, i=%d\n", z, v, i);fflush(stdout); 
        //     printf("deep=%d, S[0],r=%d, S[0].c=%d, S[1],r=%d, S[1].c=%d, wsSize=%d, t=%d, z=%d, S[0].size=%d, S[1].size=%d\n", deep, S[t].r, S[t].c, S[z].r, S[z].c, wsSize, t, z, S[0].cSize(), S[1].cSize());
        // }
        
        branch(deep+1, ne);
        S[t].c = ct;
        S[z].swapByPos(--S[z].r,--S[z].c);
    
        // if (S[0][0] == 12) {
        //     exit(-1);
        // }
        // S[t].c = c[t];
        // S[z].swapByPos(--S[z].r,--S[z].c);

        // if (deep == 2) {
        //     printf("z=%d, v=%d, done\n", z, v);fflush(stdout); 
        // }
    }

    
    // pivot vertex
    S[t].swapByPos(S[t].r++, S[t].pos(u));

    branch(deep+1, ne);

    S[t].r = r[t];
    S[t].c = c[t];
    S[z].r = r[z];
    S[z].c = c[z];

}

void MBCS:: initSize() {
auto t1 = std::chrono::steady_clock::now();
    uint32_t t = 0, z = 1;
    if (g->maxDu < g->maxDv) {t = 1; z = 0;} 
    t = 0; z = 1;
    // printf("t=%d,z=%d, sz[t]=%d, sz[z]=%d\n", t, z, sz[t], sz[z]);

    for (uint32_t t = 0; t <= 1; ++t) {
    z = t^1;
    std::vector<bool> visited(g->n[t], false);
    uint32_t its = 0;
    for (int u = std::max(g->n[t]-sz[t], uint32_t(0)); u >= 0; u--) {
        if (g->cores[t][u] < g->maxcore) break;
        if (!visited[u]) {

            S[t].r = S[t].c = 0;
            S[z].r = S[z].c = 0;
            S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

            uint32_t nextv = g->n[z]+1, mdv = 0, mdc = 0, du = g->p[t][u+1]- g->p[t][u];
            uint32_t cnts = 0;

            //test
            int oldv = -1;
            for (uint32_t j = g->p[t][u]; j < g->p[t][u+1]; ++j) {
                uint32_t v = g->e[t][j];
                uint32_t dv = g->p[z][v+1] - g->p[z][v];
                if (oldv >= 0) assert(oldv < v); oldv = v;
                if (dv >= sz[t] && mSize < dv * du) {
                    // if (mdc < g->cores[z][v]) {
                    //     nextv = v; mdc = g->cores[z][v]; mdv = dv;
                    // }
                    // else if (dv > mdv) {
                    //     nextv = v; mdv = dv;
                    // }
                    // cnts++;

                    if (dv >= mdv) {
                        nextv = v; mdv = dv;
                    }
                    S[z].swapByPos(S[z].c++, S[z].pos(v));
                }
            }

            if (nextv >= g->n[z]) assert(mdv == 0);
            if (nextv >= g->n[z]) continue;
            
            S[z].swapByPos(S[z].r++, S[z].pos(nextv));


            cnts = 0;
            for (uint32_t j = g->p[z][nextv]; j < g->p[z][nextv+1]; ++j){
                uint32_t u1 = g->e[z][j];
                if (u != u1) {
                    S[t].swapByPos(S[t].c++, S[t].pos(u1));
                    cnts++;
                }
            }

            cnts = 0;
            while (S[z].cSize() > 0 && S[t].cSize() > 0 && S[t].c *S[z].c > mSize) {
                nextv = g->n[z]+1; mdv = 0; mdc = 0;
            
                for (int i = S[z].r; i < S[z].c; ++i) {
                    uint32_t v = S[z][i];
                    uint32_t dv = g->p[z][v+1] - g->p[z][v];

                    if (dv >= sz[t] && mSize < dv * S[z].c) {
                        // if (mdc < g->cores[z][v]) {
                        //     nextv = v; mdc = g->cores[z][v]; mdv = dv;
                        // }
                        // else 
                        if (dv >= mdv) {
                            nextv = v; mdv = dv;
                        }
                    }
                    if (mSize >= dv * S[z].c) {
                        S[z].swapByPos(--S[z].c, i); i--;
                    }

                }
                if (nextv >= g->n[z]) assert(S[z].cSize() == 0);
                if (nextv >= g->n[z]) break;

                assert(nextv < g->n[z]);
                S[z].swapByPos(S[z].r++, S[z].pos(nextv));

                if (g->deg(z,nextv) <= S[t].cSize()) {
                    uint32_t stsz = S[t].r;
                    for (uint32_t j = g->p[z][nextv]; j < g->p[z][nextv+1]; ++j){
                        uint32_t u1 = g->e[z][j];
                        uint32_t pu = S[t].pos(u1);
                        if (pu >= stsz && pu < S[t].c) {
                            S[t].swapByPos(stsz++, pu);
                        }
                    }
                    S[t].c = stsz;
                }
                else {
                    for (uint32_t j = S[t].r; j < S[t].c; ++j) {
                        uint32_t u1 = S[t][j];
                        if (!g->is_nbr(z, nextv, u1)) {
                            S[t].swapByPos(--S[t].c, j); j--;
                        }

                    } 
                }

                cnts++;
                // assert(cnts < 100);
                if (S[z].r >= sz[z] && mSize < S[z].r * S[t].c) {
                    mSize = S[z].r * S[t].c;
                    for (uint32_t j = 0; j < S[t].c; ++j)
                        visited[S[t][j]] = true;
                }
            }

            // for (uint32_t j = 0; j < S[t].c; ++j)
            //     visited[S[t][j]] = true;

            // while ( S[z].cSize() > 0 && S[t].cSize() > 0 && S[t].c *S[z].c > mSize) {
            //     nextv = g->n[z]+1; mdv = 0; mdc = 0;
            //     if (u == 213855)
            //     printf("in\n"); fflush(stdout);
            //     for (int i = S[z].r; i < S[z].c; ++i) {
            //         uint32_t v = S[z][i];
            //         uint32_t dv = 0;
            //         if (g->deg(z,v) <= S[t].cSize()) {
            //             for (uint32_t j = g->p[z][v]; j < g->p[z][v+1]; ++j){
            //                 uint32_t u1 = g->e[z][j];
            //                 uint32_t pu = S[t].pos(u1);
            //                 if (pu >= S[t].r && pu < S[t].c) dv++;
            //             }
            //         }
            //         else {
            //             for (uint32_t j = S[t].r; j < S[t].c; ++j) {
            //                 uint32_t u1 = S[t][j];
            //                 if (g->is_nbr(z, v, u1)) dv++;

            //             } 
            //         }
            //         if (dv >= sz[t] && mSize < dv * S[z].c) {
            //             if (mdc < g->cores[z][v]) {
            //                 nextv = v; mdc = g->cores[z][v]; mdv = dv;
            //             }
            //             else if (dv > mdv) {
            //                 nextv = v; mdv = dv;
            //             }
            //         }
            //         else {
            //             S[z].swapByPos(--S[z].c, i); i--;
            //         }

            //     }
            //     if (nextv >= g->n[z]) assert(S[z].cSize() == 0);
            //     // if (its < 10)
            //         if (u == 213855)
            //         printf("cnts=%d, nextv=%d, S[z].r=%d,S[z].c=%d\n", cnts, nextv, S[z].r, S[z].c);
            //     if (nextv >= g->n[z]) break;
            //     assert(nextv < g->n[z]);
            //     S[z].swapByPos(S[z].r++, S[z].pos(nextv));

            //     if (g->deg(z,nextv) <= S[t].cSize()) {
            //         uint32_t stsz = S[t].r;
            //         for (uint32_t j = g->p[z][nextv]; j < g->p[z][nextv+1]; ++j){
            //             uint32_t u1 = g->e[z][j];
            //             uint32_t pu = S[t].pos(u1);
            //             if (pu >= stsz && pu < S[t].c) {
            //                 S[t].swapByPos(stsz++, pu);
            //             }
            //         }
            //         S[t].c = stsz;
            //     }
            //     else {
            //         for (uint32_t j = S[t].r; j < S[t].c; ++j) {
            //             uint32_t u1 = S[t][j];
            //             if (!g->is_nbr(z, nextv, u1)) {
            //                 S[t].swapByPos(--S[t].c, j); j--;
            //             }

            //         } 
            //     }
            //     // if (its < 10)
            //     if (u == 213855){
            //         if (nextv == 61741) printf("nextv=%d, S[t][1]=%d, g->deg(z,nextv)=%d, S[t].cSize()=%d\n", nextv, S[t][1], g->deg(z,nextv), S[t].cSize());
            //         printf("cnts=%d, nextv=%d, S[t].r=%d,S[t].c=%d\n", cnts, nextv, S[t].r, S[t].c);
            //     }
            //     cnts++;
            //     // assert(cnts < 100);
            //     if (S[z].r >= sz[z] && mSize < S[z].r * S[t].c) mSize = S[z].r * S[t].c;
            // }

            // for (int i = S[z].r; i < S[z].c && S[z].c >= sz[z]; ++i) {
            //     uint32_t v = S[z][i];
            //     bool fg = false;
            //     for (uint32_t j = S[t].r; j < S[t].c; ++j) {
            //         if (!g->is_nbr(z,v,S[t][j])) {fg = true; break;}
            //     }
            //     if (fg) {S[z].swapByPos(--S[z].c, i); i--;}
            // }

            // if (its < 10){
            //     printf("u=%d, S[z].r=%d,S[z].c=%d\n", u, S[z].r, S[z].c);
            //     printf("u=%d, S[t].r=%d,S[t].c=%d\n", u, S[t].r, S[t].c);
            // }else u = 1;
            //test
            // if (u == 213855)
            // for (int i = 0; i < S[z].r ; ++i) {
            //     uint32_t v = S[z][i];
            //     for (uint32_t j = S[t].r; j < S[t].c; ++j){
            //         if (!g->is_nbr(z,v, S[t][j])) {
            //             printf("i=%d,v=%d, j=%d, u1=%d\n", i,v,j,S[t][j]);
            //         }
            //         assert( g->is_nbr(z,v, S[t][j]) == true );
            //     }
            // }

            if (S[t].c >= sz[t] && S[z].c >= sz[t] && mSize < S[t].c *S[z].c)
                mSize = S[t].c *S[z].c;

            if (its >= 500) break;
            its++;

        }

        visited[u] = true;
        // its++;
        if (u == 0) break;
        // if (u == 213855) break;
        // if (its % 1000 == 0) {
        //     printf("u=%d, mSize=%d\n", u, mSize);
        // }
    
        // break;
    }
    }

    printf("mSize=%d\n", mSize);
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Init time:" << duration.count() << "ms" << std::endl;
}

void MBCS::gemSize(uint32_t t, uint32_t tr, uint32_t zr, uint32_t ne) {

    if (tr >= ne+zr) {
        uint32_t a = std::min(ne, tr-zr);
        ne -= a;
        tr -= (a + ne/2);
        zr -= (ne - ne/2);
    }
    else {
        uint32_t a = std::min(ne, zr-tr);
        ne -= a;
        zr -= (a + ne/2);
        tr -= (ne - ne/2);
    }
    if (mSize < tr * zr) {
        mSize = tr * zr;
        realsz[t] = tr;
        realsz[t^1] = zr;
    }
} 

void MBCS::printRes() {
    printf("U:");
    for (auto u : realRes[0])
        printf("%d ", u);
    printf("\nV:");
    for (auto u : realRes[1])
        printf("%d ", u);
    printf("\n"); fflush(stdout);
}

MBCS::MBCS(/* args */)
{
}

MBCS::MBCS(const std::string & filePath, uint32_t ls, uint32_t rs)
{
    auto t1 = std::chrono::steady_clock::now();
    g = new biGraph(filePath, ls, rs);
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    printf("load graph   %ld ms\n", duration.count());fflush(stdout);  

    S[0].resize(g->n[0]);
    S[1].resize(g->n[1]);
    sdeg[0].resize(g->n1);
    sdeg[1].resize(g->n2);
    
    bin[0].resize(g->n[0]);
    bin[1].resize(g->n[1]);

    ws.resize(std::max(g->n1,g->n2));
    vqueue[0].resize(g->n1);
    vqueue[1].resize(g->n2);
    _deg[0].resize(g->n1);
    _deg[1].resize(g->n2);
    sz[0] = ls; sz[1] = rs;
    // exit(1);
}


MBCS::~MBCS()
{
    if (g != NULL) delete g;
}

void MBCS::reSetdeg(uint32_t *c) {
    for (uint32_t t = 0; t <= 1; ++t) {
        uint32_t z = t^1;
        for (uint32_t i = S[t].c; i < c[t]; ++i) {
            uint32_t v = S[t][i];

            if(g->deg(t, v) > c[z]) {
                for(uint32_t j = 0; j < c[z]; j++) {
                    uint32_t w = S[z][j];
                    if(g->is_nbr(t, v, w))
                        _deg[z][w]++;
                }
            }
            else {
                for(uint32_t j = g->p[t][v]; j < g->p[t][v + 1]; j++) {
                    uint32_t w = g->e[t][j];
                    uint32_t pw = S[z].pos(w);
                    if(pw < c[z]) {//in C and R
                        _deg[z][w]++;
                    }
                }
            }

        }
    }
}


void MBCS::reSetdeg(uint32_t t, uint32_t *c) {
    uint32_t z = t^1;
    for (uint32_t i = S[t].c; i < c[t]; ++i) {
        uint32_t v = S[t][i];

        if(g->deg(t, v) > c[z]) {
            for(uint32_t j = 0; j < c[z]; j++) {
                uint32_t w = S[z][j];
                if(g->is_nbr(t, v, w))
                    _deg[z][w]++;
            }
        }
        else {
            for(uint32_t j = g->p[t][v]; j < g->p[t][v + 1]; j++) {
                uint32_t w = g->e[t][j];
                uint32_t pw = S[z].pos(w);
                if(pw < c[z]) {//in C and R
                    _deg[z][w]++;
                }
            }
        }

    }
}

void MBCS::branchNew(uint32_t deep, uint32_t ne) {

    if (S[0].c < sz[0] || S[1].c < sz[1]) return;
    if (mSize >= S[0].c * S[1].c) return;

    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (S[0].cSize() == 0) {
        if (mSize < S[0].r * S[1].c) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
        return;
    }
    else if (S[1].cSize() == 0) {
        if (mSize < S[0].c * S[1].r) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
        }
    
        return;
    }

    uint32_t r[2], c[2];
    r[0] = S[0].r; r[1] = S[1].r;
    c[0] = S[0].c; c[1] = S[1].c;

    for (uint32_t t = 0; t <= 1; ++t) 
    if (S[t].cSize() == 1) {
        uint32_t d = 0;
        uint32_t z = t^1;
        uint32_t u = S[t][S[t].r];
        for (uint32_t i = S[z].r; i < S[z].c; ++i) {
            uint32_t v = S[z][i];
            if (g->is_nbr(t, u, v)) d++;
        }
        if ((S[z].r+d) * (S[t].r+1) > mSize  &&  S[t].r+1 >= sz[t] && S[z].r+d >= sz[z]) {
            // mSize = (S[z].r+d) * (S[t].r+1);
            mSize = (S[z].r+d) * (S[t].r+1); realsz[t] = (S[t].r+1); realsz[z] = (S[z].r+d);
            // assert((S[z].r+d) * (S[t].r+1) != 880);
            
            // gemSize(t,S[t].r+1, S[z].r+d, ne);
        }
        else if (S[z].c * S[t].r > mSize &&  S[t].r >= sz[t] && S[z].c >= sz[z]) {
            // mSize = S[z].c * S[t].r;
            mSize = S[z].c * S[t].r; realsz[t] = S[t].r; realsz[z] = S[z].c;
            // assert(S[z].c * S[t].r != 880);
            // gemSize(t, S[t].r, S[z].c, ne);
        }
        return;
    }

    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};

    uint32_t t = 0, z = 1;
    uint32_t sube = 0;
    if (S[1].cSize() < S[0].cSize()) {t = 1; z = 0;} 

    if (S[t].c == sz[t]) {

        for (int i = S[z].r; i < S[z].c; ++i) {
            uint32_t u = S[z][i];
            for (uint32_t j = S[t].r; j < S[t].c; ++j) {
                uint32_t v = S[t][j];
                if (!g->is_nbr(z,u,v)) {
                    S[z].swapByPos(i--, --S[z].c);
                    break;
                }

            }
            if (S[t].c * S[z].c <= mSize) break;
        }

        S[t].r = S[t].c;
        branchNew(deep+1, ne);
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].c = c[z];
        return;
    }
    if (S[z].c == sz[z]) {

        for (int i = S[t].r; i < S[t].c; ++i) {
            uint32_t u = S[t][i];
            for (uint32_t j = S[z].r; j < S[z].c; ++j) {
                uint32_t v = S[z][j];
                if (!g->is_nbr(t,u,v)) {
                    S[t].swapByPos(i--, --S[t].c);
                    break;
                }

            }
            if (S[t].c * S[z].c <= mSize) break;
        }

        S[z].r = S[z].c;
        branchNew(deep+1, ne);
        S[z].r = r[z];
        S[z].c = c[z];
        S[t].c = c[t];
        return;
    }


#ifdef UB1_
        uint32_t qsz[2] = {0};
        for(uint32_t j = S[z].r; j < c[z]; j++) _deg[z][S[z][j]] = 0;
        for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
        // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
            uint32_t u = S[t][i];
            uint32_t e = 0;
            if(g->deg(t, u) > 2*S[z].cSize()) {
                for(uint32_t j = S[z].r; j < S[z].c; j++) {
                    uint32_t v = S[z][j];
                    if(g->is_nbr(t, u, v)) {
                        e++; _deg[z][v]++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    uint32_t pv = S[z].pos(v);
                    if(S[z].r <= pv && pv < S[z].c) {//in C
                        e++; _deg[z][v]++;
                    }
                }
            }
            
            _deg[t][u] = e;
            if (S[z].r+e < sz[z] || S[t].c * (S[z].r+e) <= mSize) {
                S[t].swapByPos(--S[t].c, i--);
                vqueue[t][qsz[t]++] = u;
                if (S[t].c < sz[t] || S[t].c * S[z].c <= mSize)  {
                    S[t].r = r[t];
                    S[t].c = c[t];
                    S[z].r = r[z];
                    S[z].c = c[z];
                    return;
                }
            }
        }

        for(int j = S[z].r; j < S[z].c; j++) {
            uint32_t u = S[z][j];
            uint32_t e = _deg[z][u];
            if (r[t]+e < sz[t] || S[z].c *(r[t]+e) < mSize) {
                S[z].swapByPos(--S[z].c, j--);
                vqueue[z][qsz[z]++] = u;
            }
        }

        // if (deep == 1 && S[z][0] == 26285) {
        //     printf("sz[t]=%d, sz[z]=%d\n", sz[t], sz[z]);
        //     for (uint32_t i = S[t].r; i < S[t].c; ++i) {
        //         printf("t=%d, v=%d, deg[v]=%d, _deg[v]=%d\n", t, S[t][i], _deg[t][S[t][i]]+S[z].r, S[z].c - S[z].r - _deg[t][S[t][i]]);
        //     }
        // }

        uint32_t qr[2] = {0};
        while(qr[0] < qsz[0] || qr[1] < qsz[1]) {
            for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                uint32_t zzz = ttt^1;
                for (uint32_t i = qr[ttt]; i < qsz[ttt]; ++i) {
                    uint32_t u = vqueue[ttt][i];
                    _deg[ttt][u] = 0;

                    if(g->deg(ttt, u) > 10*S[zzz].cSize()) {
                        for(int j = S[zzz].r; j < S[zzz].c; j++) {
                            uint32_t v = S[zzz][j];
                            if(g->is_nbr(ttt, u, v)) {
                                uint32_t d = --_deg[zzz][v];
                                d += S[ttt].r;
                                if (d < sz[ttt] || d * S[zzz].c <= mSize) {
                                    S[zzz].swapByPos(--S[zzz].c, j--);
                                    vqueue[zzz][qsz[zzz]++] = v;
                                    if (S[t].c < sz[t] || S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  break;
                                }
                            }
                        }
                    }
                    else 
                    for (uint32_t j = g->p[ttt][u]; j < g->p[ttt][u+1]; ++j) {
                        uint32_t v = g->e[ttt][j];
                        uint32_t pv = S[zzz].pos(v);
                        if (pv < S[zzz].c && S[zzz].r <= pv) {
                            uint32_t d = --_deg[zzz][v];
                            d += S[ttt].r;
                            if (d < sz[ttt] || d * S[zzz].c <= mSize) {
                                S[zzz].swapByPos(--S[zzz].c, pv);
                                vqueue[zzz][qsz[zzz]++] = v;
                                if (S[t].c < sz[t] || S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  break;
                            }
                        }
                    }
                    if (S[t].c < sz[t] || S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  {
                        i = qsz[ttt]; qr[zzz] = qsz[zzz];
                        break;
                    }
                }
                qr[ttt] = qsz[ttt];

            }
        }
        if (S[t].c < sz[t] || S[z].c < sz[z]|| S[0].c * S[1].c <= mSize)  {
            S[t].c = c[t]; S[t].r = r[t]; 
            S[z].c = c[z]; S[z].r = r[z];
            return;
        }

        maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
        uint32_t minI = 0, minU = 0, mnume = 0, minE = std::max(S[t].c, S[z].c)+1;

        mnume = S[t].r*S[z].c+S[z].r*(S[t].c-S[t].r);
        if (S[t].cSize() > S[z].cSize()) {
            t = z; z = t^1;
        }

        qsz[t] = 0; qsz[z] = 0;
        for(uint32_t j = S[t].r; j < S[t].c; j++) {
            uint32_t u = S[t][j];
            uint32_t e = _deg[t][u];

            if (e == S[z].c - r[z])  {
                vqueue[t][j] = vqueue[t][S[t].r];
                S[t].swapByPos(S[t].r++, j);
                continue;
            }
            else vqueue[t][j] = S[z].c - r[z] - e;


            if(e > maxE[t]) {
                maxE[t] = e;
                maxI[t] = j;
                maxUV[t] = u;
            }
            if (e < minE) {
                minE = e;
                minI = j;
                minU = u;
            }
            mnume += e;
        }

        for(int j = S[z].r; j < S[z].c; j++) {
            uint32_t u = S[z][j];
            uint32_t e = _deg[z][u]; 

            if (e == S[t].c - r[t]) {
                vqueue[z][j] = vqueue[z][S[z].r];
                S[z].swapByPos(S[z].r++, j);
                continue;
            }else vqueue[z][j] = S[t].c - r[t] - e;

            if(e > maxE[z]) {
                maxE[z] = e;
                maxI[z] = j;
                maxUV[z] = u;
            }
        }

        if (S[t].cSize() > 0 || S[z].cSize() > 0){


            // std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c, [](int a, int b){return a > b;});
            // std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c, [](int a, int b){return a > b;});

            std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c);
            std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c);

            uint32_t mxnbz = 0;
            for (uint32_t i = S[z].r; i < S[z].c; ++i) {
                uint32_t nbz = vqueue[z][i];
                bin[z][nbz]++;
                mxnbz = std::max(mxnbz, nbz);
            }

            uint32_t ubm = 0, zsz = S[z].c, nonedges[2] = {0};
            for (uint32_t i = S[t].r; i < S[t].c; ++i) {
                uint32_t _d = vqueue[t][i];
                uint32_t d = S[z].c - _d;
                nonedges[t] += _d;

                while (zsz > S[z].r && nonedges[t] > nonedges[z]) {
                    zsz--;
                    uint32_t _dz = vqueue[z][zsz];
                    nonedges[z] += _dz;
                    if (zsz == S[z].r) break;
                }
                while (zsz > S[z].r && vqueue[z][zsz-1] + i + 1 > S[t].c) {
                    nonedges[z] += vqueue[z][--zsz];
                    if (zsz == S[z].r) break;
                }

                if (zsz > S[z].r && i-S[t].r >= mxnbz) {
                    uint32_t j = i - mxnbz;
                    uint32_t tail = vqueue[t][j] + vqueue[t][j+1];
                    while (zsz > S[z].r && tail+zsz > S[z].c)
                        nonedges[z] += vqueue[z][--zsz];
                }

                uint32_t ubz = std::min(d, zsz);
                if (ubz < sz[z]) break;
                if (i+1 >= sz[t]) {
                    ubm = std::max(ubm, (i+1) * ubz);
                }
            }
        
            for (uint32_t i = 0; i <= mxnbz; ++i) bin[z][i] = 0;

            if (ubm <= mSize) {
                S[t].c = c[t]; S[t].r = r[t]; 
                S[z].c = c[z]; S[z].r = r[z];
                return;
            }

        }
#else
    // r[0] = S[0].r; r[1] = S[1].r;
    // c[0] = S[0].c; c[1] = S[1].c;

    // sube = 0;

    maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
    uint32_t minI = 0, minU = 0, minE = S[z].c+1;
    for(int i = S[t].r; i < S[t].c; i++) {//scan all vertices in C
    // for(uint32_t i = S[t].c; i <= S[t].c; i++) {
        uint32_t u = S[t][i];
        uint32_t e = 0;
        if(g->deg(t, u) > 2*S[z].cSize()) {
            for(uint32_t j = S[z].r; j < S[z].c; j++) {
                uint32_t v = S[z][j];
                if(g->is_nbr(t, u, v)) {
                    e++; _deg[z][v]++;
                }
            }
        }
        else {
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t v = g->e[t][j];
                uint32_t pv = S[z].pos(v);
                if(S[z].r <= pv && pv < S[z].c) {//in C
                    e++; _deg[z][v]++;
                }
            }
        }

        assert(e == _deg[t][u]);

        // if (e == S[z].cSize()) {
        //     S[t].swapByPos(S[t].r++, i);
        //     continue;
        // }

        if(e > maxE[t]) {
            maxE[t] = e;
            maxI[t] = i;
            maxUV[t] = u;
            // if (e >= S[z].cSize()) break;
        }
        if (e < minE) {
            minE = e;
            minI = i;
            minU = u;
        }
        // _deg[t][u] = e;
        // sube += e;
    }

    for(int i = S[z].r; i < S[z].c; i++) {
        uint32_t u = S[z][i];
        uint32_t e = _deg[z][u];
        _deg[z][u] = 0;
        // if (e == S[t].c - r[t]) {
        //     // uint32_t d = 0;
        //     // for (uint32_t j = S[t].r; j  < S[t].c; ++j)
        //     //     if (g->is_nbr(z,u, S[t][j])) d++;
        //     // assert(d == e);
        //     // assert(i == S[z].pos(u));
        //     S[z].swapByPos(S[z].r++, i);
        // }
        // else 
        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[z] = i;
            maxUV[z] = u;
        }
    }
#endif

    if (S[0].cSize() == 0) {

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);

            // for (uint32_t v = 0; v < S[0].r; ++v){
            //     for (uint32_t u = 0; u < S[1].c; ++u)
            //         assert(g->is_nbr(0, S[0][v], S[1][u]));
            // }
            // printRes();
        }
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].r = r[z];
        S[z].c = c[z];
        return;
    }
    else if (S[1].cSize() == 0) {

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            // gemSize(0, S[0].c, S[1].r, ne);
            mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
            // for (uint32_t v = 0; v < S[0].r; ++v){
            //     for (uint32_t u = 0; u < S[1].c; ++u)
            //         assert(g->is_nbr(0, S[0][v], S[1][u]));
            // }
            // printRes();
        }
  
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].r = r[z];
        S[z].c = c[z];
        return;
    }
    
    // printf("maxE[0]=%d, maxE[1]=%d\n", maxE[0],maxE[1]);
    // exit(1);
    if (maxE[t] == 0) {

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            // gemSize(0, S[0].r, S[1].c, ne);
            mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
            // printRes();
        }

        if (mSize < S[0].c * S[1].r && S[0].c >= sz[0] && S[1].r >= sz[1]) {
            gemSize(0, S[0].c, S[1].r, ne);
            // mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
            // printRes();
        }
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].r = r[z];
        S[z].c = c[z];
        return;
    }

    assert(S[0].cSize() > 0);
    assert(S[1].cSize() > 0);

    if (BIBRANCH && minE + 4 < S[z].cSize() && S[t].cSize() * 10 < S[z].cSize()) {

        uint32_t u = minU, cz = S[z].c;
        assert( S[t].pos(u) < S[t].c);
        S[t].swapByPos(S[t].r++, S[t].pos(u));

        if(g->deg(t, u) > 2*S[z].cSize()) {
            for(int j = S[z].r; j < S[z].c; j++) {
                if(!g->is_nbr(t, u, S[z][j])) {
                    S[z].swapByPos(--S[z].c, j); j--;
                }
            }
        }
        else {
            uint32_t cc = S[z].r;
            for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                uint32_t w = g->e[t][j];
                uint32_t pw = S[z].pos(w);
                if(S[z].r <= pw && pw < S[z].c) {
                    S[z].swapByPos(cc++, pw);
                }
            }
            S[z].c = cc;
        }

        branchNew(deep+1, ne);
        S[z].c = cz;
        assert(S[t][S[t].r-1] == u);
        S[t].swapByPos(--S[t].r, --S[t].c);

        branchNew(deep+1, ne);
        S[t].r = r[t];
        S[t].c = c[t];
        S[z].r = r[z];
        S[z].c = c[z];
        return;
    }

    // if (maxE[t]+1 >= S[z].cSize()) {
    //     if (maxE[t] == S[z].cSize()) {
    //         S[t].swapByPos(S[t].r++, S[t].pos(maxUV[t]));
    //         for(uint32_t i = S[z].r; i < S[z].c; i++) 
    //             subdeg[S[z][i]] = 0;
    //         branchNew(deep+1, ne);
    //         S[t].r--;
    //     }
    //     else {
            
    //         for(uint32_t i = S[z].r; i < S[z].c; i++) 
    //             subdeg[S[z][i]] = 0;
    //     }
    //     return;
    // }
    
    // for(uint32_t i = S[z].r; i < S[z].c; i++) {
    //     uint32_t u = S[z][i];
    //     uint32_t e = subdeg[u];
    //     subdeg[u] = 0;

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[z] = i;
    //         maxUV[z] = u;
    //     }
    // }

    // t = 0, z = 1;
    // if (S[z].cSize() - maxE[t] >= S[t].cSize()-maxE[z]) {
    //     t = z; z = t^1;
    // }
    if (S[z].c - r[z] - maxE[t] >= S[t].c - r[t]-maxE[z]) {
        t = z; z = t^1;
    }
    // printf("t=%d, z=%d\n", t,z);
    // exit(1);

    uint32_t u = maxUV[t], wsSize = 0;

    if(g->deg(t, u) > 2*S[z].cSize()) {
        if (S[z].cSize() > ws[deep].capacity()) ws[deep].resize(S[z].cSize()+1);
        for(uint32_t i = S[z].r; i < S[z].c; i++) {
            uint32_t w = S[z][i];
            if(!g->is_nbr(t, u, w)) {
                // if ((_deg[z][w]+r[t]) >= mSize / (_deg[z][w]+r[t]))
                    ws[deep][wsSize++] = S[z][i];
            }
        }
    }
    else {
        uint32_t cc = S[z].c;
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].r <= pv && pv < S[z].c) {
                // if ((_deg[z][v]+r[t]) >= mSize / (_deg[z][v]+r[t]))
                    S[z].swapByPos(--cc, pv);
            }
        }
        wsSize = cc - S[z].r;
        if (wsSize > ws[deep].capacity()) ws[deep].resize(wsSize+1);
        memcpy(ws[deep].data(), S[z].begin() + S[z].r, sizeof(uint32_t) * wsSize);
    }

    maxsubbranches = wsSize+1 > maxsubbranches ? wsSize+1 : maxsubbranches;
    fflush(stdout); 
    uint32_t ct = S[t].c, cz = S[z].c;
    uint32_t dv = _deg[z][ws[deep][0]] + r[t];
    // if (wsSize <= 1) subbcnt++;
    // std::sort(wsSize[deep].begin(), wsSize[deep].begin()+wsSize, [](int a, int b){return a > b;});
    std::sort(ws[deep].begin(), ws[deep].begin()+wsSize, [&](int a, int b){return _deg[z][a] < _deg[z][b];});
    for (uint32_t i = 0; i < wsSize; ++i) {

        uint32_t v = ws[deep][i];
        // if (deep == 0) {
        //     printf("z=%d, v=%d, i=%d\n", z, v, i);fflush(stdout); 
        //     printf("ws:");
        //     for ( auto i = 0; i < wsSize; ++i) 
        //         printf("%d ", ws[deep][i]);
        //     printf("\n");

        //     if (i >=12 ) exit(1);
        // }
        S[z].swapByPos(S[z].r++, S[z].pos(v));
        // if (g->deg(z, v) < 0) printf("z=%d, v=%d, g->n[z]=%d, g->deg(z, v)=%d\n", z, v, g->n[z], g->deg(z, v)); fflush(stdout);  
        if(g->deg(z, v) > S[t].cSize()) {
            for(int j = S[t].r; j < S[t].c; j++) {
                // assert(v >= 0);
                // if (v >= g->n[z]) {
                //     printf("z=%d, v=%d, g->n[z]=%d, g->deg(z, v)=%d\n", z, v, g->n[z], g->deg(z, v)); fflush(stdout);  
                //     printf("deep=%d, r=%d, wsSize=%d, c=%d, t=%d, z=%d, i=%d\n", deep, S[z].r, wsSize, S[z].c, t, z, i);

                //     printf("ws:");
                //     for ( auto i = 0; i < wsSize; ++i) 
                //         printf("%d ", ws[deep][i]);
                //     printf("\n");
                //     fflush(stdout);  
                // }
                // assert(v < g->n[z]);
                if(!g->is_nbr(z, v, S[t][j])) {
                    S[t].swapByPos(--S[t].c, j); j--;
                }
            }
        }
        else {
            uint32_t cc = S[t].r;
            for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                uint32_t w = g->e[z][j];
                uint32_t pw = S[t].pos(w);
                if(S[t].r <= pw && pw < S[t].c) {
                    S[t].swapByPos(cc++, pw);
                }
            }
            S[t].c = cc;
        }

        uint32_t szr = S[z].r;
        branchNew(deep+1, ne);
        S[t].c = ct;
        assert(S[z].c == cz - i);
        if (S[z][S[z].r-1] != v) {
            S[0].print();
            S[1].print();
            printf("z=%d, v=%d, S[z][S[z].r-1]=%d, szr=%d, S[z].r=%d, S[t].r=%d, S[t].c=%d\n", z, v, S[z][S[z].r-1], szr, S[z].r, S[t].r, S[t].c);
        }
        assert(S[z][S[z].r-1] == v);
        S[z].swapByPos(--S[z].r,--S[z].c);
    
        // if (S[0][0] == 12) {
        //     exit(-1);
        // }
        // S[t].c = c[t];
        // S[z].swapByPos(--S[z].r,--S[z].c);

        // if (deep == 2) {
        //     printf("z=%d, v=%d, done\n", z, v);fflush(stdout); 
        // }
    }

    if (wsSize == 1) {
        uint32_t x = S[z].c - 1;
        uint32_t y = mSize / x;

        if (y >= (x+1)*(S[z].c - dv)) {
            S[t].r = r[t];
            S[t].c = c[t];
            S[z].r = r[z];
            S[z].c = c[z];
            return;
        }
    }
    // pivot vertex
    S[t].swapByPos(S[t].r++, S[t].pos(u));

    branchNew(deep+1, ne);

    S[t].r = r[t];
    S[t].c = c[t];
    S[z].r = r[z];
    S[z].c = c[z];

}