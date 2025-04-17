#include "MBCS.h"

#define RUN_1

#define UB1_
#define UB_HNBR_
#define _BASIC_
// #define _REVERSE_

#define REDUCTION_

#define INIT_SIZE

#ifdef RUN_1

void MBCS::run(int alg) {
    
    float dalpha = 0.3; 

    #if defined(UB_HNBR_)
    printf("run: defined_UB_HNBR_\n");
    #elif defined(UB1_)
    printf("run: undefined_UB_HNBR_but_defined_UB1\n");
    #else
    printf("run: undefined_UB_HNBR_and_UB1\n");
    #endif

    // for(uint32_t i = 0; i < g->n[1]; ++i) {
    //     uint32_t d = g->p[1][i+1] - g->p[1][i];
    //     if (d >= 7591)
    //         printf("V v=%d, d=%d, is_nbr(12,53)=%d\n", i, d, g->is_nbr(0,12,53));
    // }
    S[0].c = g->n[0]; S[1].c = g->n[1];
    this->alg = alg;
    // mSize=10212;
    // mSize= 205180;
    // mSize = 880;

    #ifdef INIT_SIZE
    initSize();
    #endif

    // mSize = 18960;

    // if (mSize >= std::min(g->maxDu, g->maxDu) * 100){
    //     #define BIBRANCH_
    // }

    Sc[0].resize(g->n[0]);
    Sc[1].resize(g->n[1]);

    uint32_t qsz[2] = {0};
    uint32_t t = 0, z = 1;
    std::vector<uint32_t> orders;
    if (g->maxDu < g->maxDv) {t = 1; z = 0;} 
    orders.resize(g->n[t]+1);
    if (g->maxcore>=45) {
        // stp = 1.2;
        // for (uint32_t i = 0; i < g->n[t]; ++i) orders[i] = g->n[t]-1-i;
            for (uint32_t i = 0; i < g->n[t]; ++i) orders[i] = i;
    }
    else {
        stp = 2;
        for (uint32_t i = 0; i < g->n[t]; ++i) orders[i] = i;
    }

    // return;
    printf("g->n[0]=%d, g->n[1]=%d, mSize=%d, t=%d, z=%d, stp=%.3f\n", g->n[0], g->n[1], mSize, t, z, stp);
    // printf("t=%d, u=21878, realu=%d\n", t, g->realLable[t][21878]); return;

    unsigned long maxBranchduration = 0;
    uint32_t crsz[2] = {0}, tt = 0, zz = 1;

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

    // #ifdef _REVERSE_
    // for (uint32_t u = std::max((uint32_t)0, g->n[t]-crsz[t]); u >= 0; u--) {
    // #else
    // for (uint32_t u = 0; u < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); u++) {
    // #endif
    // for (uint32_t u = 61458; u < 61458+1 && u < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); u++) {

    for (uint32_t u, ui = 0; ui < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); ui++) {
        u = orders[ui];
        S[t].r = S[t].c = 0;
        S[z].r = S[z].c = 0;
        if (S[t].capacity() <= u) {
            printf("S[t].capacity()=%d, u=%d, std::max((uint32_t)0, g->n[t]-crsz[t])=%d\n", S[t].capacity(), u, std::max((uint32_t)0, g->n[t]-crsz[t]));
        }

        assert(S[t].capacity() > u);
        assert(S[t][S[t].pos(u)] == u);

        uint32_t Szs = g->deg(t,u);
        uint32_t tc = 0, zc = 0;

        if (Szs * g->mxd[z] <= mSize || Szs < crsz[z]) {
            #ifdef _REVERSE_
            if (u == 0) break;
            #endif
            continue;
        }
        S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

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

            if (Szs * nbrv <= mSize || nbrv < crsz[t]) {
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
                if (j == 0) break;
            }
        }

        tc = S[t].c; zc = S[z].c;

        ws[0].resize(0);
        // reduction

//         sz[t] = crsz[t]; sz[z] = crsz[z];
// if (alg == 4 && (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize)) {

//     uint32_t e = 0;
//     for(uint32_t v = 0; v < S[z].c; ++v)
//         e += sdeg[z][S[z][v]];
//     if (double(e / (S[t].c * S[z].c)) < 0.25) {
//         subGraph_generate();
//         branchSubG(1,0);
//     }
//     else {
//         lrs[t] = S[t].c; lrs[z] = S[z].c;
//         graph_convert();
//         csz[t] = lrs[t] - sz[t]; csz[z] = lrs[z] - sz[z];
//         branch_minCover(u,1,0);
//     }
// }
// else {

        sz[tt] = crsz[tt]; sz[zz] = g->mxd[tt] / 4;
        sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt]));

        it = 0; Sc[t].c = 0; Sc[z].c = 0;
        unsigned long old_numbranches = numbranches;
        while (sz[tt] <= g->mxd[zz] && sz[zz] >= crsz[zz]) {

            if (it > 0) {
                sz[tt] = sz[tt] + 1;
                sz[zz] = std::max(crsz[zz], mSize / sz[tt]);
            }
            if (sz[tt] > g->mxd[zz] || it >= 1) break;

            old_numbranches = numbranches;
            if (alg == 3){sz[t] = crsz[t]; sz[z] = crsz[z];}

            S[t].c = tc; S[z].c = zc;
            qsz[0] = 0; qsz[1] = 0;
            for (int j = 0; j < S[t].c; ++j){
                uint32_t v = S[t][j];
                uint32_t dv = sdeg[t][v];
                _deg[t][v] = dv;
            }
            assert(u == S[t][0]);
            for (uint32_t j = 0; j < S[z].c; ++j){
                uint32_t v = S[z][j];
                uint32_t dv = sdeg[z][v];
                _deg[z][v] = sdeg[z][v];
            }

            qsz[0] = 0; qsz[1] = 0;
            for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                uint32_t zzz = ttt^1;
                uint32_t v, dv;
                for (int j = 0 ; j < S[ttt].c; ++j) {
                    v = S[ttt][j];
                    dv = _deg[ttt][v];
                    if (dv < sz[zzz] || dv * S[ttt].c <= mSize) {
                        vqueue[ttt][qsz[ttt]++] = v;
                        if ((v == u && ttt == t) || S[ttt].c < sz[ttt] || S[t].c * S[z].c <= mSize) {
                            qsz[ttt] = 0; qsz[zzz] = 0;
                            S[ttt].c = 0;
                            S[zzz].c = 0;
                            break;
                        }
                        S[ttt].swapByPos(--S[ttt].c, j--);
                    }
                }
            }
            // if (u == 15432 && it ==1) printf("u=%d, _deg[t][u]=%d, qsz[z]=%d, S[z].c=%d\n", u ,  _deg[t][u], qsz[z], S[z].c);
            if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
                uint32_t qr[2] = {0};
                while (qr[0] < qsz[0] || qr[1] < qsz[1]) {
                    for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                        uint32_t zzz = ttt^1;
                        for (uint32_t j = qr[ttt]; j < qsz[ttt]; ++j) {
                            uint32_t v = vqueue[ttt][j];
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
                                            S[ttt].c = 0;
                                            S[zzz].c = 0;
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
                if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
                        // printf("u=%d, it=%d, S[t].r=%d, S[t].c=%d, S[z].r=%d, S[z].c=%d\n", u, it, S[t].r, S[t].c, S[z].r, S[z].c);
                        uint32_t dens = 0, cnts = 0, wsSize = 0;
                        for (uint32_t j = S[t].r; j < S[t].c; ++j) {
                            uint32_t v = S[t][j];
                            dens += _deg[t][v]; 
                            // uint32_t vp = Sc[t].pos(v);
                            // if (vp >= Sc[t].c) {
                            //     Sc[t].swapByPos(Sc[t].c++, vp); cnts++;
                            //     // ws[0].push_back(v);
                            // }
                        }
                        // graph_convert();
                        assert(gt == tt);
                        // BIBRANCH = false;
                        if (S[t].c / 100 >= S[z].c || S[z].c / 100 >= S[t].c) BIBRANCH = true;
                        else if (std::max(S[t].c, S[z].c) > 10000 && std::max(S[t].c, S[z].c) / std::min(S[t].c, S[z].c) >= 30) BIBRANCH = true;
                        else BIBRANCH = false;
                        uint32_t mSize1 = mSize;
                        if (alg == 1) {
                            subGraph_generate();
                            branchSubG(1,0);

                            // branchNew(1, 0);

                            // lrs[t] = S[t].c; lrs[z] = S[z].c; Sc[t].c = 0; Sc[z].c = 0;
                            // branch_basic(1);
                        }
                        else if (alg == 2 || alg == 3) {

                            // if (Sc[t].r > 0 && S[z].c > 0 && (double(dens)/double(Sc[t].c*S[z].c)) >= 0.8 && (double(dens)/double(Sc[t].c*S[z].c)) <= 0.9) {
                            if (S[t].c > 0 && S[z].c > 0 && (double(dens)/double(S[t].c*S[z].c)) >= dalpha ) {
                                lrs[t] = S[t].c; lrs[z] = S[z].c;
                                graph_convert();
                                csz[t] = lrs[t] - sz[t]; csz[z] = lrs[z] - sz[z];
                                branch_minCover(u,1,0);
                            }
                            else if (S[t].c > 0 && S[z].c > 0) {
                                subGraph_generate();
                                branchSubG(1,0);
                                // branchNew(1, 0);

                            //     lrs[t] = S[t].c; lrs[z] = S[z].c; Sc[t].c = 0; Sc[z].c = 0;
                            // branch_basic(1);
                            }
                        }
                        // Sc[t].r = Sc[t].c;

                        if (numbranches - old_numbranches > maxBranchduration || numbranches - old_numbranches >= 10000000) {
                            printf("it=%d, u=%d, sz[t]=%d, sz[z]=%d, ct=%d, cz=%d, maxBranchduration=%d, mSize=%d, density=%.3f\n", it, u, sz[t], sz[z], S[t].c, S[z].c, numbranches - old_numbranches, mSize, S[t].c*S[z].c<=0?0:double(dens)/double(S[t].c*S[z].c));
                        }
                        if (mSize1 < mSize) printf("it=%d, u=%d, sz[t]=%d, sz[z]=%d, ct=%d, cz=%d, maxBranchduration=%d, mSize=%d, density=%.3f\n", it, u, sz[t], sz[z], S[t].c, S[z].c, numbranches - old_numbranches, mSize, S[t].c*S[z].c<=0?0:double(dens)/double(S[t].c*S[z].c));
                        maxBranchduration = std::max(maxBranchduration, numbranches - old_numbranches);
                    }

                    // for (uint32_t j = 0; j < S[t].c; ++j)
                    //     _deg[t][S[t][j]] = 0;
                    // for (uint32_t j = 0; j < S[z].c; ++j)
                    //     _deg[z][S[z][j]] = 0;
            }

            if (sz[zz] == crsz[zz] || sz[tt] == g->mxd[zz]) break;
            it++;
            // if (it > 20) { printf("error\n"); exit(1); }
        }
// }
        for (uint32_t j = 0; j < tc; ++j){
            uint32_t v = S[t][j];
            sdeg[t][v] = 0;
        }
        for (uint32_t j = 0; j < zc; ++j){
            uint32_t v = S[z][j];
            sdeg[z][v] = 0;
        }

        #ifdef _REVERSE_
            if (u == 0) break;
        #endif
    
    }

    if (mSize / crsz[zz] <= crsz[tt] || alg == 3) {
        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        printf("mSize=%d, lsz=%d,rsz=%d, mxdeep=%d, numbranches=%ld, maxsubbranches=%ld, maxBranchduration=%d, subbcnt=%ld\n", mSize, realsz[0], realsz[1], mxdeep, numbranches, maxsubbranches, maxBranchduration, subbcnt);
        std::cout << "iteration time:" << duration.count() << "ms" << std::endl;
        printRes();
        return;  
    }

    for (uint32_t u, ui = 0; ui < std::min(g->n[t]-1, std::max((uint32_t)0,g->n[t]-crsz[t])); ui++) {
        u = orders[ui];
        S[t].r = S[t].c = 0;
        S[z].r = S[z].c = 0;
        if (S[t].capacity() <= u) {
            printf("S[t].capacity()=%d, u=%d, std::max((uint32_t)0, g->n[t]-crsz[t])=%d\n", S[t].capacity(), u, std::max((uint32_t)0, g->n[t]-crsz[t]));
        }

        assert(S[t].capacity() > u);
        assert(S[t][S[t].pos(u)] == u);

        uint32_t Szs = g->deg(t,u);
        uint32_t tc = 0, zc = 0;

        if (Szs * g->mxd[z] <= mSize || Szs < crsz[z]) {
            #ifdef _REVERSE_
            if (u == 0) break;
            #endif
            continue;
        }
        S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

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

            if (Szs * nbrv <= mSize || nbrv < crsz[t]) {
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
                if (j == 0) break;
            }
        }

        tc = S[t].c; zc = S[z].c;

        ws[0].resize(0);
        // reduction

        sz[tt] = crsz[tt]; sz[zz] = g->mxd[tt] / 4;
        sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt]));

        it = 1; Sc[t].c = 0; Sc[z].c = 0;
        unsigned long old_numbranches = numbranches;
        while (sz[tt] <= g->mxd[zz] && sz[zz] >= crsz[zz]) {

            // if (it > 0 && it < 1) {
            //     sz[tt] = sz[tt] + 1;
            //     sz[zz] = std::max(crsz[zz], mSize / sz[tt]);
            // }
            if (it == 1) {
                sz[tt] = sz[tt] + 1;
                sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt]/ stp));
            }
            else if (it > 1) {
                if (uint32_t(sz[zz] / stp) < crsz[zz]) {
                    sz[tt] = sz[tt] * sz[zz] / crsz[zz];
                    sz[zz] = crsz[zz];
                }
                else {
                    if (uint32_t(sz[tt] * stp) == sz[tt]) {
                        sz[tt] = sz[tt] + 1;
                        sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt] / stp));
                    }
                    else {
                        sz[tt] = sz[tt] * stp;
                        sz[zz] = std::max(crsz[zz], uint32_t(mSize / sz[tt] / stp));
                    }
                }
            }

            if (sz[tt] > g->mxd[zz]) break;

            old_numbranches = numbranches;

            S[t].c = tc; S[z].c = zc;
            qsz[0] = 0; qsz[1] = 0;
            for (int j = 0; j < S[t].c; ++j){
                uint32_t v = S[t][j];
                uint32_t dv = sdeg[t][v];
                _deg[t][v] = dv;
            }
            assert(u == S[t][0]);
            for (uint32_t j = 0; j < S[z].c; ++j){
                uint32_t v = S[z][j];
                uint32_t dv = sdeg[z][v];
                _deg[z][v] = sdeg[z][v];
            }

            qsz[0] = 0; qsz[1] = 0;
            for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                uint32_t zzz = ttt^1;
                uint32_t v, dv;
                for (int j = 0 ; j < S[ttt].c; ++j) {
                    v = S[ttt][j];
                    dv = _deg[ttt][v];
                    if (dv < sz[zzz] || dv * S[ttt].c <= mSize) {
                        vqueue[ttt][qsz[ttt]++] = v;
                        if ((v == u && ttt == t) || S[ttt].c < sz[ttt] || S[t].c * S[z].c <= mSize) {
                            qsz[ttt] = 0; qsz[zzz] = 0;
                            S[ttt].c = 0;
                            S[zzz].c = 0;
                            break;
                        }
                        S[ttt].swapByPos(--S[ttt].c, j--);
                    }
                }
            }
            // if (u == 15432 && it ==1) printf("u=%d, _deg[t][u]=%d, qsz[z]=%d, S[z].c=%d\n", u ,  _deg[t][u], qsz[z], S[z].c);
            if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c > mSize) {
                uint32_t qr[2] = {0};
                while (qr[0] < qsz[0] || qr[1] < qsz[1]) {
                    for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                        uint32_t zzz = ttt^1;
                        for (uint32_t j = qr[ttt]; j < qsz[ttt]; ++j) {
                            uint32_t v = vqueue[ttt][j];
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
                                            S[ttt].c = 0;
                                            S[zzz].c = 0;
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
                if (S[t].c >= sz[t] && S[z].c >= sz[z] && S[t].c * S[z].c >= mSize) {
                        // printf("u=%d, it=%d, S[t].r=%d, S[t].c=%d, S[z].r=%d, S[z].c=%d\n", u, it, S[t].r, S[t].c, S[z].r, S[z].c);
                        uint32_t dens = 0, cnts = 0, wsSize = 0;
                        for (uint32_t j = S[t].r; j < S[t].c; ++j){
                            uint32_t v = S[t][j];
                            dens += _deg[t][v]; 
                            // uint32_t vp = Sc[t].pos(v);
                            // if (vp >= Sc[t].c) {
                            //     Sc[t].swapByPos(Sc[t].c++, vp); cnts++;
                            //     // ws[0].push_back(v);
                            // }
                        }
                        // graph_convert();
                        assert(gt == tt);
                        // BIBRANCH = false;
                        if (S[t].c / 100 >= S[z].c || S[z].c / 100 >= S[t].c) BIBRANCH = true;
                        else BIBRANCH = false;
                        uint32_t mSize1 = mSize;
                        if (alg == 1) {
                            subGraph_generate();
                            branchSubG(1,0);
                            // branchNew(1, 0);

                            // lrs[t] = S[t].c; lrs[z] = S[z].c; Sc[t].c = 0; Sc[z].c = 0;
                            // branch_basic(1);
                        }
                        else if (alg == 2) {

                            if (S[t].c > 0 && S[z].c > 0 && (double(dens)/double(S[t].c*S[z].c)) >= dalpha ) {
                                lrs[t] = S[t].c; lrs[z] = S[z].c;
                                graph_convert();
                                csz[t] = lrs[t] - sz[t]; csz[z] = lrs[z] - sz[z];
                                branch_minCover(u,1,0);
                            }
                            else if (S[t].c > 0 && S[z].c > 0) {
                                subGraph_generate();
                                branchSubG(1,0);
                                // // branchNew(1, 0);

                                // lrs[t] = S[t].c; lrs[z] = S[z].c; Sc[t].c = 0; Sc[z].c = 0;
                                // branch_basic(1);
                            }
                        }

                        // Sc[t].r = Sc[t].c;

                        if (numbranches - old_numbranches > maxBranchduration || numbranches - old_numbranches >= 10000000) {
                            printf("it=%d, u=%d, sz[t]=%d, sz[z]=%d, ct=%d, cz=%d, maxBranchduration=%d, mSize=%d, density=%.3f\n", it, u, sz[t], sz[z], S[t].c, S[z].c, numbranches - old_numbranches, mSize, S[t].c*S[z].c<=0?0:double(dens)/double(S[t].c*S[z].c));
                        }
                        if (mSize1 < mSize) printf("it=%d, u=%d, sz[t]=%d, sz[z]=%d, ct=%d, cz=%d, maxBranchduration=%d, mSize=%d, density=%.3f\n", it, u, sz[t], sz[z], S[t].c, S[z].c, numbranches - old_numbranches, mSize, S[t].c*S[z].c<=0?0:double(dens)/double(S[t].c*S[z].c));
                        maxBranchduration = std::max(maxBranchduration, numbranches - old_numbranches);
                    }

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

        #ifdef _REVERSE_
            if (u == 0) break;
        #endif
    
    }

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
void MBCS::run(int alg) {

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

    #ifdef INIT_SIZE
    initSize();
    #endif

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
    uint32_t its = 0, _its = 0;
    for (int u = std::max(g->n[t]-sz[t], uint32_t(0)); u >= 0; u--) {
        if (g->cores[t][u] < g->maxcore/2) break;
        if (!visited[u]) {

            S[t].r = S[t].c = 0;
            S[z].r = S[z].c = 0;
            S[t].swapByPos(S[t].r++, S[t].pos(u)); S[t].c = S[t].r;

            uint32_t nextv = g->n[z]+1, mdv = 0, mdc = 0, du = g->p[t][u+1]- g->p[t][u];
            uint32_t cnts = 0;

            //test
            int qsz = 0, oldv = -1;
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
                    // S[z].swapByPos(S[z].c++, S[z].pos(v));
                    vqueue[z][qsz++] = v;
                    _deg[z][v] = dv;
                }
            }

            if (nextv >= g->n[z]) assert(mdv == 0);
            if (nextv >= g->n[z]) continue;
            
            // S[z].swapByPos(S[z].r++, S[z].pos(nextv));


            cnts = 0;
            for (uint32_t j = g->p[z][nextv]; j < g->p[z][nextv+1]; ++j){
                uint32_t u1 = g->e[z][j];
                if (u != u1) {
                    S[t].swapByPos(S[t].c++, S[t].pos(u1));
                    cnts++;
                }
            }

            cnts = 0;

            sort_deg(vqueue[z], z, 0, qsz);
            bool fg = true;
            for (uint32_t i = 0; i < qsz; ++i) {
                uint32_t v = vqueue[z][i];
                if (qsz * _deg[z][v] <= mSize) break;
                if (g->deg(z,v) <= S[t].cSize()) {
                    uint32_t stsz = S[t].r;
                    for (uint32_t j = g->p[z][v]; j < g->p[z][v+1]; ++j){
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
                        if (!g->is_nbr(z, v, u1)) {
                            S[t].swapByPos(--S[t].c, j); j--;
                        }

                    } 
                }
                cnts++;
                if (S[t].c < sz[t] || qsz * S[t].c <= mSize) break;
                // assert(cnts < 100);
                if ((i+1) >= sz[z] && mSize < (i+1) * S[t].c && S[t].c >= sz[t]){
                    mSize = (i+1) * S[t].c; 
                    realsz[t] = S[t].c; realsz[z] = (i+1);
                    if (fg) {
                        fg = false;
                        for (uint32_t j = 0; j < S[t].c; ++j)
                            visited[S[t][j]] = true;
                    }
                    _its = its;
                }
            }

            // while (S[z].cSize() > 0 && S[t].cSize() > 0 && S[t].c *S[z].c > mSize) {
            //     nextv = g->n[z]+1; mdv = 0; mdc = 0;
            
            //     for (int i = S[z].r; i < S[z].c; ++i) {
            //         uint32_t v = S[z][i];
            //         uint32_t dv = g->p[z][v+1] - g->p[z][v];

            //         if (dv >= sz[t] && mSize < dv * S[z].c) {
            //             // if (mdc < g->cores[z][v]) {
            //             //     nextv = v; mdc = g->cores[z][v]; mdv = dv;
            //             // }
            //             // else 
            //             if (dv >= mdv) {
            //                 nextv = v; mdv = dv;
            //             }
            //         }
            //         if (mSize >= dv * S[z].c) {
            //             S[z].swapByPos(--S[z].c, i); i--;
            //         }

            //     }
            //     if (nextv >= g->n[z]) assert(S[z].cSize() == 0);
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

            //     cnts++;
            //     // assert(cnts < 100);
            //     if (S[z].r >= sz[z] && mSize < S[z].r * S[t].c) {
            //         mSize = S[z].r * S[t].c;
            //         // for (uint32_t j = 0; j < S[t].c; ++j)
            //         //     visited[S[t][j]] = true;
            //     }
            // }

            // if (S[t].c >= sz[t] && S[z].c >= sz[t] && mSize < S[t].c *S[z].c)
            //     mSize = S[t].c *S[z].c;

            if (its >= 1800 || its > _its+600) break;
            its++;
            if (u == 0) break;
            u -= 2;
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

    printf("Init mSize=%d\n", mSize);
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Init time:" << duration.count() << "ms" << std::endl;
}


void MBCS::gemSize(uint32_t t, uint32_t tr, uint32_t zr, uint32_t ne) {

    if (tr >= zr) {
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
    if (tr >= sz[t] && zr >= sz[t^1] && mSize < tr * zr) {
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
    sdeg[0].resize(g->n1+1);
    sdeg[1].resize(g->n2+1);
    
    bin[0].resize(std::max(g->n[0], g->n[1])+1);
    bin[1].resize(std::max(g->n[0], g->n[1])+1);

    ws.resize(g->n1+g->n2+1);
    vqueue[0].resize(g->n1+1);
    vqueue[1].resize(g->n2+1);
    _deg[0].resize(g->n1+1);
    _deg[1].resize(g->n2+1);
    sz[0] = ls; sz[1] = rs;
    // exit(1);
}


MBCS::~MBCS()
{
    if (g != NULL) delete g;
    if (gs != NULL) delete gs;
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
  
    assert(S[1].r <= S[1].c);
    assert(S[0].r <= S[0].c);
    if (S[0].c < sz[0] || S[1].c < sz[1]) return;
    if (mSize >= S[0].c * S[1].c) return;

    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (S[0].cSize() == 0) {
        if (mSize < S[0].r * S[1].c) {
            gemSize(0, S[0].r, S[1].c, ne);
            // printf("mSize=%d\n",mSize);
            // mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);
            // printRes();
        }
        return;
    }
    else if (S[1].cSize() == 0) {
        if (mSize < S[0].c * S[1].r) {
            gemSize(0, S[0].c, S[1].r, ne);
            // printf("mSize=%d\n",mSize);
            // mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
            // assert(S[0].c * S[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].c; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].r; ++v)
            //     realRes[1].emplace_back(S[1][v]);
            // printRes();
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
            // mSize = (S[z].r+d) * (S[t].r+1); realsz[t] = (S[t].r+1); realsz[z] = (S[z].r+d);
            // assert((S[z].r+d) * (S[t].r+1) != 880);
            
            gemSize(t,S[t].r+1, S[z].r+d, ne);
        }
        else if (S[z].c * S[t].r > mSize &&  S[t].r >= sz[t] && S[z].c >= sz[z]) {
            // mSize = S[z].c * S[t].r;
            // mSize = S[z].c * S[t].r; realsz[t] = S[t].r; realsz[z] = S[z].c;
            // assert(S[z].c * S[t].r != 880);
            gemSize(t, S[t].r, S[z].c, ne);
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
                    else for (uint32_t j = g->p[ttt][u]; j < g->p[ttt][u+1]; ++j) {
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

        if (S[t].c < sz[t] || S[z].c < sz[z] || S[0].c * S[1].c <= mSize)  {
            S[t].c = c[t]; S[t].r = r[t]; 
            S[z].c = c[z]; S[z].r = r[z];
            return;
        }

        for(uint32_t j = S[t].r; j < S[t].c; j++) {
            uint32_t u = S[t][j];
            if (_deg[t][u] == S[z].c - r[z])  
                S[t].swapByPos(S[t].r++, j);
        }
        for(int j = S[z].r; j < S[z].c; j++) {
            uint32_t u = S[z][j];
            if (_deg[z][u] == S[t].c - r[t]) 
                S[z].swapByPos(S[z].r++, j);
        }

        maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};
        uint32_t minI = 0, minU = 0, mnume = 0, minE = std::max(S[t].c, S[z].c)+1;

        mnume = S[t].r*S[z].c+S[z].r*(S[t].c-S[t].r);
        if (S[t].cSize() > S[z].cSize()) {
            t = z; z = t^1;
        }

        qsz[t] = 0; qsz[z] = 0;
        for(uint32_t j = S[t].r; j < S[t].c; j++) {
            uint32_t e, u = S[t][j];
            _deg[t][u] -= (S[z].r - r[z]);
            e = _deg[t][u];
            vqueue[t][j] = S[z].c - S[z].r - e;

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
            uint32_t e = 0, u = S[z][j];
            _deg[z][u] -= (S[t].r - r[t]);
            e = _deg[z][u]; 
            vqueue[z][j] = S[t].c - S[t].r - e;

            if(e > maxE[z]) {
                maxE[z] = e;
                maxI[z] = j;
                maxUV[z] = u;
            }
        }

        if (S[t].cSize() > 0 && S[z].cSize() > 0){


            // for (uint32_t i = S[t].r; i < S[t].c; ++i) {
            //     uint32_t u = S[t][i];
            //     uint32_t d = 0;
            //     for (uint32_t j = S[z].r; j < S[z].c; ++j) {
            //         uint32_t v = S[z][j];
            //         if (g->is_nbr(t,u,v)) d++;
            //     }
            //     assert(d == _deg[t][u]);
            // }

            // for (uint32_t i = S[z].r; i < S[z].c; ++i) {
            //     uint32_t u = S[z][i];
            //     uint32_t d = 0;
            //     for (uint32_t j = S[t].r; j < S[t].c; ++j) {
            //         uint32_t v = S[t][j];
            //         if (g->is_nbr(z,u,v)) d++;
            //     }
            //     assert(d == _deg[z][u]);
            // }

            // std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c, [](int a, int b){return a > b;});
            // std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c, [](int a, int b){return a > b;});
            // if (vqueue[z].size() > S[z].c) {
            //     printf("t=%d, z=%d, vqueue[z].size=%d, S[z].c=%d\n", t,z, vqueue[z].size(), S[z].c);
            // }
 
            std::sort(vqueue[t].begin()+S[t].r, vqueue[t].begin()+S[t].c);
            std::sort(vqueue[z].begin()+S[z].r, vqueue[z].begin()+S[z].c);

            uint32_t mxnbz = 0;
            for (uint32_t i = S[z].r; i < S[z].c; ++i) {
                uint32_t nbz = vqueue[z][i];
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
        if(e > maxE[z]) {
            maxE[z] = e;
            maxI[z] = i;
            maxUV[z] = u;
        }
    }
#endif

    if (S[0].cSize() == 0) {

        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            gemSize(0, S[0].r, S[1].c, ne);
            // mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
            // assert(S[0].r * S[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < S[0].r; ++v)
            //     realRes[0].emplace_back(S[0][v]);
            // for (uint32_t v = 0; v < S[1].c; ++v)
            //     realRes[1].emplace_back(S[1][v]);

            for (uint32_t v = 0; v < S[0].r; ++v){
                for (uint32_t u = 0; u < S[1].c; ++u)
                    assert(g->is_nbr(0, S[0][v], S[1][u]));
            }
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
            gemSize(0, S[0].c, S[1].r, ne);
            // mSize = S[0].c * S[1].r; realsz[0] = S[0].c; realsz[1] = S[1].r;
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
        assert(maxE[z] == 0);
        if (mSize < S[0].r * S[1].c && S[0].r >= sz[0] && S[1].c >= sz[1]) {
            gemSize(0, S[0].r, S[1].c, ne);
            // mSize = S[0].r * S[1].c; realsz[0] = S[0].r; realsz[1] = S[1].c;
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
        uint32_t cc = S[z].r;
        for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
            uint32_t v = g->e[t][j];
            uint32_t pv = S[z].pos(v);
            if(S[z].r <= pv && pv < S[z].c) {
                // if ((_deg[z][v]+r[t]) >= mSize / (_deg[z][v]+r[t]))
                    S[z].swapByPos(cc++, pv);
            }
        }
        wsSize = S[z].c - cc;
        if (wsSize > ws[deep].capacity()) ws[deep].resize(wsSize+1);
        memcpy(ws[deep].data(), S[z].begin() + cc, sizeof(uint32_t) * wsSize);
    }

    maxsubbranches = wsSize+1 > maxsubbranches ? wsSize+1 : maxsubbranches;
    fflush(stdout); 
    uint32_t ct = S[t].c, cz = S[z].c;
    uint32_t dv = _deg[z][ws[deep][0]] + r[t];

    // if (wsSize == 1 && dv + 1 == S[t].c) {
    //     uint32_t v = ws[deep][0];
    //     S[t].swapByPos(S[t].r++, S[t].pos(u));
    //     S[z].swapByPos(S[z].r++, S[z].pos(v));
    //     branchNew(deep, ne+1);
    //     S[t].r = r[t];
    //     S[t].c = c[t];
    //     S[z].r = r[z];
    //     S[z].c = c[z];
    //     return;
    // }
 
    // std::sort(wsSize[deep].begin(), wsSize[deep].begin()+wsSize, [](int a, int b){return a > b;});
    std::sort(ws[deep].begin(), ws[deep].begin()+wsSize, [&](int a, int b){return _deg[z][a] < _deg[z][b];});

    uint32_t wsSize1 = wsSize;
    if (wsSize > 1) {
        uint32_t x = 1;
        while (x+1 < wsSize) {
            uint32_t vl = ws[deep][wsSize-x];
            uint32_t qsz = 0;
            if(g->deg(z, vl) > 2*S[t].cSize()) {
                for(uint32_t i = S[t].r; i < S[t].c; i++) {
                    uint32_t w = S[t][i];
                    if(!g->is_nbr(z, vl, w)) {
                        vqueue[t][qsz++] = w;
                    }
                }
            }
            else {
                uint32_t cc = S[t].r;
                for(uint32_t j = g->p[z][vl]; j < g->p[z][vl + 1]; j++) {
                    uint32_t v = g->e[z][j];
                    uint32_t pv = S[t].pos(v);
                    if(S[t].r <= pv && pv < S[t].c) {
                        S[t].swapByPos(cc++, pv);
                    }
                }
                memcpy(vqueue[t].data(), S[t].begin() + cc, sizeof(uint32_t) * (S[t].c - cc));
                qsz = (S[t].c - cc);
            }
            wsSize1 = wsSize; wsSize = 0;   
            for (uint32_t i = 0; i+x < wsSize1; ++i) {
                uint32_t v = ws[deep][i];
                uint32_t nbrs = 0;
                for (uint32_t j = 0; j < qsz; ++j) {
                    uint32_t u1 = vqueue[t][j];
                    if (v >= g->n[z]) {
                        printf("u1=%d, x=%d, wsSize1=%d, qsz=%d, u=%d, deep=%d, v=%d, t=%d\n", u1, x, wsSize1, qsz, u, deep, v, t);
                    }
                    assert(u1 < g->n[t]);
                    assert(v < g->n[z]);
                    assert(v > 0 );
                    if (g->is_nbr(z, v, u1)) {
                        nbrs++;
                        if (nbrs > 1) break;
                    }
                }
                if (nbrs < 1) subbcnt++;
                if (nbrs == 0) continue;
                else ws[deep][wsSize++] = v;
            }
            for (uint32_t i = wsSize1-x; i < wsSize1; ++i) ws[deep][wsSize++] = ws[deep][i];
            x++; break;
        }
    }

    uint32_t ii = 0; 
    uint32_t rt = S[t].r, rz = S[z].r;
    for (uint32_t i = 0; i < wsSize; ++i) {
        uint32_t v = ws[deep][i];

        S[z].swapByPos(S[z].r++, S[z].pos(v));
        if(g->deg(z, v) > 2*S[t].cSize()) {
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

        uint32_t szr = S[z].r;
        branchNew(deep+1, ne);
        S[t].c = ct;
        assert(S[z].c == cz - ii); ii++;
        if (S[z][S[z].r-1] != v) {
            S[0].print();
            S[1].print();
            printf("z=%d, v=%d, S[z][S[z].r-1]=%d, szr=%d, S[z].r=%d, S[t].r=%d, S[t].c=%d\n", z, v, S[z][S[z].r-1], szr, S[z].r, S[t].r, S[t].c);
        }
        assert(S[z][S[z].r-1] == v);
        S[z].swapByPos(--S[z].r,--S[z].c);
    }

    // if (wsSize == 1) {
    //     uint32_t x = S[z].c - 1;
    //     uint32_t y = mSize / x;

    //     if (y >= (x+1)*(S[z].c - dv)) {
    //         S[t].r = r[t];
    //         S[t].c = c[t];
    //         S[z].r = r[z];
    //         S[z].c = c[z];
    //         return;
    //     }
    // }
    // pivot vertex
    S[t].swapByPos(S[t].r++, S[t].pos(u));

    branchNew(deep+1, ne);

    S[t].r = r[t];
    S[t].c = c[t];
    S[z].r = r[z];
    S[z].c = c[z];

}

void MBCS::branchSubG(uint32_t deep, uint32_t ne) {
    if (Sc[0].c < sz[0] || Sc[1].c < sz[1]) return;
    if (mSize >= Sc[0].c * Sc[1].c) return;
    assert(Sc[1].r <= Sc[1].c);
    if (Sc[0].r > Sc[0].c) printf("Sc[0].r=%d, Sc[0].c=%d, deep=%d, mSize=%d\n", Sc[0].r, Sc[0].c, deep, mSize);
    assert(Sc[0].r <= Sc[0].c);

    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (Sc[0].cSize() == 0) {
        if (mSize < Sc[0].r * Sc[1].c) {
            gemSize(0, Sc[0].r, Sc[1].c, ne);
            // printf("mSize=%d\n",mSize);
            // mSize = Sc[0].r * Sc[1].c; realsz[0] = Sc[0].r; realsz[1] = Sc[1].c;
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].r; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].c; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // printRes();
        }
        return;
    }
    else if (Sc[1].cSize() == 0) {
        if (mSize < Sc[0].c * Sc[1].r) {
            gemSize(0, Sc[0].c, Sc[1].r, ne);
            // printf("mSize=%d\n",mSize);
            // mSize = Sc[0].c * Sc[1].r; realsz[0] = Sc[0].c; realsz[1] = Sc[1].r;
            // assert(Sc[0].c * Sc[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].c; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].r; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // printRes();
        }
    
        return;
    }

    uint32_t r[2], c[2];
    r[0] = Sc[0].r; r[1] = Sc[1].r;
    c[0] = Sc[0].c; c[1] = Sc[1].c;

    #ifdef REDUCTION_
    for (uint32_t t = 0; t <= 1; ++t) 
    if (Sc[t].cSize() == 1) {
        uint32_t d = 0;
        uint32_t z = t^1;
        uint32_t u = Sc[t][Sc[t].r];
        for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
            uint32_t v = Sc[z][i];
            if (is_nbr_sub(t, u, v)) d++;
        }
        if ((Sc[z].r+d) * (Sc[t].r+1) > mSize  &&  Sc[t].r+1 >= sz[t] && Sc[z].r+d >= sz[z]) {
            // mSize = (Sc[z].r+d) * (Sc[t].r+1);
            // mSize = (Sc[z].r+d) * (Sc[t].r+1); realsz[t] = (Sc[t].r+1); realsz[z] = (Sc[z].r+d);
            // assert((Sc[z].r+d) * (Sc[t].r+1) != 880);
            
            gemSize(t,Sc[t].r+1, Sc[z].r+d, ne);
        }
        else if (Sc[z].c * Sc[t].r > mSize &&  Sc[t].r >= sz[t] && Sc[z].c >= sz[z]) {
            // mSize = Sc[z].c * Sc[t].r;
            // mSize = Sc[z].c * Sc[t].r; realsz[t] = Sc[t].r; realsz[z] = Sc[z].c;
            // assert(Sc[z].c * Sc[t].r != 880);
            gemSize(t, Sc[t].r, Sc[z].c, ne);
        }
        return;
    }
    #endif

    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};

    uint32_t t = 0, z = 1;
    if (Sc[1].cSize() < Sc[0].cSize()) {t = 1; z = 0;} 

    #ifdef REDUCTION_
    if (Sc[t].c == sz[t]) {

        for (int i = Sc[z].r; i < Sc[z].c; ++i) {
            uint32_t u = Sc[z][i];
            for (uint32_t j = Sc[t].r; j < Sc[t].c; ++j) {
                uint32_t v = Sc[t][j];
                if (!is_nbr_sub(z,u,v)) {
                    Sc[z].swapByPos(i--, --Sc[z].c);
                    break;
                }

            }
            if (Sc[t].c * Sc[z].c <= mSize) break;
        }

        Sc[t].r = Sc[t].c;
        branchSubG(deep+1, ne);
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].c = c[z];
        return;
    }
    if (Sc[z].c == sz[z]) {

        for (int i = Sc[t].r; i < Sc[t].c; ++i) {
            uint32_t u = Sc[t][i];
            for (uint32_t j = Sc[z].r; j < Sc[z].c; ++j) {
                uint32_t v = Sc[z][j];
                if (!is_nbr_sub(t,u,v)) {
                    Sc[t].swapByPos(i--, --Sc[t].c);
                    break;
                }

            }
            if (Sc[t].c * Sc[z].c <= mSize) break;
        }

        Sc[z].r = Sc[z].c;
        branchSubG(deep+1, ne);
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        Sc[t].c = c[t];
        return;
    }
    #endif

    uint32_t qsz[2] = {0};
    uint32_t qr[2] = {0};
    uint32_t minI = 0, minU = 0, mnume = 0, minE = std::max(Sc[t].c, Sc[z].c)+1;

        for(uint32_t j = Sc[z].r; j < c[z]; j++) _deg[z][Sc[z][j]] = 0;
        for(int i = Sc[t].r; i < Sc[t].c; i++) {//scan all vertices in C
        // for(uint32_t i = Sc[t].c; i <= Sc[t].c; i++) {
            uint32_t u = Sc[t][i];
            uint32_t e = 0;
            if(gs->deg(t, u) > 2*Sc[z].cSize()) {
                for(uint32_t j = Sc[z].r; j < Sc[z].c; j++) {
                    uint32_t v = Sc[z][j];
                    if(is_nbr_sub(t, u, v)) {
                        e++; _deg[z][v]++;
                    }
                }
            }
            else {
                for(uint32_t j = gs->p[t][u]; j < gs->p[t][u + 1]; j++) {
                    uint32_t v = gs->e[t][j];
                    uint32_t pv = Sc[z].pos(v);
                    if(Sc[z].r <= pv && pv < Sc[z].c) {//in C
                        e++; _deg[z][v]++;
                    }
                }
            }
            
            _deg[t][u] = e;
            #ifdef UB1_
            if (Sc[z].r+e < sz[z] || Sc[t].c * (Sc[z].r+e) <= mSize) {
                Sc[t].swapByPos(--Sc[t].c, i--);
                vqueue[t][qsz[t]++] = u;
                if (Sc[t].c < sz[t] || Sc[t].c * Sc[z].c <= mSize)  {
                    Sc[t].r = r[t];
                    Sc[t].c = c[t];
                    Sc[z].r = r[z];
                    Sc[z].c = c[z];
                    return;
                }
            }
            #endif
        }

        #ifdef UB1_
        for(int j = Sc[z].r; j < Sc[z].c; j++) {
            uint32_t u = Sc[z][j];
            uint32_t e = _deg[z][u];
            if (r[t]+e < sz[t] || Sc[z].c *(r[t]+e) < mSize) {
                Sc[z].swapByPos(--Sc[z].c, j--);
                vqueue[z][qsz[z]++] = u;
            }
        }
        #endif

       // if (deep == 1 && gs->realLable[1][Sc[1][0]] == 7141 && sz[1] == 13) {
       //      printf("t=%d, z=%d, qsz[t]=%d, qsz[z]=%d\n", t, z, qsz[t], qsz[z]);
       //  }
#ifdef UB1_
        while(qr[0] < qsz[0] || qr[1] < qsz[1]) {
            for (uint32_t ttt = 0; ttt <= 1; ++ttt) {
                uint32_t zzz = ttt^1;
                for (uint32_t i = qr[ttt]; i < qsz[ttt]; ++i) {
                    uint32_t u = vqueue[ttt][i];
                    _deg[ttt][u] = 0;

                    if(gs->deg(ttt, u) > 5*Sc[zzz].cSize()) {
                        for(int j = Sc[zzz].r; j < Sc[zzz].c; j++) {
                            uint32_t v = Sc[zzz][j];
                            if(is_nbr_sub(ttt, u, v)) {
                                uint32_t d = --_deg[zzz][v];
                                d += Sc[ttt].r;
                                if (d < sz[ttt] || d * Sc[zzz].c <= mSize) {
                                    Sc[zzz].swapByPos(--Sc[zzz].c, j--);
                                    vqueue[zzz][qsz[zzz]++] = v;
                                    if (Sc[t].c < sz[t] || Sc[z].c < sz[z] || Sc[0].c * Sc[1].c <= mSize)  break;
                                }
                            }
                        }
                    }
                    else for (uint32_t j = gs->p[ttt][u]; j < gs->p[ttt][u+1]; ++j) {
                        uint32_t v = gs->e[ttt][j];
                        uint32_t pv = Sc[zzz].pos(v);
                        if (pv < Sc[zzz].c && Sc[zzz].r <= pv) {
                            uint32_t d = --_deg[zzz][v];
                            d += Sc[ttt].r;
                            if (d < sz[ttt] || d * Sc[zzz].c <= mSize) {
                                Sc[zzz].swapByPos(--Sc[zzz].c, pv);
                                vqueue[zzz][qsz[zzz]++] = v;
                                if (Sc[t].c < sz[t] || Sc[z].c < sz[z] || Sc[0].c * Sc[1].c <= mSize)  break;
                            }
                        }
                    }
                    if (Sc[t].c < sz[t] || Sc[z].c < sz[z] || Sc[0].c * Sc[1].c <= mSize)  {
                        i = qsz[ttt]; qr[zzz] = qsz[zzz];
                        break;
                    }
                }
                qr[ttt] = qsz[ttt];

            }
        }

        if (Sc[t].c < sz[t] || Sc[z].c < sz[z] || Sc[0].c * Sc[1].c <= mSize)  {
            Sc[t].c = c[t]; Sc[t].r = r[t]; 
            Sc[z].c = c[z]; Sc[z].r = r[z];
            return;
        }
#endif
        // if (deep == 7 && gs->realLable[0][Sc[0][0]] == 36360){
        //     printf("deep1=%d, _deg[z][20]=%d\n", deep, _deg[z][20]);
        // }

    #ifdef REDUCTION_
    // #ifdef UB_HNBR_
        for(uint32_t j = Sc[t].r; j < Sc[t].c; j++) {
            uint32_t u = Sc[t][j];
            if (_deg[t][u] == Sc[z].c - r[z]) 
                Sc[t].swapByPos(Sc[t].r++, j);
        }
        for(int j = Sc[z].r; j < Sc[z].c; j++) {
            uint32_t u = Sc[z][j];
            if (_deg[z][u] == Sc[t].c - r[t]) 
                Sc[z].swapByPos(Sc[z].r++, j);
        }
    // #endif
    #endif

        //  if (deep == 7 && gs->realLable[0][Sc[0][0]] == 36360){
        //     printf("deep2=%d, _deg[z][20]=%d, bool=%d\n", deep,_deg[z][20], Sc[t].cSize() > 0 && Sc[z].cSize() > 0);
        // }

        if (Sc[t].cSize() > Sc[z].cSize()) {
            t = z; z = t^1;
        }

        // if (deep == 7 && gs->realLable[0][Sc[0][0]] == 36360){
        //     printf("deep2=%d, _deg[z][20]=%d, bool=%d\n", deep,_deg[z][20], Sc[t].cSize() > 0 && Sc[z].cSize() > 0);
        // }
        // if (deep == 10) printf("1111 Sc[0].r=%d, Sc[0].c=%d, deep=%d, mSize=%d\n", Sc[0].r, Sc[0].c, deep, mSize);
        maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {0};

        mnume = Sc[t].r*Sc[z].c+Sc[z].r*(Sc[t].c-Sc[t].r);
        qsz[t] = 0; qsz[z] = 0;
        for(uint32_t j = Sc[t].r; j < Sc[t].c; j++) {
            uint32_t e = 0, u = Sc[t][j];
            _deg[t][u] -= (Sc[z].r-r[z]);
         
            e = _deg[t][u];
            vqueue[t][j] = Sc[z].c - Sc[z].r - e;

            if(e > maxE[t]) {
                maxE[t] = e; maxI[t] = j; maxUV[t] = u;
            }
            if (e < minE) {
                minE = e; minI = j; minU = u;
            }
            mnume += e;
        }
        #ifdef REDUCTION_
        // #ifndef UB_HNBR_
            if (maxE[t] == Sc[z].cSize() && Sc[t].r < Sc[t].c) {
                uint32_t u = maxUV[t];
                Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u));
                branchSubG(deep+1, ne);
                Sc[t].r = r[t];
                Sc[t].c = c[t];
                Sc[z].r = r[z];
                Sc[z].c = c[z];
                return;
            }
        // #endif
        #endif

        for(int j = Sc[z].r; j < Sc[z].c; j++) {
            uint32_t e = 0, u = Sc[z][j];
            _deg[z][u] -= (Sc[t].r-r[t]);
            e = _deg[z][u]; 
            vqueue[z][j] = Sc[t].c - Sc[t].r - e;

            if(e > maxE[z]) {
                maxE[z] = e; maxI[z] = j; maxUV[z] = u;
            }
        }

        #ifdef REDUCTION_
        // #ifndef UB_HNBR_
            if (maxE[z] == Sc[t].cSize() && Sc[z].r < Sc[z].c) {
                uint32_t u = maxUV[z];
                Sc[z].swapByPos(Sc[z].r++, Sc[z].pos(u));
                branchSubG(deep+1, ne);
                Sc[t].r = r[t];
                Sc[t].c = c[t];
                Sc[z].r = r[z];
                Sc[z].c = c[z];
                return;
            } 
        // #endif
        #endif

        #ifdef UB_HNBR_
        if (Sc[t].cSize() > 0 && Sc[z].cSize() > 0){

        //     if (deep == 7 && gs->realLable[0][Sc[0][0]] == 36360){
        //     printf("deep3=%d, _deg[z][20]=%d, dertr[0]=%d, dertr[1]=%d\n", deep,_deg[z][20], Sc[0].r - r[0], Sc[1].r- r[1]);
        // }
            // for (uint32_t i = Sc[t].r; i < Sc[t].c; ++i) {
            //     uint32_t u = Sc[t][i];
            //     uint32_t d = 0;
            //     for (uint32_t j = Sc[z].r; j < Sc[z].c; ++j) {
            //         uint32_t v = Sc[z][j];
            //         if (is_nbr_sub(t,u,v)) d++;
            //     }
            //     if (d != _deg[t][u]) {
            //         printf("t=%d, z=%d, i=%d, u=%d, d=%d, _deg[z][u]=%d\n", t,z, i,Sc[t][i],d, _deg[t][u]);
            //         printf("deep=%d, Sc[0][0]=%d, Sc[0].r=%d, Sc[0].c=%d, Sc[1].r=%d, Sc[1].c=%d\n", deep, gs->realLable[0][Sc[0][0]], Sc[0].r, Sc[0].c, Sc[1].r, Sc[1].c);
            //     }
            //     assert(d == _deg[t][u]);
            // }

            // for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
            //     uint32_t u = Sc[z][i];
            //     uint32_t d = 0;
            //     for (uint32_t j = Sc[t].r; j < Sc[t].c; ++j) {
            //         uint32_t v = Sc[t][j];
            //         if (is_nbr_sub(z,u,v)) d++;
            //     }
            //     if (d != _deg[z][u]) {
            //         printf("t=%d, z=%d, i=%d, u=%d, d=%d, _deg[z][u]=%d\n", t,z, i,Sc[z][i],d, _deg[z][u]);
            //         printf("deep=%d, Sc[0][0]=%d, Sc[0].r=%d, Sc[0].c=%d, Sc[1].r=%d, Sc[1].c=%d\n", deep, gs->realLable[0][Sc[0][0]], Sc[0].r, Sc[0].c, Sc[1].r, Sc[1].c);
            //     }
            //     assert(d == _deg[z][u]);
            // }

            // std::sort(vqueue[t].begin()+Sc[t].r, vqueue[t].begin()+Sc[t].c, [](int a, int b){return a > b;});
            // std::sort(vqueue[z].begin()+Sc[z].r, vqueue[z].begin()+Sc[z].c, [](int a, int b){return a > b;});
            // if (vqueue[z].size() > Sc[z].c) {
            //     printf("t=%d, z=%d, vqueue[z].size=%d, Sc[z].c=%d\n", t,z, vqueue[z].size(), Sc[z].c);
            // }
 
            // std::sort(vqueue[t].begin()+Sc[t].r, vqueue[t].begin()+Sc[t].c);
            // std::sort(vqueue[z].begin()+Sc[z].r, vqueue[z].begin()+Sc[z].c);
            sort_in(vqueue[t], t, Sc[t].r, Sc[t].c);
            sort_in(vqueue[z], z, Sc[z].r, Sc[z].c);

            uint32_t mxnbz = 0;
            for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
                uint32_t nbz = vqueue[z][i];
                mxnbz = std::max(mxnbz, nbz);
            }

            uint32_t ubm = 0, zsz = Sc[z].c, nonedges[2] = {0};
            for (uint32_t i = Sc[t].r; i < Sc[t].c; ++i) {
                uint32_t _d = vqueue[t][i];
                uint32_t d = Sc[z].c - _d;
                nonedges[t] += _d;

                while (zsz > Sc[z].r && nonedges[t] > nonedges[z]) {
                    zsz--;
                    uint32_t _dz = vqueue[z][zsz];
                    nonedges[z] += _dz;
                    if (zsz == Sc[z].r) break;
                }
                while (zsz > Sc[z].r && vqueue[z][zsz-1] + i + 1 > Sc[t].c) {
                    nonedges[z] += vqueue[z][--zsz];
                    if (zsz == Sc[z].r) break;
                }

                if (zsz > Sc[z].r && i-Sc[t].r >= mxnbz) {
                    uint32_t j = i - mxnbz;
                    uint32_t tail = vqueue[t][j] + vqueue[t][j+1];
                    while (zsz > Sc[z].r && tail+zsz > Sc[z].c)
                        nonedges[z] += vqueue[z][--zsz];
                }

                uint32_t ubz = std::min(d, zsz);
                if (ubz < sz[z]) break;
                if (i+1 >= sz[t]) {
                    ubm = std::max(ubm, (i+1) * ubz);
                }
            }
        
            if (ubm <= mSize) {
                Sc[t].c = c[t]; Sc[t].r = r[t]; 
                Sc[z].c = c[z]; Sc[z].r = r[z];
                return;
            }
        }
        #endif

    #ifdef REDUCTION_
    if (Sc[0].cSize() == 0) {
        assert(Sc[0].cSize() == 0);
        if (mSize < Sc[0].r * Sc[1].c && Sc[0].r >= sz[0] && Sc[1].c >= sz[1]) {
            gemSize(0, Sc[0].r, Sc[1].c, ne);
            // mSize = Sc[0].r * Sc[1].c; realsz[0] = Sc[0].r; realsz[1] = Sc[1].c;
            // assert(Sc[0].r * Sc[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].r; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].c; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // if (mSize == Sc[0].r * Sc[1].c)
            // for (uint32_t v = 0; v < Sc[0].r; ++v){
            //     for (uint32_t u = 0; u < Sc[1].c; ++u){
            //         if (!is_nbr_sub(0, Sc[0][v], Sc[1][u])) {
            //             printf("mSize=%d, i=%d, j=%d, Sc[0][v]=%d, Sc[1][u]=%d, ux=%d, Sc[0][0]=%d\n", mSize, v,u, Sc[0][v], Sc[1][u], gs->realLable[0][Sc[0][0]], Sc[0][0]);
            //             printf("mSize=%d, sz[t]=%d, sz[z]=%d\n", mSize, sz[t],sz[z]);
            //         }
            //         assert(is_nbr_sub(0, Sc[0][v], Sc[1][u]));
            //     }
            // }
            // printRes();
        }
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        return;
    }
    else if (Sc[1].cSize() == 0) {

        if (mSize < Sc[0].c * Sc[1].r && Sc[0].c >= sz[0] && Sc[1].r >= sz[1]) {
            gemSize(0, Sc[0].c, Sc[1].r, ne);
            // mSize = Sc[0].c * Sc[1].r; realsz[0] = Sc[0].c; realsz[1] = Sc[1].r;
            // assert(Sc[0].c * Sc[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].c; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].r; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // for (uint32_t v = 0; v < Sc[0].r; ++v){
            //     for (uint32_t u = 0; u < Sc[1].c; ++u)
            //         assert(is_nbr_sub(0, Sc[0][v], Sc[1][u]));
            // }
            // printRes();
        }
  
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        return;
    }
    
    // printf("maxE[0]=%d, maxE[1]=%d\n", maxE[0],maxE[1]);
    // exit(1);
    // #ifdef REDUCTION_
    if (maxE[z] == 0) assert(maxE[t] == 0);
    if (maxE[t] == 0) {
        assert(maxE[z] == 0);
        if (mSize < Sc[0].r * Sc[1].c && Sc[0].r >= sz[0] && Sc[1].c >= sz[1]) {
            gemSize(0, Sc[0].r, Sc[1].c, ne);
            // mSize = Sc[0].r * Sc[1].c; realsz[0] = Sc[0].r; realsz[1] = Sc[1].c;
            // assert(Sc[0].r * Sc[1].c != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].r; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].c; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // printRes();
        }

        if (mSize < Sc[0].c * Sc[1].r && Sc[0].c >= sz[0] && Sc[1].r >= sz[1]) {
            gemSize(0, Sc[0].c, Sc[1].r, ne);
            // mSize = Sc[0].c * Sc[1].r; realsz[0] = Sc[0].c; realsz[1] = Sc[1].r;
            // assert(Sc[0].c * Sc[1].r != 880);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].c; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].r; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
            // printRes();
        }
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        return;
    }
    #endif

    if (Sc[0].r >= sz[0] && Sc[0].r * Sc[1].c > mSize)
        gemSize(0, Sc[0].r, Sc[1].c, ne);
    else if (Sc[1].r >= sz[1] && Sc[0].c * Sc[1].r > mSize)
        gemSize(0, Sc[0].c, Sc[1].r, ne);

    // assert(Sc[0].cSize() > 0);
    // assert(Sc[1].cSize() > 0);

    // if (deep == 1 && gs->realLable[1][Sc[1][0]] == 7141 && sz[1] == 13) {
    //     printf("mSize=%d, sz[t]=%d, sz[z]=%d, szt=%d, szz=%d, minU=%d, maxUV[t]=%d, maxUV[z]=%d\n", mSize, sz[t], sz[z], Sc[0].c, Sc[1].c, minU, maxUV[t], maxUV[z]);
    //     printf("mSize=%d, Sc[0].r=%d, Sc[0].c=%d, Sc[1].r=%d, Sc[1].c=%d, maxE[0]=%d, maxE[1]=%d, BIBRANCH=%d\n", mSize, Sc[0].r, Sc[0].c, Sc[1].r, Sc[1].c, maxE[0], maxE[1], BIBRANCH);

    // //     for (uint32_t i = 0; i < Sc[0].c; ++i) {
    // //         uint32_t u = Sc[0][i];
    // //         if (u == 27) {
    // //             printf("u=%d, is_nbr(27, 0), Sc[1][0]=%d\n", u, is_nbr_sub(0, u, Sc[1][0]));

    // //         }
    // //     }

    // //     // exit(1);
    // }

    // if (deep == 10) printf("222 Sc[0].r=%d, Sc[0].c=%d, deep=%d, mSize=%d\n", Sc[0].r, Sc[0].c, deep, mSize);
    if (BIBRANCH && minE + 4 < Sc[z].cSize() && Sc[t].cSize() * 20 < Sc[z].cSize()) {

        uint32_t u = minU, cz = Sc[z].c;
        // assert( Sc[t].pos(u) < Sc[t].c);
        Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u));

        if(gs->deg(t, u) > 2*Sc[z].cSize()) {
            for(int j = Sc[z].r; j < Sc[z].c; j++) {
                if(!is_nbr_sub(t, u, Sc[z][j])) {
                    Sc[z].swapByPos(--Sc[z].c, j); j--;
                }
            }
        }
        else {
            uint32_t cc = Sc[z].r;
            for(uint32_t j = gs->p[t][u]; j < gs->p[t][u + 1]; j++) {
                uint32_t w = gs->e[t][j];
                uint32_t pw = Sc[z].pos(w);
                if(Sc[z].r <= pw && pw < Sc[z].c) {
                    Sc[z].swapByPos(cc++, pw);
                }
            }
            Sc[z].c = cc;
        }

        branchSubG(deep+1, ne);
        Sc[z].c = cz;
        // assert(Sc[t][Sc[t].r-1] == u);
        Sc[t].swapByPos(--Sc[t].r, --Sc[t].c);

        branchSubG(deep+1, ne);
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        return;
    }

    // if (maxE[t]+1 >= Sc[z].cSize()) {
    //     if (maxE[t] == Sc[z].cSize()) {
    //         Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(maxUV[t]));
    //         for(uint32_t i = Sc[z].r; i < Sc[z].c; i++) 
    //             subdeg[Sc[z][i]] = 0;
    //         branchSubG(deep+1, ne);
    //         Sc[t].r--;
    //     }
    //     else {
            
    //         for(uint32_t i = Sc[z].r; i < Sc[z].c; i++) 
    //             subdeg[Sc[z][i]] = 0;
    //     }
    //     return;
    // }
    
    // for(uint32_t i = Sc[z].r; i < Sc[z].c; i++) {
    //     uint32_t u = Sc[z][i];
    //     uint32_t e = subdeg[u];
    //     subdeg[u] = 0;

    //     if(e > maxE[z]) {
    //         maxE[z] = e;
    //         maxI[z] = i;
    //         maxUV[z] = u;
    //     }
    // }

    // t = 0, z = 1;
    if (Sc[z].cSize() - maxE[t] >= Sc[t].cSize()-maxE[z]) {
        t = z; z = t^1;
    }
    // if (Sc[z].c - r[z] - maxE[t] >= Sc[t].c - r[t]-maxE[z]) {
    //     t = z; z = t^1;
    // }
    uint32_t u = maxUV[t], wsSize = 0;

    if(gs->deg(t, u) > 2*Sc[z].cSize()) {
        if (Sc[z].cSize() > ws[deep].capacity()) ws[deep].resize(Sc[z].cSize()+1);
        for(uint32_t i = Sc[z].r; i < Sc[z].c; i++) {
            uint32_t w = Sc[z][i];
            if(!is_nbr_sub(t, u, w)) {
                // if ((_deg[z][w]+r[t]) >= mSize / (_deg[z][w]+r[t]))
                    ws[deep][wsSize++] = Sc[z][i];
            }
        }
    }
    else {
        uint32_t cc = Sc[z].r;
        for(uint32_t j = gs->p[t][u]; j < gs->p[t][u + 1]; j++) {
            uint32_t v = gs->e[t][j];
            uint32_t pv = Sc[z].pos(v);
            if(Sc[z].r <= pv && pv < Sc[z].c) {
                // if ((_deg[z][v]+r[t]) >= mSize / (_deg[z][v]+r[t]))
                    Sc[z].swapByPos(cc++, pv);
            }
        }
        wsSize = Sc[z].c - cc;
        if (wsSize > ws[deep].capacity()) ws[deep].resize(wsSize+1);
        memcpy(ws[deep].data(), Sc[z].begin() + cc, sizeof(uint32_t) * wsSize);
    }

    // if (deep == 1 && gs->realLable[1][Sc[1][0]] == 7141 && sz[1] == 13) {
    //     for (uint32_t i = 0; i < wsSize; ++i)
    //         printf(" %d", ws[deep][i]);
    //     printf("\n");
    //     if (gs->deg(t, u) <= 2*Sc[z].cSize()) {
    //     for (uint32_t i = Sc[z].c - wsSize; i < Sc[z].c; ++i)
    //             printf(" %d", Sc[z][i]);
    //         printf("\n");
    //     }
    //     else {
    //         for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
    //             if(!is_nbr_sub(t, u, Sc[z][i]))
    //                 printf(" %d", Sc[z][i]);
    //         }
    //         printf("\n");
    //     }
    //     printf("mSize=%d, t=%d, z=%d, u=%d, b=%d\n", mSize, t, z, u, gs->deg(t, u) > 2*Sc[z].cSize());

    // //     exit(1);
    // }
    // if (deep == 10) printf("3333333 Sc[0].r=%d, Sc[0].c=%d, deep=%d, mSize=%d\n", Sc[0].r, Sc[0].c, deep, mSize);
    maxsubbranches = wsSize+1 > maxsubbranches ? wsSize+1 : maxsubbranches;
    uint32_t ct = Sc[t].c, cz = Sc[z].c;
    uint32_t dv = wsSize > 0 ? _deg[z][ws[deep][0]] + Sc[t].r : 0;

    #ifdef REDUCTION_
    if (wsSize == 1 && dv + 1 == Sc[t].c) {
        uint32_t v = ws[deep][0];
        Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u));
        Sc[z].swapByPos(Sc[z].r++, Sc[z].pos(v));
        branchSubG(deep, ne+1);
        Sc[t].r = r[t];
        Sc[t].c = c[t];
        Sc[z].r = r[z];
        Sc[z].c = c[z];
        return;
    }
    #endif
 
    // std::sort(wsSize[deep].begin(), wsSize[deep].begin()+wsSize, [](int a, int b){return a > b;});
    // std::sort(ws[deep].begin(), ws[deep].begin()+wsSize, [&](int a, int b){return _deg[z][a] < _deg[z][b];});

    uint32_t front = 0, tail = 0;
    if (wsSize > 1) {
        std::sort(ws[deep].begin(), ws[deep].begin()+wsSize, [&](int a, int b){return _deg[z][a] < _deg[z][b];});
        uint32_t qcsz = 0, nw = 0, v = ws[deep][wsSize-1];
        uint32_t dmx = 0;
        // for (uint32_t j = 0; j < wsSize; ++j) {
        //     uint32_t w = ws[deep][j];
        //     uint32_t dw = _deg[z][w];
        //     if (dmx < dw) {
        //         v = w; dmx = dw;
        //     }
        // }

        if(gs->deg(z,v) > 2*Sc[t].cSize()) {
            for(uint32_t j = Sc[t].r; j < Sc[t].c; j++) {
                uint32_t w = Sc[t][j];
                if(!is_nbr_sub(z, v, w)) {
                    vqueue[t][qcsz++] = w;
                }
            }
        }
        else {
            uint32_t cc = Sc[t].r;
            for(uint32_t j = gs->p[z][v]; j < gs->p[z][v + 1]; j++) {
                uint32_t w = gs->e[z][j];
                uint32_t pw = Sc[t].pos(w);
                if(Sc[t].r <= pw && pw < Sc[t].c) {
                    Sc[t].swapByPos(cc++, pw);
                }
            }
            qcsz = (Sc[t].c - cc);
            memcpy(vqueue[t].data(), Sc[t].begin() + cc, sizeof(uint32_t) * qcsz);
        }
        
        tail = wsSize;
        for (int i = 0; i < tail; ++i) {
            uint32_t w = ws[deep][i];
            if (w == v) {
                --tail;
                ws[deep][i] = ws[deep][tail];
                ws[deep][tail] = w; i--;
            } else {
                uint32_t nbrs = 0;
                for (uint32_t j = 0; j < qcsz; ++j) {
                    uint32_t u1 = vqueue[t][j];
                    // assert(u1 < gs->n[t]);
                    if (is_nbr_sub(z, w, u1)) {
                        nbrs++; nw = u1;
                        if (nbrs > 1) break;
                    }
                }
                if (nbrs == 1) {
                    ws[deep][i] = ws[deep][front];
                    ws[deep][front] = w;
                    bin[t][front++] = nw;
                }
                else if (nbrs > 1) {
                    --tail;
                    ws[deep][i] = ws[deep][tail];
                    ws[deep][tail] = w; i--;
                }
            }
            // if (nbrs <= 1) subbcnt++;
            // if (nbrs == 0) continue;
            // else ws[deep][wsSize++] = v;
        }

        if (wsSize+front > ws[deep].capacity()) ws[deep].resize(wsSize+front+1);
        for (uint32_t i = 0; i < front; ++i) {
            ws[deep][wsSize+i] = bin[t][i];
            bin[t][i] = 0;
        }
        if (tail+1 < wsSize) std::sort(ws[deep].begin()+tail, ws[deep].begin()+wsSize, [&](int a, int b){return _deg[z][a] < _deg[z][b];});
    }
    if (wsSize == 1 && dv < mSize / dv) {
        Sc[z].swapByPos(--Sc[z].c, Sc[z].pos(ws[deep][0]));
    }
    else{
            // uint32_t rt = Sc[t].r, rz = Sc[z].r;
        for (uint32_t i = 0; i < front; ++i) {
            uint32_t v = ws[deep][i];
            uint32_t u1 = ws[deep][i+wsSize];
            Sc[z].swapByPos(Sc[z].r++, Sc[z].pos(v));
            Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u1));

            if(gs->deg(z, v) > 2*Sc[t].cSize()) {
                for(int j = Sc[t].r; j < Sc[t].c; j++) {
                    if(!is_nbr_sub(z, v, Sc[t][j])) {
                        Sc[t].swapByPos(--Sc[t].c, j); j--;
                    }
                }
            }
            else {
                uint32_t cc = Sc[t].r;
                for(uint32_t j = gs->p[z][v]; j < gs->p[z][v + 1]; j++) {
                    uint32_t w = gs->e[z][j];
                    uint32_t pw = Sc[t].pos(w);
                    if(Sc[t].r <= pw && pw < Sc[t].c) {
                        Sc[t].swapByPos(cc++, pw);
                    }
                }
                Sc[t].c = cc;
            }

            if(gs->deg(t, u1) > 2*Sc[z].cSize()) {
                for(int j = Sc[z].r; j < Sc[z].c; j++) {
                    if(!is_nbr_sub(t, u1, Sc[z][j])) {
                        Sc[z].swapByPos(--Sc[z].c, j); j--;
                    }
                }
            }
            else {
                uint32_t cc = Sc[z].r;
                for(uint32_t j = gs->p[t][u1]; j < gs->p[t][u1 + 1]; j++) {
                    uint32_t w = gs->e[t][j];
                    uint32_t pw = Sc[z].pos(w);
                    if(Sc[z].r <= pw && pw < Sc[z].c) {
                        Sc[z].swapByPos(cc++, pw);
                    }
                }
                Sc[z].c = cc;
            }

            uint32_t szr = Sc[z].r;
            branchSubG(deep+1, ne);
            Sc[t].c = ct; Sc[z].c = cz;
            // if (Sc[z][Sc[z].r-1] != v) {
            //     Sc[0].print();
            //     Sc[1].print();
            //     printf("z=%d, v=%d, Sc[z][Sc[z].r-1]=%d, szr=%d, Sc[z].r=%d, Sc[t].r=%d, Sc[t].c=%d\n", z, v, Sc[z][Sc[z].r-1], szr, Sc[z].r, Sc[t].r, Sc[t].c);
            // }
            // assert(Sc[z][Sc[z].r-1] == v);
            Sc[z].swapByPos(--Sc[z].r,--Sc[z].c); Sc[t].r--;
        }

        for (uint32_t i = tail; i < wsSize; ++i) {
            uint32_t v = ws[deep][i];

            Sc[z].swapByPos(Sc[z].r++, Sc[z].pos(v));
            if(gs->deg(z, v) > 2*Sc[t].cSize()) {
                for(int j = Sc[t].r; j < Sc[t].c; j++) {
                    if(!is_nbr_sub(z, v, Sc[t][j])) {
                        Sc[t].swapByPos(--Sc[t].c, j); j--;
                    }
                }
            }
            else {
                uint32_t cc = Sc[t].r;
                for(uint32_t j = gs->p[z][v]; j < gs->p[z][v + 1]; j++) {
                    uint32_t w = gs->e[z][j];
                    uint32_t pw = Sc[t].pos(w);
                    if(Sc[t].r <= pw && pw < Sc[t].c) {
                        Sc[t].swapByPos(cc++, pw);
                    }
                }
                Sc[t].c = cc;
            }

            uint32_t szr = Sc[z].r;
            branchSubG(deep+1, ne);
            Sc[t].c = ct;
            // if (Sc[z][Sc[z].r-1] != v) {
            //     Sc[0].print();
            //     Sc[1].print();
            //     printf("z=%d, v=%d, Sc[z][Sc[z].r-1]=%d, szr=%d, Sc[z].r=%d, Sc[t].r=%d, Sc[t].c=%d\n", z, v, Sc[z][Sc[z].r-1], szr, Sc[z].r, Sc[t].r, Sc[t].c);
            // }
            // assert(Sc[z][Sc[z].r-1] == v);
            Sc[z].swapByPos(--Sc[z].r,--Sc[z].c);
        }

        for (uint32_t i = front; i < tail; ++i) {
            uint32_t v = ws[deep][i];
            uint32_t pv = Sc[z].pos(v);
            // assert(pv >= Sc[z].r);
            // assert(pv < Sc[z].c);
            // if (pv < Sc[z].c)
            Sc[z].swapByPos(--Sc[z].c, pv);
        }
    }
    // pivot vertex
    Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u));

    branchSubG(deep+1, ne);

    Sc[t].r = r[t];
    Sc[t].c = c[t];
    Sc[z].r = r[z];
    Sc[z].c = c[z];

}

void MBCS::branch_minCover(uint32_t ux, uint32_t deep, uint32_t ne) {

    assert(Sc[1].r <= Sc[1].c);
    assert(Sc[0].r <= Sc[0].c);
    if (Sc[0].r > csz[0] || Sc[1].r > csz[1]) return;
    if (mSize >= (lrs[0]-Sc[0].r) * (lrs[1]- Sc[1].r)) return;
    // uint32_t xxxxxxxx = 0;
    // if (ux == 3520 ) {
    //     xxxxxxxx = 1;
    //     printf("deep1=%d, Sc[0].r=%d, Sc[0].c=%d, Sc[1].r=%d, Sc[1].c=%d\n", deep, Sc[0].r, Sc[0].c, Sc[1].r, Sc[1].c);
    //     printf("   U:");
    //     for (uint32_t i = 0; i < Sc[0].r; ++i)
    //         printf(" %d", gs->realLable[0][Sc[0][i]]);
    //     printf("\n");
    //     printf("   V:");
    //     for (uint32_t i = 0; i < Sc[1].r; ++i)
    //         printf(" %d", gs->realLable[1][Sc[1][i]]);
    //     printf("\n");
    // }

    numbranches ++;
    mxdeep = std::max(mxdeep, deep);
    if (Sc[0].cSize() == 0 || Sc[1].cSize() == 0) {
        if (mSize < (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].r)) {
            mSize = (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].r); realsz[0] = lrs[0]-Sc[0].r; realsz[1] = lrs[1]-Sc[1].r;
            // printf("mSize=%d\n", mSize);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].r; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].c; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
        }
        return;
    }

    uint32_t t = 0, z = 1;
    uint32_t r[2], c[2];
    r[0] = Sc[0].r; c[0] = Sc[0].c;
    r[1] = Sc[1].r; c[1] = Sc[1].c;

    if (Sc[1].cSize() < Sc[0].cSize()) {t = 1; z = 0;} 

    for(uint32_t i = Sc[z].r; i < Sc[z].c; ++i) _deg[z][Sc[z][i]] = 0;
    for(int i = Sc[t].r; i < Sc[t].c; i++) { //scan all vertices in C
        uint32_t u = Sc[t][i];
        uint32_t e = 0;
        if(gs->deg(t, u) > 4*Sc[z].cSize()) {
            for(uint32_t j = Sc[z].r; j < Sc[z].c; j++) {
                uint32_t v = Sc[z][j];
                if(is_nbr_cvert(t, u, v)) {
                    e++; _deg[z][v]++;
                }
            }
        }
        else {
            for(uint32_t j = gs->p[t][u]; j < gs->p[t][u + 1]; j++) {
                uint32_t v = gs->e[t][j];
                uint32_t pv = Sc[z].pos(v);
                if(Sc[z].r <= pv && pv < Sc[z].c) {//in C
                    e++; _deg[z][v]++;
                }
            }
        }
        _deg[t][u] = e;
    }

    #ifdef UB_HNBR_
    while (true) {
        uint32_t qr[2] = {0};
        for (uint32_t ttt = 0; ttt <= 1; ttt++) {
            uint32_t zzz = ttt^1;
            for(int i = Sc[ttt].r; i < Sc[ttt].c; i++) {
                uint32_t v = Sc[ttt][i];
                uint32_t dv = _deg[ttt][v];
                if (dv == 0) Sc[ttt].swapByPos(--Sc[ttt].c, i--);
                else if (dv + Sc[zzz].r > csz[zzz]) {
                    Sc[ttt].swapByPos(Sc[ttt].r++, i); qr[ttt]++;

                    if (Sc[ttt].r > csz[ttt]) {
                            Sc[0].r = r[0];
                            Sc[1].r = r[1];
                            Sc[0].c = c[0];
                            Sc[1].c = c[1];
                            return; 
                    }

                    if (gs->deg(ttt, v) > 2*Sc[zzz].cSize())
                    for (uint32_t j = Sc[zzz].r; j < Sc[zzz].c; ++j) {
                        uint32_t w = Sc[zzz][j];
                        if (is_nbr_cvert(ttt, v, w)) _deg[zzz][w]--;
                    }
                    else
                    for (uint32_t j = gs->p[ttt][v]; j < gs->p[ttt][v+1]; ++j) {
                        uint32_t w = gs->e[ttt][j];
                        uint32_t wp = Sc[zzz].pos(w);
                        if (wp >= Sc[zzz].r && wp < Sc[zzz].c) {
                            _deg[zzz][w]--;
                        }
                    }
                }
            }
        }

        if (Sc[0].r > csz[0] || Sc[1].r > csz[1]) {
            Sc[0].r = r[0];
            Sc[1].r = r[1];
            Sc[0].c = c[0];
            Sc[1].c = c[1];
            return; 
        }
        if (qr[0] == 0 && qr[1] == 0) break;
    }
    #endif

    uint32_t minE[2], maxE[2] = {0}, maxU[2] = {0}, minU[2] = {0}, numes = 0;
    minE[0] = gs->n[0]; minE[1] = gs->n[1];
    for(int i = Sc[t].r; i < Sc[t].c; i++) {
        uint32_t v = Sc[t][i];
        uint32_t dv = _deg[t][v];
        numes += dv;
        if (dv > maxE[t]) {maxE[t] = dv; maxU[t] = v;}
        if (dv < minE[t]) {minE[t] = dv; minU[t] = v;}
    }

    for(int i = Sc[z].r; i < Sc[z].c; i++) {
        uint32_t v = Sc[z][i];
        uint32_t dv = _deg[z][v];
        if (dv > maxE[z]) { maxE[z] = dv; maxU[z] = v;}
        if (dv < minE[t]) { minE[z] = dv; minU[z] = v;}
    }

    if (mSize < (lrs[0]-Sc[0].r) * (lrs[1]- Sc[1].c) && Sc[1].c <= csz[1]) {
        mSize = (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].c); realsz[0] = lrs[0]-Sc[0].r; realsz[1] = lrs[1]-Sc[1].c;
        // printRes();
        // printf("mSize=%d, %d %d\n", mSize, lrs[0]-Sc[0].r, lrs[1]-Sc[1].c);
    }
    else if (mSize < (lrs[0]-Sc[0].c) * (lrs[1]- Sc[1].r) && Sc[0].c <= csz[0]) {
        mSize = (lrs[0]-Sc[0].c) * (lrs[1]-Sc[1].r); realsz[0] = lrs[0]-Sc[0].c; realsz[1] = lrs[1]-Sc[1].r;
        // printf("mSize=%d, %d %d\n", mSize, lrs[0]-Sc[0].c, lrs[1]-Sc[1].r);
        // printRes();
    }

    #ifdef REDUCTION_
    if (maxE[0] == 0) {
        assert(maxE[1] == 0);
        if (mSize < (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].r)) {
            mSize = (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].r); realsz[0] = lrs[0]-Sc[0].r; realsz[1] = lrs[1]-Sc[1].r;
            // printf("mSize=%d\n", mSize);
            // realRes[0].clear(); realRes[1].clear();
            // for (uint32_t v = 0; v < Sc[0].r; ++v)
            //     realRes[0].emplace_back(Sc[0][v]);
            // for (uint32_t v = 0; v < Sc[1].c; ++v)
            //     realRes[1].emplace_back(Sc[1][v]);
        }
        Sc[0].r = r[0];
        Sc[1].r = r[1];
        Sc[0].c = c[0];
        Sc[1].c = c[1];
        return;
    }

    if (maxE[0] <= 1 || maxE[1] <= 1) {

        if (maxE[t] < maxE[z]) {
            uint32_t qsz[2] = {0};
            for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
                uint32_t d = _deg[z][Sc[z][i]];
                vqueue[z][qsz[z]++] = d;
                // assert(d >= 1);
            }

            // sort(vqueue[z].begin(), vqueue[z].begin()+qsz[z]);
            sort_in(vqueue[z], z, 0, qsz[z]);

            uint32_t qr = 0, nums = 0;
            uint32_t a[2] = {0};
            a[t] = lrs[t] - Sc[t].r;
            a[z] = lrs[z] - Sc[z].c;

            if (mSize < a[t] * a[z] && a[t] >= sz[t] && a[z] >= sz[z]) {
                mSize = a[t] * a[z];
                realsz[t] = a[t];
                realsz[z] = a[z];
            }

            while (qr < qsz[z]) {
                uint32_t d = vqueue[z][qr++];
                a[z] ++;
                a[t] -= d;
                nums += d;
                if (mSize < a[t] * a[z] && a[t] >= sz[t] && a[z] >= sz[z]) {
                    mSize = a[t] * a[z];
                    realsz[t] = a[t];
                    realsz[z] = a[z];
                }
            }
            // assert(nums+Sc[t].r == Sc[t].c);
        }
        else if (maxE[t] > maxE[z]) {
            uint32_t qsz[2] = {0};
            for (uint32_t i = Sc[t].r; i < Sc[t].c; ++i) {
                uint32_t d = _deg[t][Sc[t][i]];
                vqueue[t][qsz[t]++] = d;
                // assert(d >= 1);
            }

            // sort(vqueue[t].begin(), vqueue[t].begin()+qsz[t]);
            sort_in(vqueue[t], t, 0, qsz[t]);

            uint32_t qr = 0, nums = 0;
            uint32_t a[2] = {0};
            a[t] = lrs[t] - Sc[t].c;
            a[z] = lrs[z] - Sc[z].r;

            if (mSize < a[t] * a[z] && a[t] >= sz[t] && a[z] >= sz[z]) {
                mSize = a[t] * a[z];
                realsz[t] = a[t];
                realsz[z] = a[z];
            }

            while (qr < qsz[t]) {
                uint32_t d = vqueue[t][qr++];
                a[t] ++;
                a[z] -= d;
                nums += d;
                if (mSize < a[t] * a[z] && a[t] >= sz[t] && a[z] >= sz[z]) {
                    mSize = a[t] * a[z];
                    realsz[t] = a[t];
                    realsz[z] = a[z];
                }
            }
            // assert(nums+Sc[z].r == Sc[z].c);
        }
        else {
            uint32_t a[2] = {0}, sne = 0;
            a[t] = lrs[t] - Sc[t].r;
            a[z] = lrs[z] - Sc[z].r;

            for(int i = Sc[t].r; i < Sc[t].c; i++) {

                uint32_t v = Sc[t][i];
                uint32_t dv = _deg[t][v];
                if (dv == 1) sne++;
            }
            uint32_t e = sne;
            if (a[t] >= a[z]) {
                uint32_t dr = std::min(sne, a[t]-a[z]);
                sne -= dr;
                a[t] -= (dr + sne/2);
                a[z] -= (sne - sne/2);
            }
            else {
                uint32_t dr = std::min(sne, a[z]-a[t]);
                sne -= dr;
                a[z] -= (dr + sne/2);
                a[t] -= (sne - sne/2);
            }
            if (a[t] < sz[t]) {a[z] -= (sz[t]-a[t]); a[t] = sz[t];}
            else if (a[z] < sz[z]) {a[t] -= (sz[z]-a[z]); a[z] = sz[z];}

            if (mSize < a[t] * a[z] && a[t] >= sz[t] && a[z] >= sz[z]) {
                mSize = a[t] * a[z];
                realsz[t] = a[t];
                realsz[z] = a[z];
            }
        }

        Sc[0].r = r[0];
        Sc[1].r = r[1];
        Sc[0].c = c[0];
        Sc[1].c = c[1];
        return;
    }
    #endif
    // if (numes > 2*(csz[t]-Sc[t].r)*(csz[z]-Sc[t].r) || Sc[t].cSize() + Sc[z].cSize() > 4*(csz[t]-Sc[t].r)*(csz[z]-Sc[t].r) ) {
    //     Sc[0].r = r[0];
    //     Sc[1].r = r[1];
    //     Sc[0].c = c[0];
    //     Sc[1].c = c[1];
    //     return;  
    // }

    #ifdef UB_HNBR_

        if (Sc[z].c < Sc[t].c) { t = z; z = t^1; }

        uint32_t qsz[2] = {0};
        uint32_t ub[2] = {0}, sne[2] = {0};
        // uint32_t mSize_ = mSize == 0 ? (lrs[t]-sz[t]) * (lrs[z]-sz[z]) : mSize;;
        ub[t] = lrs[t] - Sc[t].r; ub[z] = lrs[z] - Sc[z].r;

        bool nextb = false;
        if ((ub[t] >= sz[t] && ub[z]-Sc[z].cSize() >= sz[z] && mSize < ub[t] * (ub[z]-Sc[z].cSize())) || (ub[t]-Sc[t].cSize() >= sz[t] && ub[z] >= sz[z] && mSize < (ub[t]-Sc[t].cSize()) * ub[z]))
            nextb = true;

        if (!nextb) {
            uint32_t nume[2] = {0};
            for (uint32_t i = Sc[t].r; i < Sc[t].c; ++i) {
                uint32_t d = _deg[t][Sc[t][i]];
                vqueue[t][qsz[t]++] = d;
                nume[t] += d;
                // assert(d >= 1);
            }
            for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
                uint32_t d = _deg[z][Sc[z][i]];
                vqueue[z][qsz[z]++] = d;
                nume[z] += d;
                // assert(d >= 1);
            }
            // assert(nume[t] == nume[z]);
            // sort(vqueue[t].begin(), vqueue[t].begin()+qsz[t], [](uint32_t a, uint32_t b){return a > b;});
            // sort(vqueue[z].begin(), vqueue[z].begin()+qsz[z], [](uint32_t a, uint32_t b){return a > b;});
            sort_de(vqueue[t], t, 0, qsz[t]);
            sort_de(vqueue[z], z, 0, qsz[z]);
            int jz = 0;
            for (uint32_t i = 0; i < qsz[t]; ++i) {
                uint32_t d = vqueue[t][i];
                // if (i > 0) assert(vqueue[t][i-1] >= vqueue[t][i]);
                if (Sc[t].r+i+1 > csz[t]) break;

                sne[t] += d;
                ub[t]--;
                if (sne[z] + sne[t] < numes) {
                    for (uint32_t j = 0; j < qsz[z]; ++j) {
                        sne[z] += vqueue[z][j];
                        ub[z]--; jz++;
                        if (sne[t]+sne[z] >= numes) break;
                    }
                }
                else {
                    for (int j = jz-1; j >= 0 ; --j) {
                        sne[z] -= vqueue[z][j];
                        ub[z]++; jz--;
                        if (j == 0) break;
                        if (sne[t]+sne[z] < numes) { 
                            sne[z] += vqueue[z][j];
                            ub[z]--; jz++; 
                            break;
                        }
                        if ((sne[t]+sne[z] >= numes) && (sne[t]+sne[z]-vqueue[z][j-1] < numes)) break;
                    }
                }
               

                if (ub[t] >= sz[t] && ub[z] >= sz[z] && mSize < ub[t] * ub[z]) {
                    uint32_t uclr = 0;

                    nextb = true; break;
                }
            }

            if (!nextb) {
                Sc[0].r = r[0];
                Sc[1].r = r[1];
                Sc[0].c = c[0];
                Sc[1].c = c[1];
                return;
            }

            // t = z; z = t^1;
            // // qsz[0] = 0; qsz[1] = 0;
            // ub[0] = 0; ub[1] = 0;
            // sne[0] = 0; sne[1] = 0;
            // // nume[0] = 0; nume[1] = 0;
            // ub[t] = lrs[t] - Sc[t].r; ub[z] = lrs[z] - Sc[z].r;

            // nextb = false;
            // if ((ub[t] >= sz[t] && ub[z]-Sc[z].cSize() >= sz[z] && mSize < ub[t] * (ub[z]-Sc[z].cSize())) || (ub[t]-Sc[t].cSize() >= sz[t] && ub[z] >= sz[z] && mSize < (ub[t]-Sc[t].cSize()) * ub[z]))
            //     nextb = true;

            // if (!nextb) {
            //     // for (uint32_t i = Sc[t].r; i < Sc[t].c; ++i) {
            //     //     uint32_t d = _deg[t][Sc[t][i]];
            //     //     vqueue[t][qsz[t]++] = d;
            //     //     nume[t] += d;
            //     //     assert(d >= 1);
            //     // }
            //     // for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
            //     //     uint32_t d = _deg[z][Sc[z][i]];
            //     //     vqueue[z][qsz[z]++] = d;
            //     //     nume[z] += d;
            //     //     assert(d >= 1);
            //     // }
            //     int jz = 0;
            //     for (uint32_t i = 0; i < qsz[t]; ++i) {
            //         uint32_t d = vqueue[t][i];
            //         if (i > 0) assert(vqueue[t][i-1] >= vqueue[t][i]);
            //         if (Sc[t].r+i+1 > csz[t]) break;

            //         sne[t] += d;
            //         ub[t]--;
            //         if (sne[z] + sne[t] < numes) {
            //             for (uint32_t j = 0; j < qsz[z]; ++j) {
            //                 sne[z] += vqueue[z][j];
            //                 ub[z]--; jz++;
            //                 if (sne[t]+sne[z] >= numes) break;
            //             }
            //         }
            //         else {
            //             for (int j = jz-1; j >= 0 ; --j) {
            //                 sne[z] -= vqueue[z][j];
            //                 ub[z]++; jz--;
            //                 if (j == 0) break;
            //                 if (sne[t]+sne[z] < numes) { 
            //                     sne[z] += vqueue[z][j];
            //                     ub[z]--; jz++; 
            //                     break;
            //                 }
            //                 if ((sne[t]+sne[z] >= numes) && (sne[t]+sne[z]-vqueue[z][j-1] < numes)) break;
            //             }
            //         }

            //         // if (ux == 8 && deep == 14) {
            //         //     printf("i=%d, sne[t]=%d, sne[z]=%d, ub[t]=%d, ub[z]=%d, jz=%d\n", i, sne[t], sne[z], ub[t], ub[z], jz);
            //         // }

            //         if (ub[t] >= sz[t] && ub[z] >= sz[z] && mSize < ub[t] * ub[z]) {nextb = true; break;}
            //     }

            //     if (!nextb) {
            //         Sc[0].r = r[0];
            //         Sc[1].r = r[1];
            //         Sc[0].c = c[0];
            //         Sc[1].c = c[1];
            //         return;
            //     }
            // }
        }
    #endif

    if (mSize < (lrs[0]-Sc[0].r) * (lrs[1]- Sc[1].c) && (lrs[0]-Sc[0].r) >= sz[0] && (lrs[1]-Sc[1].c) >= sz[1]) {
        mSize = (lrs[0]-Sc[0].r) * (lrs[1]-Sc[1].c); realsz[0] = lrs[0]-Sc[0].r; realsz[1] = lrs[1]-Sc[1].c;
    }
    else if (mSize < (lrs[0]-Sc[0].c) * (lrs[1]- Sc[1].r) && (lrs[0]-Sc[0].c) >= sz[0] && (lrs[1]-Sc[1].r) >= sz[1]) {
        mSize = (lrs[0]-Sc[0].c) * (lrs[1]-Sc[1].r); realsz[0] = lrs[0]-Sc[0].c; realsz[1] = lrs[1]-Sc[1].r;
    }

    // if (maxE[t] <= 2 && maxE[z] <= 2) {
    //     uint32_t lines[4] = {0};

    //     Sc[0].r = r[0];
    //     Sc[1].r = r[1];
    //     Sc[0].c = c[0];
    //     Sc[1].c = c[1];
    //     return; 
    // }

    // if (csz[t] > csz[z]){t = z; z = t^1;}

    if (maxE[t] < maxE[z]) {t = z; z = t^1;}
    // if (Sc[t].cSize() - maxE[z] > Sc[z].cSize() - maxE[t]) {t = z; z = t^1;}

    uint32_t u = maxU[t], Str = Sc[t].r;
    uint32_t ud = _deg[t][u];
    Sc[t].swapByPos(Sc[t].r++, Sc[t].pos(u));

    branch_minCover(ux, deep+1, ne);
    Sc[t].r = Str; 
    // Sc[t].swapByPos(--Sc[t].c, Sc[t].pos(u));
    if (gs->deg(t, u) > 4*Sc[z].cSize()) {
        for (uint32_t i = Sc[z].r; i < Sc[z].c; ++i) {
            uint32_t v = Sc[z][i];
            if (is_nbr_cvert(t, u, v)) {
                Sc[z].swapByPos(Sc[z].r++, Sc[z].pos(v));
            }
        }
    }
    else {
        for (uint32_t i = gs->p[t][u]; i < gs->p[t][u+1]; ++i) {
            uint32_t v = gs->e[t][i];
            uint32_t vp = Sc[z].pos(v);
            if (vp >= Sc[z].r && vp < Sc[z].c) {
                Sc[z].swapByPos(Sc[z].r++, vp);
            }
        }
    }
    branch_minCover(ux,deep+1, ne);
    Sc[0].r = r[0];
    Sc[1].r = r[1];
    Sc[0].c = c[0];
    Sc[1].c = c[1];
}

bool MBCS::is_nbr_sub(uint32_t t, uint32_t u, uint32_t v) {

    return g->is_nbr(t,gs->realLable[t][u], gs->realLable[t^1][v]);
}

bool MBCS::is_nbr_cvert(uint32_t t, uint32_t u, uint32_t v) {
    return g->is_nbr(t,gs->realLable[t][u], gs->realLable[t^1][v])^1;
}

void MBCS::graph_convert()
{
    if (gs == NULL) {
        gs = new biGraph();
        gs->n[0] = g->n[0];
        gs->n[1] = g->n[1];
        gs->m = g->m;

        gs->p[0].resize(gs->n[0] + 5);
        gs->p[1].resize(gs->n[1] + 5);
        gs->e[0].resize(gs->m+5);
        gs->e[1].resize(gs->m+5);
        gs->realLable[0].resize(gs->n[0] + 5);
        gs->realLable[1].resize(gs->n[1] + 5);
    }

    for (uint32_t t = 0; t <= 1; ++t) {
        uint32_t z = t^1;
        for (uint32_t i = 0; i < S[t].c; ++i){
            uint32_t v = S[t][i];
            uint32_t d = _deg[t][v];
            // if (d < sz[z]) {
            //     printf("v =%d, d = %d,  _deg[t][v]=%d\n", v, d, _deg[t][v]);
            //     printf("S[t].c=%d, S[z].c=%d, S[t][0]=%d, S[z][0]=%d\n", S[t].c, S[z].c, S[t][0], S[z][0]);
            // }
            // assert(d >= sz[z]);
            uint32_t x = d;
            d = S[z].c - x;
            if (d > S[z].c) {
                printf("v =%d, d = %d, S[z].c=%d, _deg[t][v]=%d\n", v, d, S[z].c, _deg[t][v]);
                printf("S[t].r=%d, S[z].r=%d, S[t][0]=%d, S[z][0]=%d\n", S[t].r, S[z].r, S[t][0], S[z][0]);
            }
            assert(d <= S[z].c);
            gs->realLable[t][i] = v;
        }
    }
    uint32_t t = 0, z = 1;
    gs->p[t][0] = 0;
    for (uint32_t i = 0; i < S[t].c; ++i){
        uint32_t v = S[t][i];
        uint32_t d = _deg[t][v];
        uint32_t x = d;
        d = S[z].c - x;
        gs->p[t][i+1] = gs->p[t][i] + d;
        uint32_t cr = 0;
        for (uint32_t j = 0; j < S[z].c; ++j) {
            uint32_t u = S[z][j];
            if (!g->is_nbr(t,v, u)) gs->e[t][gs->p[t][i]+cr++] = j;
        }
        if (cr != d) {
            printf("v=%d, d=%d, cr=%d\n", v, d, cr);
        }
        assert(cr == d);
    }

    gs->p[z][0] = 0;
    for (uint32_t i = 0; i < S[z].c; ++i){
        uint32_t v = S[z][i];
        uint32_t d = _deg[z][v];
        uint32_t x = d;
        d = S[t].c - x;
        gs->p[z][i+1] = gs->p[z][i] + d;
    }

    for (uint32_t i = 0; i < S[t].c; ++i){
        for (uint32_t j = gs->p[t][i]; j < gs->p[t][i+1]; ++j) {
            uint32_t u = gs->e[t][j];
            gs->e[z][gs->p[z][u]++] = i;
        }
    }
    for (uint32_t i = S[z].c; i > 0 ; --i){
        gs->p[z][i] = gs->p[z][i-1];
    }
    gs->p[z][0] = 0;

    for (uint32_t i = 0; i < S[z].c; ++i){
        uint32_t v = S[z][i];
        uint32_t d = _deg[z][v];
        uint32_t x = d;
        d = S[t].c - x;
        if (gs->p[z][i+1] != (gs->p[z][i] + d)) {
            printf("v=%d, d=%d, gs->p[z][i+1]=%d, gs->p[z][i]=%d\n", v, d, gs->p[z][i+1], gs->p[z][i]);
            for (uint32_t i = 0; i < 10; ++i)
                printf("i=%d, p[z][i]=%d, d=%d\n", i, gs->p[z][i], d);
        }
        assert(gs->p[z][i+1] == (gs->p[z][i] + d));
    }

    for (uint32_t i = 0; i < S[t].c; ++i){
        Sc[t].swapByPos(i, Sc[t].pos(i));
    }
    for (uint32_t i = 0; i < S[z].c; ++i){
        Sc[z].swapByPos(i, Sc[z].pos(i));
    }
    Sc[t].c = S[t].c; Sc[z].c = S[z].c;
    Sc[t].r = 0; Sc[z].r = 0;

    // //test
    // for (uint32_t i = 0; i < Sc[t].c; ++i){
    //     uint32_t v = gs->realLable[t][i];
    //     for (uint32_t j = 0; j < Sc[z].c; ++j){
    //         uint32_t u = gs->realLable[z][j];
    //         if (g->is_nbr(t, v, u)) assert(!is_nbr_cvert(t, i,j));
    //         else assert(is_nbr_cvert(t, i,j));
    //     }
    //     for (uint32_t j = gs->p[t][i]; j < gs->p[t][i+1]; ++j) {
    //         uint32_t u = gs->realLable[z][gs->e[t][j]];
    //         assert(is_nbr_cvert(t, i, gs->e[t][j]));
    //     }
    // }
}

void MBCS::subGraph_generate()
{
    if (gs == NULL) {
        gs = new biGraph();
        gs->n[0] = g->n[0];
        gs->n[1] = g->n[1];
        gs->m = g->m;

        gs->p[0].resize(gs->n[0] + 5);
        gs->p[1].resize(gs->n[1] + 5);
        gs->e[0].resize(gs->m+5);
        gs->e[1].resize(gs->m+5);
        gs->realLable[0].resize(gs->n[0] + 5);
        gs->realLable[1].resize(gs->n[1] + 5);
    }

    for (uint32_t t = 0; t <= 1; ++t) {
        uint32_t z = t^1;
        for (uint32_t i = 0; i < S[t].c; ++i){
            uint32_t v = S[t][i];
            uint32_t d = _deg[t][v];
            assert(d >= 0);
            // if (d > S[z].c) {
            //     printf("v =%d, d = %d, S[z].c=%d, _deg[t][v]=%d\n", v, d, S[z].c, _deg[t][v]);
            //     printf("S[t].r=%d, S[z].r=%d, S[t][0]=%d, S[z][0]=%d\n", S[t].r, S[z].r, S[t][0], S[z][0]);
            //     assert(d <= S[z].c);
            // }
            gs->realLable[t][i] = v;
        }
    }
    uint32_t t = 0, z = 1, e = 0;
    gs->p[t][0] = 0;
    for (uint32_t i = 0; i < S[t].c; ++i){
        uint32_t v = S[t][i];
        uint32_t d = _deg[t][v];;
        uint32_t cr = 0;
        gs->p[t][i+1] = gs->p[t][i] + d;
        for (uint32_t j = 0; j < S[z].c; ++j) {
            uint32_t u = S[z][j];
            if (g->is_nbr(t,v, u)) gs->e[t][gs->p[t][i]+cr++] = j;
        }
        // if (cr != d) {
        //     printf("v=%d, d=%d, cr=%d\n", v, d, cr);
        //     assert(cr == d);
        // }
        e += cr;
    }

    gs->p[z][0] = 0;
    for (uint32_t i = 0; i < S[z].c; ++i){
        uint32_t v = S[z][i];
        uint32_t d = _deg[z][v];;
        gs->p[z][i+1] = gs->p[z][i] + d;
    }
    assert(e == gs->p[z][S[z].c]);
    for (uint32_t i = 0; i < S[t].c; ++i){
        for (uint32_t j = gs->p[t][i]; j < gs->p[t][i+1]; ++j) {
            uint32_t u = gs->e[t][j];
            gs->e[z][gs->p[z][u]++] = i;
        }
    }
    for (uint32_t i = S[z].c; i > 0 ; --i){
        gs->p[z][i] = gs->p[z][i-1];
    }
    gs->p[z][0] = 0;

    for (uint32_t i = 0; i < S[z].c; ++i){
        uint32_t v = S[z][i];
        uint32_t d = _deg[z][v];;
        if (gs->p[z][i+1] != (gs->p[z][i] + d)) {
            printf("v=%d, d=%d, gs->p[z][i+1]=%d, gs->p[z][i]=%d\n", v, d, gs->p[z][i+1], gs->p[z][i]);
            for (uint32_t i = 0; i < 10; ++i)
                printf("i=%d, p[z][i]=%d, d=%d\n", i, gs->p[z][i], d);
        }
        assert(gs->p[z][i+1] == (gs->p[z][i] + d));
    }

    for (uint32_t i = 0; i < S[t].c; ++i){
        Sc[t].swapByPos(i, Sc[t].pos(i));
    }
    for (uint32_t i = 0; i < S[z].c; ++i){
        Sc[z].swapByPos(i, Sc[z].pos(i));
    }
    // Sc[t].c = S[t].c; Sc[t].r = 0;
    // Sc[z].c = S[z].c; Sc[z].r = 0;
    Sc[t].c = S[t].c; Sc[t].r = S[t].r;
    Sc[z].c = S[z].c; Sc[z].r = S[z].r;

    // //test
    // for (uint32_t i = 0; i < Sc[t].c; ++i){
    //     uint32_t v = gs->realLable[t][i];
    //     for (uint32_t j = 0; j < Sc[z].c; ++j){
    //         uint32_t u = gs->realLable[z][j];
    //         if (g->is_nbr(t, v, u)) assert(is_nbr_sub(t, i,j));
    //         else assert(!is_nbr_sub(t, i,j));
    //     }
    //     for (uint32_t j = gs->p[t][i]; j < gs->p[t][i+1]; ++j) {
    //         uint32_t u = gs->realLable[z][gs->e[t][j]];

    //         assert(is_nbr_sub(t, i, gs->e[t][j]));
    //     }
    // }
}

void MBCS::sort_deg(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t) {

    uint32_t maxd = 0, z = l^1;
    uint32_t vs = 0, len = 0;
    vs = len = t-s; 
    assert(t > s >= 0);
    if (ws[0].size() < t - s) ws[0].resize(t-s);
    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = _deg[l][vec[i]];
        maxd = std::max(d, maxd);
        ws[0][--len] = vec[i];
    }

    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = _deg[l][vec[i]];
        if (d >= bin[z].size()) printf("d = %d, l=%d, bin[l].size()=%d\n", d, l , bin[l].size());
        assert(d < bin[z].size());
        bin[z][d]++;
    }

    // for (uint32_t i = 0; i <= maxd; ++i) {
    //     uint32_t c = bin[z][i];
    //     bin[z][i] = s;
    //     s += c;
    // }

    for (uint32_t i = maxd; i >= 0; --i) {
        uint32_t c = bin[z][i];
        bin[z][i] = s;
        s += c;
        if (i == 0) break;
    }

    for (uint32_t i = 0; i < vs; ++i) {
        uint32_t v = ws[0][i];
        uint32_t pos = bin[z][_deg[l][v]]++;
        vec[pos] = v;
    }

    for (uint32_t i = 0; i <= maxd; ++i) bin[z][i] = 0;
}

void MBCS::sort_de(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t) {

    uint32_t maxd = 0, z = l^1;
    uint32_t vs = 0, len = 0;
    vs = len = t-s; 
    assert(t > s >= 0);
    if (ws[0].size() < t - s) ws[0].resize(t-s);
    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = vec[i];
        maxd = std::max(d, maxd);
        ws[0][--len] = vec[i];
    }

    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = vec[i];
        if (d >= bin[z].size()) printf("d = %d, l=%d, bin[l].size()=%d\n", d, l , bin[l].size());
        assert(d < bin[z].size());
        bin[z][d]++;
    }

    // for (uint32_t i = 0; i <= maxd; ++i) {
    //     uint32_t c = bin[z][i];
    //     bin[z][i] = s;
    //     s += c;
    // }

    for (uint32_t i = maxd; i >= 0; --i) {
        uint32_t c = bin[z][i];
        bin[z][i] = s;
        s += c;
        if (i == 0) break;
    }

    for (uint32_t i = 0; i < vs; ++i) {
        uint32_t v = ws[0][i];
        uint32_t pos = bin[z][v]++;
        vec[pos] = v;
    }

    for (uint32_t i = 0; i <= maxd; ++i) bin[z][i] = 0;
}

void MBCS::sort_in(std::vector<uint32_t> &vec, uint32_t l, uint32_t s, uint32_t t) {

    uint32_t maxd = 0, z = l^1;
    uint32_t vs = 0, len = 0;
    vs = len = t-s; 
    assert(t > s >= 0);
    if (ws[0].size() < t - s) ws[0].resize(t-s);
    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = vec[i];
        maxd = std::max(d, maxd);
        ws[0][--len] = vec[i];
    }

    for (uint32_t i = s; i < t; ++i) {
        uint32_t d = vec[i];
        if (d >= bin[z].size()) printf("d = %d, l=%d, bin[l].size()=%d\n", d, l , bin[l].size());
        assert(d < bin[z].size());
        bin[z][d]++;
    }

    for (uint32_t i = 0; i <= maxd; ++i) {
        uint32_t c = bin[z][i];
        bin[z][i] = s;
        s += c;
    }

    // for (uint32_t i = maxd; i >= 0; --i) {
    //     uint32_t c = bin[z][i];
    //     bin[z][i] = s;
    //     s += c;
    //     if (i == 0) break;
    // }

    for (uint32_t i = 0; i < vs; ++i) {
        uint32_t v = ws[0][i];
        uint32_t pos = bin[z][v]++;
        vec[pos] = v;
    }

    for (uint32_t i = 0; i <= maxd; ++i) bin[z][i] = 0;
}