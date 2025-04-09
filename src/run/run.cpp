#include "../MBCS/MBCS.h"
#include "../tools/getArgs.hpp"

int main(int argc, char* argv[]) {
    argsController ac(argc, argv);

    if(!ac.exist("-f")) {
        printf("file path: -f\n");
        return 0;
    }

    // if(!ac.exist("-d")) {
    //     printf("graph order: -d(core, two)\n");
    //     return 0;
    // }

    if(!ac.exist("-a")) {
        printf("run algorithm : -a(1, 2)\n");
        return 0;
    }


    uint32_t ls = 1, rs = 1;
    if(ac.exist("-l")) ls = std::stoi(ac["-l"]);
    if(ac.exist("-r")) rs = std::stoi(ac["-r"]);

    std::string fPath = ac["-f"];
    // std::string order = ac["-d"];
    std::string runAlg = ac["-a"];

    std::cout << "file: " << fPath << " Alg: " << runAlg << std::endl;
    std::cout << "l: " << ls << " r: " << rs << std::endl;
    // std::cout << "r " << rs << std::endl;

  
    MBCS mbs(fPath, ls, rs);
    // if (ac.exist("-o")) bce.setCounts(std::stoll(ac["-o"]));

    if(ac.exist("-s")) mbs.setStp(std::stof(ac["-s"]));

    auto t1 = std::chrono::steady_clock::now();
    
    uint32_t alg = atoi(runAlg.c_str());
    if (alg == 1 || alg == 0) 
        mbs.run_basic(alg);
    else {
        uint32_t alg = atoi(runAlg.c_str());
        assert(alg == 2 || alg == 3 || alg == 4);
        mbs.run(alg);
    }

    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Total time:" << duration.count() << "ms" << std::endl;
    

    return 0;
}