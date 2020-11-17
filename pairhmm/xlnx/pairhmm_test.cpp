#include <vector>
#include <string>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <cmath>
#include <random>
#include "host/FalconPairHMM.h"


double avg_numReads;
double avg_numHaps;
int case_counter;




char GenRandBase(){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,3);
    int number = distribution(generator);
    if(number == 0)
        return 'A';
    else if(number == 1)
        return 'T';
    else if (number == 2)
        return 'C';
    else
        return 'G';
}

int GenQuals(){
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(30.0,5.0);
    int quals = distribution(generator);
    if(quals < 6)
        quals = 6;
    return quals;
}

int GenInDel(){
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(40.0,1.0);
    int quals = distribution(generator);
    if(quals < 1)
        quals = 1;
    return quals;
}

int GenLen(int limit){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(limit / 4,limit);
    int number = distribution(generator);
    return number;
}

int GenInputs(pairhmmInput* in, int size){
    in->reads.clear();
    in->haps.clear();
    in->reads.resize(16 * (size + 1));
    in->haps.resize((size + 1));
    for(int i = 0; (size_t)i < in->reads.size(); i++){
        Read& curRead = in->reads[i];
        for(int j = 0; j < GenLen(MAX_READ_LEN); j++){
            curRead.bases.push_back(GenRandBase());
            curRead._q.push_back(GenQuals());
            curRead._i.push_back(GenInDel());
            curRead._d.push_back(GenInDel());
            curRead._c.push_back(10);
        }
    }
    for(int i = 0; (size_t)i < in->haps.size(); i++){
        Hap& curHap = in->haps[i];
        for(int j = 0; j < GenLen(MAX_HAP_LEN); j++){
            curHap.bases.push_back(GenRandBase()); 
        }
    }
    return 0;
}

int GenOutputs(pairhmmInput* in, pairhmmOutput* out){
    FalconPairHMM* falcon = new FalconPairHMM();
    falcon->computePairhmmAVX(in, out, false);
    delete falcon;
    return 0;
}

int GetInputs(pairhmmInput* in, std::string filename){
    std::ifstream ifs(filename.c_str(), std::ifstream::in);
    int numReads;
    int numHaplotypes;
    //thie first line is number of reads and number of haplotypes
    char lineBuf[1024];
    if(!ifs.good())
        printf("bad file name %s\n", filename.c_str());
    ifs.getline(lineBuf, 1024);
    char* token;
    token = strtok(lineBuf, " ");
    token = strtok(NULL, " ");
    numReads = atoi(token);
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    numHaplotypes = atoi(token);
   //start from the second line are all the reads
    for(int i = 0; i < numReads; ++i){
        int curReadLen;
        Read curRead;
        ifs.getline(lineBuf, 1024);
        curReadLen = atoi(lineBuf);
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curReadLen; j++){
            if(j == 0)
                token = strtok(lineBuf, " ");
            else
                token = strtok(NULL, " ");
            curRead.bases.push_back((char)atoi(token));
        }
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curReadLen; j++){
            if(j == 0)
                token = strtok(lineBuf, " ");
            else
                token = strtok(NULL, " ");
            curRead._q.push_back((char)atoi(token));
        }
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curReadLen; j++){
            if(j == 0)
                token = strtok(lineBuf, " ");
            else
                token = strtok(NULL, " ");
            curRead._i.push_back((char)atoi(token));
        }
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curReadLen; j++){
            if(j == 0)
                token = strtok(lineBuf, " ");
            else
                token = strtok(NULL, " ");
            curRead._d.push_back((char)atoi(token));
        }
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curReadLen; j++){
            if(j == 0)
                token = strtok(lineBuf, " ");
            else
                token = strtok(NULL, " ");
            curRead._c.push_back((char)atoi(token));
        }
        in->reads.push_back(curRead);
    }
    ifs.getline(lineBuf, 1024);
    for(int i = 0; i < numHaplotypes; ++i){
        int curHapLen;
        Hap curHap;
        ifs.getline(lineBuf, 1024);
        curHapLen = atoi(lineBuf);
        ifs.getline(lineBuf, 1024);
        ifs.getline(lineBuf, 1024);
        for(int j = 0; j < curHapLen; j++){
            curHap.bases.push_back(lineBuf[j]);
        }
        in->haps.push_back(curHap);
    }
    ifs.close();
    return 0;
}

int GetOutputs(pairhmmOutput* out, int outputSize, std::string filename){
    std::ifstream ifs(filename.c_str(), std::ifstream::in);
    if(!ifs.good())
        printf("bad file name %s\n", filename.c_str());
    for(int i = 0; i < outputSize; ++i){
        double ref;
        ifs >> ref;
        union{
            long long i;
            double d;
        } value;
        ifs >> value.i;
        out->likelihoodData.push_back(value.d);
    }
    return 0;
}
/*
void print_input(pairhmmInputDataOpt& input){
    std::ofstream ofs("test.txt", std::ofstream::out);
    ofs << "readListSize " << input.numRead << " numHaplotypes " << input.numHap << "\n";
    for(int i = 0; i < input.numRead; ++i){
        ofs << input.dataPack.readDataLen[i] << "\n";
        ofs << "readDataArray[" << i << "].readBases[" << input.dataPack.readDataLen[i] << "]: " << "\n";
        for(int j = 0; j < input.dataPack.readDataLen[i]; j++){
            ofs << (int)input.dataPack.readData[i].readBases[j] << " ";
        }
        ofs << "\n";
        ofs << "readDataArray[" << i << "].readQuals[" << input.dataPack.readDataLen[i] << "]: " << "\n";
        for(int j = 0; j < input.dataPack.readDataLen[i]; j++){
            ofs << (int)input.dataPack.readData[i].readQuals[j] << " ";
        }
        ofs << "\n";
        ofs << "readDataArray[" << i << "].insertionGOP[" << input.dataPack.readDataLen[i] << "]: " << "\n";
        for(int j = 0; j < input.dataPack.readDataLen[i]; j++){
            ofs << (int)input.dataPack.readData[i].insertionGOP[j] << " ";
        }
        ofs << "\n";
        ofs << "readDataArray[" << i << "].deletionGOP[" << input.dataPack.readDataLen[i] << "]: " << "\n";
        for(int j = 0; j < input.dataPack.readDataLen[i]; j++){
            ofs << (int)input.dataPack.readData[i].deletionGOP[j] << " ";
        }
        ofs << "\n";
        ofs << "readDataArray[" << i << "].overallGCP[" << input.dataPack.readDataLen[i] << "]: " << "\n";
        for(int j = 0; j < input.dataPack.readDataLen[i]; j++){
            ofs << (int)input.dataPack.readData[i].overallGCP[j] << " ";
        }
        ofs << "\n";
    }
    ofs << "\n";
    for(int i = 0; i < input.numHap; ++i){
        ofs << input.dataPack.hapDataLen[i] << "\n";
        ofs << "mHaplotypeDataArray[" << i << "].haplotypeBases[" << input.dataPack.hapDataLen[i] << "]: \n";
        for(int j = 0; j < input.dataPack.hapDataLen[i]; j++){
            ofs << input.dataPack.hapData[i].haplotypeBases[j];
        }
        ofs << "\n";
    }
    ofs.close();
}*/

int cmp(pairhmmOutput* target, pairhmmOutput* golden, int size, int test_id, bool exact_match, double& total_error_count, double& largest_error){
    int error_count = 0;
    for(int i = 0; i < size; i++){
        if(exact_match){
            if(target->likelihoodData[i] != golden->likelihoodData[i]){ 
                printf("errors in %d th data of %d test, target is %f, golden is %f\n", i, test_id, target->likelihoodData[i], golden->likelihoodData[i]);
                error_count++;
            }
        }
        else{
            if(std::isnan(target->likelihoodData[i])){
                printf("error, target is nan\n");
                error_count++;
            }
            double cur_error = fabs((target->likelihoodData[i] - golden->likelihoodData[i]) / golden->likelihoodData[i]);
            if(cur_error > largest_error){
                largest_error = cur_error;
            }
            if(cur_error > 5e-3){ 
                printf("%dth test: %dth result has significant error, golden=%f, target=%f\n",test_id, i, golden->likelihoodData[i], target->likelihoodData[i]);
                error_count++;
            }
        }
    }
    if(error_count > 0){
        printf("%d out of %d have significant error\n", error_count, size);
    }

    total_error_count = total_error_count + error_count;
    return 0;
}
/*
void analyzeInput(pairhmmInputDataOpt* input, float* read_distrbtn, float* avg_read_length, float* read_count){
    if(worthFPGA(input)){
        for(int i = 0; i < input->numRead; i++){
            read_distrbtn[input->dataPack.readDataLen[i]] = read_distrbtn[input->dataPack.readDataLen[i]] + 1;
            *avg_read_length = ((*read_count) * (*avg_read_length) + input->dataPack.readDataLen[i]) / ((*read_count) + 1);
            *read_count = *read_count + 1;
        }
    }
}*/


double countCells(pairhmmInput* input){    
    double numCell = 0;
    for(int i = 0; (size_t)i < input->reads.size(); i++){
        for(int j = 0; (size_t)j < input->haps.size(); j++){
            double cur_numCell = input->reads[i].bases.size() * input->haps[j].bases.size();
            numCell += cur_numCell;
        }
    }
    return numCell;
}


void printHelp(){
    printf("./test_bin -v or ./test_bin --version: get the compatible platform for current host binary\n");
    printf("./test_bin -h or ./test_bin --help: print the help information\n");
    printf("./test_bin [bitstream filename] --real [real cases folder]: run real tests\n"); 
    printf("./test_bin [bitstream filename] --syn [syn cases number]: run synthetic tests\n");
}

void printVersion(){
    printf("Host code for Xilinx %s \n", PLATFORM);
    printf("DSA version is %s \n", DSA);
    printf("slr0 has %d PEs, slr1 has %d PEs, slr2 has %d PEs\n", SLR0_PE_NUM, SLR1_PE_NUM, SLR2_PE_NUM);

}



int main(int argc, char* argv[]){
    if(argc == 2) {
        if(strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
            printVersion();
            return 0;
        }
        else if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
            printHelp();
            return 0;
        }
        else{
            printf("Invalid argument list\n");
            printHelp();
            return EXIT_FAILURE;
        }
    }

    if(argc != 4){
        printf("Invalid argument list\n");
        printHelp();
        return EXIT_FAILURE;
    }
    
    std::string input_common_prefix;
    std::string output_common_prefix;
    DIR *dir;
    struct dirent *ent;

    bool synthetic = false;
    int test_num = 0;
#ifdef HLS_SIM
        input_common_prefix = std::string("/curr/jysheng/pairHmm_testcases_WGS/input");
        output_common_prefix = std::string("/curr/jysheng/pairHmm_testcases_WGS/output");
#else
    if(strcmp(argv[2], "--real") == 0){
        if ((dir = opendir (argv[3])) != NULL) {
            while ((ent = readdir (dir)) != NULL) {
                test_num++;
            }
            closedir (dir);
        } else {
            printf("cannot find dir %s\n", argv[3]);
            return EXIT_FAILURE;
        }
        printf("find %d files in dir %s\n", test_num, argv[3]);

        input_common_prefix = std::string(argv[3]) + std::string("input");
        output_common_prefix = std::string(argv[3]) + std::string("output");
        test_num = (test_num - 2) / 2;
    }
    else if(strcmp(argv[2], "--syn") == 0){
        try{
            test_num = std::stoi(std::string(argv[3]));
        }
        catch(const std::invalid_argument& ia){
            std::cout << "Invalid synthetic cases number " << ia.what() << '\n';
            return EXIT_FAILURE;
        }
        synthetic = true;
    }
    else{
        printHelp();
        return EXIT_FAILURE;
    }
#endif
    
    //int totalTestNum = 123300;
#ifdef HLS_SIM 
    int totalTestNum = 1;
#else
#ifdef SIM
    int totalTestNum = test_num;
#else
    int totalTestNum = test_num;
#endif
#endif

    avg_numReads = 0;
    avg_numHaps = 0;
    int total_numHaps = 0;
    case_counter = 0;
    pairhmmInput* input;
    pairhmmOutput* golden_output;
    pairhmmOutput* target_output;
    input = new pairhmmInput();
    golden_output = new pairhmmOutput();
    target_output = new pairhmmOutput();
    struct timespec time1, time2, time_diff;
    double total_avx_count = 0;
    double total_fpga_count = 0;
    double total_avx_time = 0;
    double total_fpga_time = 0;
    double total_kernel_time = 0;
    double cur_time = 0;
    double total_avx_cells = 0;
    double total_fpga_cells = 0;
    double total_fpga_results = 0;
    double total_avx_results = 0;
    double current_cells = 0;
    double error_count = 0;
    double largest_error = 0;
    float peak_GCUPS = 0;
#ifdef HLS_SIM
    FalconPairHMM* falcon = new FalconPairHMM(NULL);
#else
#ifdef FPGA
    FalconPairHMM* falcon = new FalconPairHMM(argv[1]);
#else
    FalconPairHMM* falcon = new FalconPairHMM();
#endif
#endif

#ifdef HLS_SIM
    for(int i = 0; i < totalTestNum; ++i){
#else
        //495
    for(int i = 0; i < totalTestNum; ++i){
    //for(int i = 0; i < totalTestNum; ){
#endif
        if(synthetic == false){
            std::string test_id = std::to_string(i);
            std::string input_filename = input_common_prefix + test_id;
            GetInputs(input, input_filename);
            std::string output_filename = output_common_prefix + test_id;
            GetOutputs(golden_output, input->reads.size() * input->haps.size(), output_filename);
        }
        else{
            GenInputs(input, i);
            GenOutputs(input, golden_output);
        }


        bool usedFPGA = false;
        current_cells = countCells(input);
        clock_gettime(CLOCK_REALTIME, &time1);
        falcon->computePairhmm(input, target_output, usedFPGA);
        clock_gettime(CLOCK_REALTIME, &time2);
        time_diff = diff_time(time1, time2);
        cur_time = (long)(time_diff.tv_sec*1e9 + time_diff.tv_nsec);
        if(current_cells / cur_time > peak_GCUPS)
            peak_GCUPS = current_cells / cur_time;
        if(usedFPGA){
            total_fpga_time += cur_time;
            total_fpga_cells += current_cells;
            total_fpga_results += input->reads.size() * input->haps.size();
            total_fpga_count += 1;
            total_numHaps += input->haps.size();
        }
        else{
            total_avx_time += cur_time;   
            total_avx_cells += current_cells;
            total_avx_results += input->reads.size() * input->haps.size();
            total_avx_count += 1;
        }
        cmp(target_output, golden_output, input->reads.size() * input->haps.size(), i, false, error_count, largest_error);
        printf("%d test passed, overall time is %f secs, overall GCUPS is %f, current GCUPS is %f, peak GCUPS is %f\n", i, (total_avx_time + total_fpga_time) * (1e-9), (total_fpga_cells + total_avx_cells) / (total_fpga_time + total_avx_time), current_cells / cur_time, peak_GCUPS);
        input->reads.clear();
        input->haps.clear();
        golden_output->likelihoodData.clear();
        target_output->likelihoodData.clear();
    }
#ifdef FPGA
    total_kernel_time = falcon->get_kernel_time();
    delete falcon;
    printf("%f out of %d tests use FPGA, use %f secs, FPGA GCUPs is %f, avg numHaps = %f\n", total_fpga_count, totalTestNum, total_fpga_time * 1e-9, total_fpga_cells / total_fpga_time, (float)total_numHaps /(float)total_fpga_count); 
    printf("pure kernel GCUPs is %f\n", total_fpga_cells / total_kernel_time); 
#endif
    delete input;
    delete golden_output;
    delete target_output;
    printf("%f out of %d tests use AVX, use %f secs, AVX GCUPs is %f\n", total_avx_count, totalTestNum, total_avx_time * 1e-9, total_avx_cells / total_avx_time);
    printf("%e out of %e FPGA run results have significant errors, total number of results is %e, error rate is %e, largest error is %e\n", error_count, total_fpga_results, (total_fpga_results + total_avx_results), error_count / total_fpga_results, largest_error);
    return 0;
}
