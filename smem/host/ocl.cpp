#include "host/ocl.h"
#include <stdio.h>
#include <iostream>
#include <vector>


bool init = false;
struct smemOclCtx ctx;

// Checks OpenCL error codes
void check(cl_int err_code) {
    if (err_code != CL_SUCCESS) {
        printf("ERROR: %d\n", err_code);
        exit(EXIT_FAILURE);
    }
}

// An event callback function that prints the operations performed by the OpenCL
// runtime.
void event_cb(cl_event event, cl_int cmd_status, void *data) {
    cl_command_type command;
    clGetEventInfo(event, CL_EVENT_COMMAND_TYPE, sizeof(cl_command_type),
                 &command, nullptr);
    cl_int status;
    clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int),
                 &status, nullptr);
    const char *command_str;
    const char *status_str;
    switch (command) {
        case CL_COMMAND_READ_BUFFER:
            command_str = "buffer read";
        break;
        case CL_COMMAND_WRITE_BUFFER:
            command_str = "buffer write";
        break;
        case CL_COMMAND_NDRANGE_KERNEL:
            command_str = "kernel";
        break;
    }
    switch (status) {
        case CL_QUEUED:
            status_str = "Queued";
        break;
        case CL_SUBMITTED:
            status_str = "Submitted";
        break;
        case CL_RUNNING:
            status_str = "Executing";
        break;
        case CL_COMPLETE:
            status_str = "Completed";
        break;
    }
    printf("%s %s %s\n", status_str, reinterpret_cast<char *>(data), command_str);
    fflush(stdout);
}

// Sets the callback for a particular event
void set_callback(cl_event event, const char *queue_name) {
    clSetEventCallback(event, CL_COMPLETE, event_cb, (void *)queue_name);
}

int load_file_to_memory(const char *filename, char **result) { 
    size_t size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) { 
        printf("ERROR : Kernel binary %s not exist!\n", filename);
        *result = NULL;
        return -1; // -1 means file opening fail 
    } 
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if ((int)size != (int)fread(*result, sizeof(char), size, f)) { 
        free(*result);
        return -2; // -2 means file reading fail 
    } 
    fclose(f);
    (*result)[size] = 0;
    return size;
}

bool init_FPGA(char* bitstream){
    char cl_platform_vendor[1001];
    char cl_platform_name[1001];
    cl_platform_vendor[0] = 0;
    cl_platform_name[0] = 0;
    int err;
    
    err = clGetPlatformIDs(1, &ctx.platform_id, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: Failed to find an OpenCL platform!Code %i\n", err);
        return false;
    }
    printf("Successfully create platform %d\n", ctx.platform_id);
    err = clGetPlatformInfo(ctx.platform_id, CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed! Code %i\n", err);
        return false;
    }
    printf("CL_PLATFORM_VENDOR %s\n", cl_platform_vendor);
    err = clGetPlatformInfo(ctx.platform_id, CL_PLATFORM_NAME, 1000, (void *)cl_platform_name, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: clGetPlatformInfo(CL_PLATFORM_NAME) failed! Code %i\n", err);
        return false;
    }
    printf("CL_PLATFORM_NAME %s\n", cl_platform_name);
    
    err = clGetDeviceIDs(ctx.platform_id, CL_DEVICE_TYPE_ACCELERATOR, 1, &ctx.device_id, NULL);
    if (err != CL_SUCCESS) {
        printf("Warning: Failed to create a device group! Code %i\n", err);
        return false;
    }
    printf("Successfully create device %d\n", ctx.device_id);

    ctx.context = clCreateContext(0, 1, &ctx.device_id, NULL, NULL, &err);
    if (!ctx.context) {
        printf("Warning: Failed to create a compute context! Code %i\n", err);
        return false;
    }
    printf("Successfully create context \n");

    unsigned char *kernelbinary;

    int n_i = 0;
    n_i = load_file_to_memory(bitstream, (char **) &kernelbinary);
    if (n_i < 0) {
        printf("Warning : failed to load kernel from binary: %s\n", bitstream);
        return false;
    }
    printf("Successfully load kernel from binary: %s\n", bitstream);
    
    int status;
    size_t n = n_i;
    ctx.program = clCreateProgramWithBinary(ctx.context, 1, &ctx.device_id, &n, \
            (const unsigned char **) &kernelbinary, &status, &err);
    if ((!ctx.program) || (err!=CL_SUCCESS)) {
        printf("Warning: Failed to create compute program from binary! Code %d\n", err);
        return false;
    }
    printf("Success to create compute program from binary! \n");

    // Build the program executable
    err = clBuildProgram(ctx.program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[2048];
        printf("Warning: Failed to build program executable!\n");
        clGetProgramBuildInfo(ctx.program, ctx.device_id, CL_PROGRAM_BUILD_LOG, \
                sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        return false;
    }
    printf("Sucess to build program executable!\n");
    ctx.command = clCreateCommandQueue(ctx.context, ctx.device_id,
            CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE,
            &err);
    if (err != CL_SUCCESS){
        std::cout << "Error: Failed to create a command queue!" << std::endl;
        std::cout << "Test failed" << std::endl;
        return false;
    }
    return true;

}


void ocl_init(char* btsm, const uint32_t* bwt, const uint64_t* bwt_para, uint64_t bwt_size, bwtintv_t* host_mem, int batch_size){
    int err;
    std::string kernel_comm = "mem_collect_intv_core";
    std::string kernel_dram_comm = "bwt_request_core";
    printf("bwt_size inside = %lu words, %lu 512 bits\n", bwt_size, bwt_para[6]);
    int kernel_max_batch_size = BATCH_SIZE / 4;
    printf("kernel_max_batch_size = %d\n", kernel_max_batch_size); 
    cl_mem_ext_ptr_t bank_ext[BANK_NUM];
#if BANK_NUM==4
    unsigned bankID[BANK_NUM] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
#else
    unsigned bankID[BANK_NUM] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
#endif
    if(!init){
        init_FPGA(btsm);
        for(int i = 0; i < BANK_NUM; i++){
            double total_dram_size = 0;
            bank_ext[i].flags = bankID[i];
            bank_ext[i].param = 0;
            bank_ext[i].obj = nullptr;
#if BANK_NUM==4
            std::string cur_kernel = kernel_comm + std::to_string(i);
#else
            std::string cur_kernel = kernel_comm + std::to_string(i < 1 ? i : i + 1);
#endif
            ctx.kernel_smem[i] = clCreateKernel(ctx.program, cur_kernel.c_str(), &err);
            if (!ctx.kernel_smem[i]) {
                printf("failed to create mem_collect_intv kernel %d\n, err code is %d\n", i, err);
            }
            else 
                printf("create kernel %d\n", i);
#ifdef CL
#if BANK_NUM==4
            cur_kernel = kernel_dram_comm + std::to_string(i);
#else
            cur_kernel = kernel_dram_comm + std::to_string(i < 1 ? i : i + 1);
#endif
            ctx.kernel_dram[i] = clCreateKernel(ctx.program, cur_kernel.c_str(), &err);
            if (!ctx.kernel_dram[i]) {
                printf("failed to create bwt_request kernel %d\n, err code is %d\n", i, err);
            }
            else 
                printf("create kernel %d\n", i);
#endif
            ctx.bwt_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    bwt_size * sizeof(uint32_t), &bank_ext[i], &err);
            if (!ctx.bwt_buffer[i]) {
                printf("failed to create bwt buffer %d\n, err code is %d\n", i, err);
            }
            else{
                printf("create bwt buffer %d\n", i);
                total_dram_size += (double)bwt_size * sizeof(uint32_t); 
            }
            int cl_status; 
            cl_status = clEnqueueWriteBuffer(ctx.command, ctx.bwt_buffer[i], CL_TRUE, 0,\
                    sizeof(uint32_t) * bwt_size, bwt, 0, NULL, NULL);
            if (cl_status != CL_SUCCESS) {
                printf("failed to write bwt buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("write bwt buffer %d\n", i);
            }

            ctx.bwt_para_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    7 * sizeof(uint64_t), &bank_ext[i], &err);
            if (!ctx.bwt_para_buffer[i]) {
                printf("failed to create bwt para buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("create bwt para buffer %d\n", i);
                total_dram_size += 7 * sizeof(uint64_t); 
            }
            cl_status = clEnqueueWriteBuffer(ctx.command, ctx.bwt_para_buffer[i], CL_TRUE,\
                    0, sizeof(uint64_t) * 7, bwt_para, 0, NULL, NULL);
            if (cl_status != CL_SUCCESS) {
                printf("failed to write bwt para buffer %d\n, err code is %d\n", i, err);
            }
            else
                printf("write bwt para buffer %d\n", i);
            ctx.seq_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    kernel_max_batch_size * SEQ_LENGTH * sizeof(uint8_t), &bank_ext[i], &err);
            if (!ctx.seq_buffer[i]) {
                printf("failed to create seq buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("create seq buffer %d\n", i);
                total_dram_size += kernel_max_batch_size * SEQ_LENGTH * sizeof(uint8_t); 
            }
            
            ctx.mem_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    kernel_max_batch_size * MAX_INTV_ALLOC * sizeof(bwtintv_t), &bank_ext[i], &err);
            if (!ctx.mem_buffer[i]) {
                printf("failed to create mem buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("create mem buffer %d, err = %d\n", i, err);
                total_dram_size += (double)kernel_max_batch_size * MAX_INTV_ALLOC * sizeof(bwtintv_t); 
            }

            ctx.seq_len_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    kernel_max_batch_size * sizeof(uint8_t), &bank_ext[i], &err);
            if (!ctx.seq_len_buffer[i]) {
                printf("failed to create seq len buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("create seq len buffer %d\n", i);
                total_dram_size += kernel_max_batch_size * sizeof(uint8_t); 
            }
            ctx.mem_num_buffer[i] = clCreateBuffer(ctx.context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX,\
                    kernel_max_batch_size * sizeof(int), &bank_ext[i], &err);
            if (!ctx.mem_num_buffer[i]) {
                printf("failed to create mem_num buffer %d\n, err code is %d\n", i, err);
            }
            else {
                printf("create mem num buffer %d\n", i);
                total_dram_size += kernel_max_batch_size * sizeof(int); 
            }
            printf("bank %d used %e GB\n", i, total_dram_size * 1e-9);
            printf("Besides bwt reference %f Bytes are used\n", total_dram_size - (double)bwt_size * sizeof(uint32_t));
            //ctx.host_mem_buffer[i] = &host_mem[MAX_INTV_ALLOC * i * kernel_max_batch_size];
            //bank_ext[i].obj = ctx.host_mem_buffer[i];
        }
        init = true;
    }
}


void ocl_kernel_invoke(uint8_t* seq, uint8_t* seq_len, bwtintv_t* mem_output, int* mem_num, int batch_size, double cur_kernel_time[BANK_NUM], uint64_t bwt_size){
    int err;
    int kernel_batch_size = batch_size / BANK_NUM;
    printf("before invoke kernel kernel_batch_size = %d\n", kernel_batch_size);
    int kernel_max_batch_size = BATCH_SIZE / BANK_NUM;
    cl_event write_event[2 * BANK_NUM];
    for(int i = 0; i < BANK_NUM; i++){
        err = clEnqueueWriteBuffer(ctx.command, ctx.seq_buffer[i], CL_TRUE, 0, \
                sizeof(uint8_t) * SEQ_LENGTH * kernel_batch_size, seq + i * SEQ_LENGTH * kernel_batch_size, 0, NULL, &write_event[2 * i]);
        err = clEnqueueWriteBuffer(ctx.command, ctx.seq_len_buffer[i], CL_TRUE, 0, \
                sizeof(uint8_t) * kernel_batch_size, seq_len + i * kernel_batch_size, 0, NULL, &write_event[2 * i + 1]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 0, sizeof(cl_mem), &ctx.bwt_para_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 1, sizeof(cl_mem), &ctx.seq_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 2, sizeof(cl_mem), &ctx.seq_len_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 3, sizeof(cl_mem), &ctx.mem_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 4, sizeof(cl_mem), &ctx.mem_num_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_smem[i], 5, sizeof(int), &kernel_batch_size);
#ifdef CL
        err |= clSetKernelArg(ctx.kernel_dram[i], 0, sizeof(cl_mem), &ctx.bwt_buffer[i]);
        err |= clSetKernelArg(ctx.kernel_dram[i], 1, sizeof(uint64_t), &bwt_size);
#else
        err |= clSetKernelArg(ctx.kernel_smem[i], 6, sizeof(cl_mem), &ctx.bwt_buffer[i]);
#endif
    }
    if(err != CL_SUCCESS) {
        fprintf(stderr, "Error: Failed to set kernel arguments %d\n", err);
        return;
    }
    clWaitForEvents(2 * BANK_NUM, write_event);
    
    cl_event kernel_event[BANK_NUM];
    for(int i = 0; i < BANK_NUM; i++){
        err |= clEnqueueTask(ctx.command, ctx.kernel_smem[i], 0, NULL, &kernel_event[i]);
#ifdef CL
        err |= clEnqueueTask(ctx.command, ctx.kernel_dram[i], 0, NULL, NULL);
#endif
    }
    if(err != CL_SUCCESS) {
        fprintf(stderr, "Error: Failed to enqueue tasks, error = %d\n", err);
        return;
    }
    clWaitForEvents(BANK_NUM, kernel_event);

    for (int i = 0; i < BANK_NUM; i++) {
        cl_ulong start, end;
        clGetEventProfilingInfo(kernel_event[i], CL_PROFILING_COMMAND_START, sizeof(start), &start, NULL);
        clGetEventProfilingInfo(kernel_event[i], CL_PROFILING_COMMAND_END, sizeof(end), &end, NULL);
        cur_kernel_time[i] = (end - start);
    }

    for(int i = 0; i < BANK_NUM; i++){
        err |= clEnqueueReadBuffer(ctx.command, ctx.mem_buffer[i], CL_TRUE, 0,\
                sizeof(bwtintv_t) * MAX_INTV_ALLOC * kernel_batch_size, &mem_output[i * MAX_INTV_ALLOC * kernel_batch_size], 0, NULL, NULL);

        err |= clEnqueueReadBuffer(ctx.command, ctx.mem_num_buffer[i], CL_TRUE, 0,\
                sizeof(int) * kernel_batch_size, &mem_num[i * kernel_batch_size], 0, NULL, NULL);
    }

    //err |= clEnqueueMigrateMemObjects(ctx.command, BANK_NUM, ctx.mem_buffer, CL_MIGRATE_MEM_OBJECT_HOST, 0, NULL, NULL);
    if(err != CL_SUCCESS) {
        fprintf(stderr, "Error: Failed to enqueue read buffers, error = %d\n", err);
        return;
    }
    clFinish(ctx.command);
    for (int i = 0; i < 2 * BANK_NUM; i++) {
        clReleaseEvent(write_event[i]);
    }
}

bool ocl_cmp(char* fname_golden, bwtintv_v& mem, int mem_num[BATCH_SIZE]){
    bool pass = true;
    for (int k = 0; k < BATCH_SIZE; k++) {
        mem.n = mem_num[k];
        FILE* fgolden = fopen(fname_golden, "rb");
        int n = 0;
        bwtint_t golden_x0;
        bwtint_t golden_x1;
        bwtint_t golden_x2;
        bwtint_t golden_info;
        fread(&n, sizeof(int), 1, fgolden);
        if (n != mem.n) {
            fclose(fgolden); printf("In batch %d, there are %d intv in the result but %d in the golden.\n", k, mem.n, n);
            return false;
        }
        else {
            for (int i=0; i<n; i++) {
                fread(&golden_x0, sizeof(bwtint_t), 1, fgolden);
                fread(&golden_x1, sizeof(bwtint_t), 1, fgolden);
                fread(&golden_x2, sizeof(bwtint_t), 1, fgolden);
                fread(&golden_info, sizeof(bwtint_t), 1, fgolden);
                if (golden_x0 != mem.a[k * MAX_INTV_ALLOC + i].x[0]) {
                    printf("x0 is different for intv %d\n", i);
                    fclose(fgolden);
                    return false;
                }
                if (golden_x1 != mem.a[k * MAX_INTV_ALLOC + i].x[1]) {
                    printf("x1 is different for intv %d\n", i);
                    fclose(fgolden);
                    return false;
                }
                if (golden_x2 != mem.a[k * MAX_INTV_ALLOC + i].x[2]) {
                    printf("x2 is different for intv %d\n", i);
                    fclose(fgolden);
                    return false;
                }
                if(golden_info != mem.a[k * MAX_INTV_ALLOC + i].info) {
                    printf("info is different for intv %d\n", i);
                    fclose(fgolden);
                    return false;
                }
            }
            fclose(fgolden);
        }
    }
    return true;
}

int smem_ocl(char* btsm, const uint32_t* bwt, const uint64_t* bwt_para, uint8_t* seq, uint8_t* seq_len, \
        int batch_size, bwtintv_t* mem_output, int* mem_num, double kernel_time[BANK_NUM]){
    //ocl_init(btsm, bwt, bwt_para, batch_size);
    ocl_kernel_invoke(seq, seq_len, mem_output, mem_num, batch_size, kernel_time, bwt_para[6]);
    return 0;
}


/*
int ocl_run(const uint32_t *bwt, char* btsm, char* fname_input, char* fname_golden, uint64_t& time_total) {
    const uint64_t bwt_para[7] = {BWT_PRIMARY, L2[0], L2[1], L2[2], L2[3], L2[4], BWT_SIZE};
    FILE* finput = fopen(fname_input, "rb");
    int len;
    uint8_t *seq;
    fread(&len, sizeof(int), 1, finput);
    if(!(seq = (uint8_t*)malloc(len * sizeof(uint8_t) * BATCH_SIZE))){
        printf("failed to malloc space for seq in ocl_run()\n");
        exit(EXIT_FAILURE);
    };
    fread(seq, sizeof(uint8_t), len, finput);
    fclose(finput);

    for(int i = 1; i < BATCH_SIZE; i++) {
        memcpy(&seq[i * SEQ_LENGTH], seq, len * sizeof(uint8_t));
    }
    bwtintv_v mem; 
    uint8_t seq_len[BATCH_SIZE];
    int mem_num[BATCH_SIZE];
    mem.a = (bwtintv_t*)malloc(sizeof(bwtintv_t)*MAX_INTV_ALLOC*BATCH_SIZE);
    for(int i = 0; i < BATCH_SIZE; i++) {
        seq_len[i] = len;
    }
    double kernel_time[BATCH_SIZE];
    ocl_init(btsm, bwt, bwt_para, mem.a, BATCH_SIZE);
    ocl_kernel_invoke(seq, seq_len, mem.a, mem_num, BATCH_SIZE, kernel_time, BWT_SIZE);
    if (ocl_cmp(fname_golden, mem, mem_num)) {
        return 0; 
    }
    else {
        return -1; 
    }
}*/
