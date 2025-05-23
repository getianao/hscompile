#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "hs.h"
#include "hs_compile_mnrl.h"
#include "ht.h"
#include "read_input.h"
#include <time.h>
#include <unistd.h>
#include <stdint.h>


uint64_t read_rapl_energy() {
    FILE *file = fopen("/sys/class/powercap/intel-rapl/intel-rapl:0/energy_uj", "r");
    if (!file) {
        perror("Failed to open RAPL energy file");
        return 0;
    }
    uint64_t energy;
    fscanf(file, "%lu", &energy);
    fclose(file);
    return energy;
}
typedef struct run_ctx_t {
    r_map *report_map;
    hs_database_t *database;
    hs_scratch_t *scratch;
    char *inputData;
    size_t length;
    unsigned int inp_off; //for the location in filename list
    unsigned int db_off; // for the location in filename list
    unsigned int *counts; // for support values, if used
} run_ctx;

/**
 * This is the function that will be called for each match that occurs. @a ctx
 * is to allow you to have some application-specific state that you will get
 * access to for each match. In our simple example we're just going to use it
 * to pass in the pattern that was being searched for so we can print it out.
 */
static int eventHandler(unsigned int id, unsigned long long from,
                        unsigned long long to, unsigned int flags, void *ctx) {
    
    r_map *report_map = (r_map *) ctx;
    
    r_map *m = find_mapping(id, &report_map);
    if(m == NULL) {
        printf("couldn't find mapping: %u\n", id);
        return 1;
    }    
    printf("Match at id::code::offset %s::%s::%llu\n", m->name, m->report, to);
    return 0;
}

static int supportEventHandler(unsigned int id, unsigned long long from,
                               unsigned long long to, unsigned int flags, void *ctx) {
    
    unsigned int *support = (unsigned int *) ctx;
    
    #pragma omp atomic
    support[id]++;
    
    return 0;
}


static void usage(char *prog) {
    fprintf(stderr, "Usage: %s [-t NUM_TREADS] [--support] <hs databases> <input files>\n", prog);
    fprintf(stderr, "     -t NUM_THREADS     use no more than NUM_THREADS threads\n");
    fprintf(stderr, "     -d          duplicate input stream\n");
    fprintf(stderr, "     -i          input length\n");
    fprintf(stderr, "     -v          quick validation\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        usage(argv[0]);
		return 1;
    }

    printf("Command: ");
    for (int i = 0; i < argc; i++)
        printf("%s ", argv[i]);
    printf("\n");
    
    unsigned int num_dbs = 0;
    unsigned int num_inputs = 0;
    
    bool support = false;
    
    int num_threads = 0;
    int input_len = 2000000;
    int duplicate_input_stream = 1;
    int quick_validation = -1;

    char *db_fns[argc];
    char *input_fns[argc];
    
    for ( int i=1; i<argc; i++ ) {
        
        if ( strcmp("-t", argv[i]) == 0 ) {
            // setting number of threads
            
            if ( i+1 <= argc ) {
                i++;
                num_threads = atoi(argv[i]);
            } else {
                usage(argv[0]);
                return 44;  
            }
            continue;
        }

        if ( strcmp("-i", argv[i]) == 0 ) {
            
            if ( i+1 <= argc ) {
                i++;
                input_len = atoi(argv[i]);
            } else {
                usage(argv[0]);
                return 44;  
            }
            continue;
        }

        if ( strcmp("-d", argv[i]) == 0 ) {
            
            if ( i+1 <= argc ) {
                i++;
                duplicate_input_stream = atoi(argv[i]);
            } else {
                usage(argv[0]);
                return 44;  
            }
            continue;
        }

        if ( strcmp("-v", argv[i]) == 0 ) {
            
            if ( i+1 <= argc ) {
                i++;
                quick_validation = atoi(argv[i]);
            } else {
                usage(argv[0]);
                return 44;  
            }
            continue;
        }
        
        if ( strcmp("--support", argv[i]) == 0 ) {
            // turn on support reporting
            
            support = true;
            
            continue;
        }
        
        size_t len = strlen(argv[i]);
        if ( len < 4 ) {
            input_fns[num_inputs] = argv[i];
            num_inputs += 1;
        } else {
            if ( argv[i][len-1] == 's' && argv[i][len-2] == 'h' && argv[i][len-3] == '.' )  {
                db_fns[num_dbs] = argv[i];
                num_dbs += 1;
            } else {
                input_fns[num_inputs] = argv[i];
                num_inputs += 1;
            }
        }
    }
    
    if ( num_dbs == 0 || num_inputs == 0 ) {
        usage(argv[0]);
        return 45;
    }
    bool isDup = false;
    if (duplicate_input_stream > 1 && num_inputs == 1) {
        num_inputs = duplicate_input_stream;
        isDup = true;
    }
    
    run_ctx contexts[num_dbs*num_inputs];
    
    size_t inputs_length[num_inputs];
    printf("num_inputs = %u\n", num_inputs);
    printf("num_dbs = %u\n", num_dbs);


    // for cleanup
    hs_database_t *dbs_to_delete[num_dbs];
    char *inputs_to_delete[num_inputs];
    r_map *rmaps_to_delete[num_dbs];
    unsigned int *supports_to_delete[num_dbs];
    
    //loop through the inputs to get input data
    if (isDup) {
        printf("duplicate_input_stream = %d\n", duplicate_input_stream);
        for (int j = 0; j < num_inputs; j++) {
            char *inputFN = input_fns[0];
            size_t length;

            /* Next, we read the input data file into a buffer. */
            char *inputData;
            inputData = readInputData(inputFN, &length, input_len);
            if (!inputData) {
                fprintf(stderr,
                        "ERROR: Unable to read input data '%s'. Exiting.\n",
                        inputFN);
                // hs_free_database(database);
                for (int i = 0; i < j; i++) {
                  free(inputs_to_delete[i]);
                }
                return 4;
            }

            inputs_to_delete[j] = inputData;
            inputs_length[j] = length;
            printf("inputs_length_%d = %lu\n", j, inputs_length[j]);
        }
    } else {
        for (int j = 0; j < num_inputs; j++) {
            char *inputFN = input_fns[j];
            size_t length;

            /* Next, we read the input data file into a buffer. */
            char *inputData;
            inputData = readInputData(inputFN, &length, input_len);
            if (!inputData) {
                fprintf(stderr,
                        "ERROR: Unable to read input data '%s'. Exiting.\n",
                        inputFN);
                // hs_free_database(database);
                for (int i = 0; i < j; i++) {
                  free(inputs_to_delete[i]);
                }
                return 4;
            }

            inputs_to_delete[j] = inputData;
            inputs_length[j] = length;
            printf("inputs_length_%d = %lu\n", j, inputs_length[j]);
        }
    }

    // loop through the dbs
    for ( int i=0; i < num_dbs; i++ ) {
        
        char *hsDB = db_fns[i];
        
        // First, read in the database
        size_t length;
        char *hsSerDB;
        char *hsmSer;
        
        hsmSer = readInputData(hsDB, &length, UINT_MAX);
        
        if(!hsmSer) {
            return 2;
        }
        
        // extract the mapping
        r_map *report_map = NULL;
        
        size_t map_length;
        unserialize_mapping (hsmSer, &map_length, &report_map);
        
        rmaps_to_delete[i] = report_map;
        
        // redo the database pointer 
        hsSerDB = hsmSer + map_length ;
        length -= map_length ;
        
        // Next, we try to deserialize
        hs_database_t *database;
        hs_compile_error_t *compile_err;
    
        if(hs_deserialize_database(hsSerDB, length, &database) != HS_SUCCESS) {
            fprintf(stderr, "ERROR: Unable to load HyperScan database file \"%s\": %s. Exiting.\n",
                    hsDB, compile_err->message);
            free(hsmSer);
            
            delete_all(&report_map);
            hs_free_compile_error(compile_err);
            
            // deallocate inputs
            for(int j=0; j<num_inputs; j++) {
                free(inputs_to_delete[j]);
            }
            
            
            for(int j=0; j<i; j++) {
                // deallocate previous databases
                free(dbs_to_delete[j]);
                
                if(support){
                    free(supports_to_delete[j]);
                }
                
                // delete report map
                delete_all(&(rmaps_to_delete[j]));
                
                // kill off all the scratch space that was stored previously
                for(int k=0; k<num_inputs; k++) {
                    free(contexts[j*num_inputs+k].scratch);
                }
            }
            
            return 3;
        }
        
        // keep track of the database
        dbs_to_delete[i] = database;
        
        // make a support array if needed
        if(support) {
            supports_to_delete[i] = (unsigned int*) malloc(sizeof(unsigned int) * count_mapping(&report_map));
            memset(supports_to_delete[i], 0, sizeof(unsigned int) * count_mapping(&report_map));
        } else {
            supports_to_delete[i] = NULL;
        }
        
        //printf("Allocating scratch...\n");
        hs_scratch_t *db_scratch = NULL;
        if (hs_alloc_scratch(database, &db_scratch) != HS_SUCCESS) {
            fprintf(stderr, "ERROR: Unable to allocate scratch space for database '%s'. Exiting.\n", hsDB);
            free(hsmSer);
            delete_all(&report_map);
            hs_free_database(database);
            
            // deallocate inputs
            for(int j=0; j<num_inputs; j++) {
                free(inputs_to_delete[j]);
            }
            
            if(support) {
                free(supports_to_delete[i]);
            }
            
            
            for(int j=0; j<i; j++) {
                // deallocate previous databases
                free(dbs_to_delete[j]);
                
                if(support) {
                    free(supports_to_delete[j]);
                }
                
                // delete report map
                delete_all(&(rmaps_to_delete[j]));
                
                // kill off all the scratch space that was stored previously
                for(int k=0; k<num_inputs; k++) {
                    free(contexts[j*num_inputs+k].scratch);
                }
            }
            
            return 5;
        }
        
        
        //loop through the inputs
        for ( int j=0; j < num_inputs; j++ ) {
            char *inputData;
            inputData = inputs_to_delete[j];
            
            /* Finally, we issue a call to hs_scan, which will search the input buffer
             * for the pattern represented in the bytecode. Note that in order to do
             * this, scratch space needs to be allocated with the hs_alloc_scratch
             * function. In typical usage, you would reuse this scratch space for many
             * calls to hs_scan, but as we're only doing one, we'll be allocating it
             * and deallocating it as soon as our matching is done.
             *
             * When matches occur, the specified callback function (eventHandler in
             * this file) will be called. Note that although it is reminiscent of
             * asynchronous APIs, Hyperscan operates synchronously: all matches will be
             * found, and all callbacks issued, *before* hs_scan returns.
             *
             * In this example, we provide the input pattern as the context pointer so
             * that the callback is able to print out the pattern that matched on each
             * match event.
             */
            
            hs_scratch_t *scratch = NULL;
            if ( hs_clone_scratch(db_scratch, &scratch) != HS_SUCCESS ) {
                printf("ERROR: Unable to allocate cloned scratch space. Exiting.\n");
                free(inputData);
                free(hsmSer);
                delete_all(&report_map);
                hs_free_database(database);
                
                // deallocate inputs
                for(int j=0; j<num_inputs; j++) {
                    free(inputs_to_delete[j]);
                }
                
                if(support) {
                    free(supports_to_delete[i]);
                }
                
                
                for(int j=0; j<i; j++) {
                    // deallocate previous databases
                    free(dbs_to_delete[j]);
                    
                    if(support) {
                        free(supports_to_delete[j]);
                    }
                    
                    // delete report map
                    delete_all(&(rmaps_to_delete[j]));
                    
                    // kill off all the scratch space that was stored previously
                    for(int k=0; k<num_inputs; k++) {
                        free(contexts[j*num_inputs+k].scratch);
                    }
                }
                return 7;
            }
            
            /* Store all of the context information */
            contexts[i*num_inputs+j].report_map = report_map;
            contexts[i*num_inputs+j].database = database;
            contexts[i*num_inputs+j].scratch = scratch;
            contexts[i*num_inputs+j].inputData = inputData;
            contexts[i*num_inputs+j].length = inputs_length[j];
            contexts[i*num_inputs+j].db_off = i;
            contexts[i*num_inputs+j].inp_off = j;
            contexts[i*num_inputs+j].counts = supports_to_delete[i];
            
            
            /* Scanning is complete, any matches have been handled, so now we just
             * clean up and exit.
             */
        } // input loop
        free(hsmSer);
        // hs_free_database(database);
        
        // free up db_scratch
        hs_free_scratch(db_scratch);
    } // database loop
 
    //printf("Simulating graph on input data with Hyperscan...\n");
    
    //okay do the scanning
    if(num_threads > 0) {
        omp_set_dynamic(1);
        omp_set_num_threads(num_threads);
    }

    struct timespec start, end;
    uint64_t start_energy = read_rapl_energy();
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    
    #pragma omp parallel for
    for ( int i=0; i<num_inputs*num_dbs; i++ ) {
        run_ctx ctx = contexts[i];
        
        if(!support) {
            // scan each input and report runtime
            if (hs_scan(ctx.database, ctx.inputData, ctx.length, 0, ctx.scratch, eventHandler,
                        ctx.report_map) != HS_SUCCESS) {
                fprintf(stderr, "ERROR: Unable to scan input buffer '%s' with database '%s'.\n", input_fns[ctx.inp_off], db_fns[ctx.db_off]);
                /*
                 * No need to stop, just keep trying
                hs_free_scratch(scratch);
                free(inputData);
                free(hsmSer);
                hs_free_database(database);
                return 6;
                */
            }
        } else {
            if (hs_scan(ctx.database, ctx.inputData, ctx.length, 0, ctx.scratch, supportEventHandler,
                        ctx.counts) != HS_SUCCESS) {
                fprintf(stderr, "ERROR: Unable to scan input buffer '%s' with database '%s'.\n", input_fns[ctx.inp_off], db_fns[ctx.db_off]);
            }
        }
        
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t end_energy = read_rapl_energy();
    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    
    printf("delta_us = %lu\n", delta_us);
    double total_bytes = 0;
    for (int i = 0; i < num_inputs; i++) {
        total_bytes += inputs_length[i];
    }

    double second = (0.0+delta_us) / 1000000;
    printf("time_second = %.2lf\n", second);

    double throughput_bytepersecond = total_bytes / second;
    double throughput_MBs = throughput_bytepersecond / 1000.0 / 1000.0;
    printf("%.2lf MB/s\n", throughput_MBs);

    // print out supports
    if(support) {
        unsigned long long int total = 0;
        printf("File, ID, Report ID, Count\n");
        for ( int i=0; i<num_dbs; i++ ) {
            r_map *mapping = rmaps_to_delete[i];
            for ( int j=0; j<count_mapping(&mapping); j++) {
                r_map *m = find_mapping(j, &mapping);
                printf("%s, %s, %s, %u\n", db_fns[i], m->name, m->report, supports_to_delete[i][j]);
                total+=supports_to_delete[i][j];
            }
        }
        printf("total_count = %llu\n", total);
        if (quick_validation >= 0) {
            if (total >=
                ((long long int)(quick_validation)) * duplicate_input_stream) {
                if (total == ((long long int)(quick_validation)) *
                                 duplicate_input_stream)
                    printf("Quick Validation PASS! (perfect)\n");
                else
                    printf("Quick Validation PASS! (not perfect)\n");
            } else {
                printf("Quick Validation FAILED! \n");
            }
        } else {
            printf("Fake Validation PASS! \n");
        }
    } else {
        printf("Fake Validation PASS! \n");
    }
    printf("throughput = %lf MB/s\n", throughput_MBs);
    double power = (end_energy - start_energy) / 1000000.0;
    printf("Power: %lf J\n", power);
    double watts = power / second;
    printf("Watts: %lf W\n", watts);
    double power_efficiency = total_bytes / 1000000.0 / power;
    printf("Power_Efficiency: %lf MB/J\n", power_efficiency);
    
    // cleanup
    for ( int i=0; i<num_inputs*num_dbs; i++ ) {
        run_ctx ctx = contexts[i];
        hs_free_scratch(ctx.scratch);
    }
    
    for ( int i=0; i<num_dbs; i++ ) {
        free(dbs_to_delete[i]);
        delete_all(&(rmaps_to_delete[i]));
        free(supports_to_delete[i]);
    }
    
    for ( int i=0; i<num_inputs; i++ ) {
        free(inputs_to_delete[i]);
    }
    printf("FINISHED\n");
    
    return 0;
}
