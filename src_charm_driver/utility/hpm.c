/*----------------------------------------------------------*/
/* Simple interface to the hardware performance counters.   */
/* Copied from BG/P /soft/apps/UPC/examples/                */
/* This code, along with timebase.c and timebase.h are NOT  */
/* covered by the OpenAtom license. We make no claims as to */
/* authorship or copyright wrt to these files.              */
/* They have been modified from the BG/P distributed source */
/* to remove MPI so they can be used from Charm++.    The   */
/* revised version is included for distribution convenience */
/* purposes only. IBM included no license with these files  */
/* and the files contained no claimed author or copyplate.  */
/* Contact them for more information, but my guess is they  */
/* were hacked up by an intern and are distributed as is,   */
/* caveat emptor.                                           */
/*----------------------------------------------------------*/
/* C / C++:                                                 */
/*    HPM_Init();          initialize the counters          */
/*    HPM_Start("label");  start counting                   */
/*    HPM_Stop("label");   stop  counting                   */
/*    HPM_Print();  print counter values and labels         */
/*----------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <spi/UPC.h>

#define NUM_COUNTERS 256

#define NUM_CLOCK_X1 72

#define LABEL_LEN  80

void HPM_Init(int);
void HPM_Start(char *,int);
void HPM_Stop(char *,int);
void HPM_Print(int, int);

#include "timebase.h"

int index_from_label(char *);

static int initialized = 0;
static int code_block  = 0;

#define MAX_CODE_BLOCKS 20

static unsigned long long timebase_in[MAX_CODE_BLOCKS];
static unsigned long long timebase_sum[MAX_CODE_BLOCKS];
static unsigned long long counter_in[MAX_CODE_BLOCKS][NUM_COUNTERS];
static unsigned long long counter_sum[MAX_CODE_BLOCKS][NUM_COUNTERS];
static char code_block_label[MAX_CODE_BLOCKS][LABEL_LEN];
static int block_starts[MAX_CODE_BLOCKS];
static int block_stops[MAX_CODE_BLOCKS];

void set_labels(int);

static char label[NUM_COUNTERS][LABEL_LEN];


static int trigger_method = BGP_UPC_CFG_LEVEL_HIGH;

static int counter_mode = 0;

struct CounterStruct {
  int32_t rank;                       // Rank
  int32_t core;                       // Core
  int32_t upc_number;                 // UPC Number
  int32_t number_processes_per_upc;   // Number of processes per UPC unit
  BGP_UPC_Mode_t mode;                // User mode
  int32_t number_of_counters;         // Number of counter values returned
  char location[24];                  // Location
  int64_t elapsed_time;               // Elapsed time
  uint32_t reserved_1;                // Reserved for alignment
  uint32_t reserved_2;                // Reserved for alignment
  int64_t values[256];                // Counter values
} counter_data;



/*==========================*/
/* Initialize the counters. */
/*==========================*/
void HPM_Init(int local_rank)
{
   int i, j;
   int counter_trigger;
   char * ptr;

   if (!initialized) 
   {
       initialized = 1;

       // set the initial cumulative counter values to zero 
       for (j=0; j<MAX_CODE_BLOCKS; j++)
          for (i=0; i<NUM_COUNTERS; i++)
              counter_sum[j][i] = (unsigned long long) 0;

       for (j=0; j<MAX_CODE_BLOCKS; j++) timebase_sum[j] = (unsigned long long) 0;

       // keep track of code block starts and stops 
       for (j=0; j<MAX_CODE_BLOCKS; j++)
       {
           block_starts[j] = 0;
           block_stops[j]  = 0;
       }

       // check to see if a counter mode has been specified
       ptr = getenv("BGP_COUNTER_MODE");

       if (ptr == NULL) counter_mode = 0;
       else             counter_mode = atoi(ptr);

       // force the counter mode to be -1,0,1,2,3
       if (counter_mode != 0 && counter_mode != 1 &&
           counter_mode != 2 && counter_mode != 3) counter_mode = 0;

       ptr = getenv("BGP_COUNTER_TRIGGER");

       if (ptr == NULL) counter_trigger = BGP_UPC_CFG_LEVEL_HIGH;
       else
       {
          if (strncasecmp(ptr,"edge", 4) == 0)  counter_trigger = BGP_UPC_CFG_EDGE_RISE;

          if (strncasecmp(ptr,"high", 4) == 0)  counter_trigger = BGP_UPC_CFG_LEVEL_HIGH;
       }

       trigger_method = counter_trigger;

       // every rank calls BGP_UPC_Initialize()
       BGP_UPC_Initialize();

       // only one rank sets the config and zeroes the counters
       if (local_rank == 0)
       {
//	   BGP_UPC_Initialize();
          BGP_UPC_Stop();

          if (counter_mode == 0) BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_0, counter_trigger);
          if (counter_mode == 1) BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_1, counter_trigger);
          if (counter_mode == 2) BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_2, counter_trigger);
          if (counter_mode == 3) BGP_UPC_Initialize_Counter_Config(BGP_UPC_MODE_3, counter_trigger);

          BGP_UPC_Zero_Counter_Values();

          BGP_UPC_Start(0);
       }
   }
   return;
}

/*=============================================*/
/* Start counting for a particular code block. */
/*=============================================*/
void HPM_Start(char * this_label, int local_rank)
{
   int i, j;
   unsigned long long tb;
   if( local_rank == 0){
       tb = timebase();

       BGP_UPC_Stop();
       BGP_UPC_Read_Counter_Values(&counter_data, sizeof(struct CounterStruct), BGP_UPC_READ_EXCLUSIVE);

       j = index_from_label(this_label);

       block_starts[j] += 1;

       timebase_in[j] = tb;

       for (i=0; i<NUM_COUNTERS; i++) counter_in[j][i] = counter_data.values[i];

       BGP_UPC_Start(0);
   }
   return;
}



/*============================================*/
/* Stop counting for a particular code block. */
/*============================================*/
void HPM_Stop(char * this_label, int local_rank)
{
    int i, j;
    unsigned long long tb;

    if (code_block >= MAX_CODE_BLOCKS) return;
    if (local_rank == 0)
    {
	tb = timebase();

	BGP_UPC_Stop();
	BGP_UPC_Read_Counter_Values(&counter_data, sizeof(struct CounterStruct), BGP_UPC_READ_EXCLUSIVE);

	j = index_from_label(this_label);

	block_stops[j] += 1;

	timebase_sum[j] += tb - timebase_in[j];

	for (i=0; i<NUM_COUNTERS; i++) counter_sum[j][i] += (counter_data.values[i] - counter_in[j][i]);

	BGP_UPC_Start(0);
    }
    return;
}



/*============================================*/
/* Print the counter values with event labels */
/*============================================*/
void HPM_Print(int node_id, int local_rank)
{
    int i, j, nblocks;
    int tx, ty, tz;
    unsigned long long counts;
    char trigger[16];
    char filename[132];
    FILE * fp;

    _BGP_Personality_t personality;
    if(local_rank == 0)
    {
    set_labels(counter_mode);

    Kernel_GetPersonality(&personality, sizeof(personality));

    tx = personality.Network_Config.Xcoord;
    ty = personality.Network_Config.Ycoord;
    tz = personality.Network_Config.Zcoord;

    sprintf(filename, "hpm_data.%d", node_id);
    fp = fopen(filename, "w");
    if (fp == NULL) fp = stderr;

    if (trigger_method == BGP_UPC_CFG_EDGE_RISE)  sprintf(trigger, "edge rise");
    if (trigger_method == BGP_UPC_CFG_LEVEL_HIGH) sprintf(trigger, "level high");

    fprintf(fp, "\n");
    fprintf(fp, "----------------------------------------------------------------------\n");
    fprintf(fp, "Hardware counter report for BGP node %d, coordinates <%d,%d,%d>.\n", node_id, tx, ty, tz);
    fprintf(fp, "BGP counter mode = %d, trigger = %s.\n", counter_mode, trigger);
    fprintf(fp, "----------------------------------------------------------------------\n");
    if (code_block >= MAX_CODE_BLOCKS) nblocks = MAX_CODE_BLOCKS;
    else                               nblocks = code_block;

    for (j=0; j<nblocks; j++)
    { 
	if (block_starts[j] == block_stops[j])
	{
	    fprintf(fp, "%s, call count = %d, cycles = %lld :\n", code_block_label[j], block_starts[j], timebase_sum[j]);
	    for (i=0; i<NUM_COUNTERS; i++)
	    {
		if (i<NUM_CLOCK_X1)  counts = counter_sum[j][i];
		else                 counts = counter_sum[j][i] / 2LL;
		fprintf(fp, "%3d %14lld  %s\n", i, counts, label[i]);
	    }
	    fprintf(fp, "\n");
	}
	else
	{
	    fprintf(fp, "mismatch in starts/stops for code block '%s'\n", code_block_label[j]);
	    fprintf(fp, "  starts = %d\n", block_starts[j]);
	    fprintf(fp, "  stops  = %d\n", block_stops[j]);
	}
    }
    fprintf(fp, "\n");
    }
    return;
}


/*===========================================*/
/* Find the code-block number from the label.*/
/*===========================================*/
int index_from_label(char * this_label)
{
   int i, match;
   char * ptr;

   if (code_block < MAX_CODE_BLOCKS)
   {
       match = 0;
       for (i=code_block-1; i>=0; i--)
       {
           if (0 == strcmp(code_block_label[i], this_label))
           {
               match = 1;
               break;
           }
       }
    
       if (match == 0)
       {
           i = code_block;
           ptr = strcpy(code_block_label[i], this_label);
           if (ptr == NULL) code_block_label[i][0] = '\0';
           code_block ++;
       }
   }

   return i;

}

/*====================================================*/
/* routine to set labels for each of the BGP counters */
/*====================================================*/
void set_labels(int mode)
{
  
   if (mode == 0)
   {
      strcpy(label[0], "BGP_PU0_JPIPE_INSTRUCTIONS");
      strcpy(label[1], "BGP_PU0_JPIPE_ADD_SUB");
      strcpy(label[2], "BGP_PU0_JPIPE_LOGICAL_OPS");
      strcpy(label[3], "BGP_PU0_JPIPE_SHROTMK");
      strcpy(label[4], "BGP_PU0_IPIPE_INSTRUCTIONS");
      strcpy(label[5], "BGP_PU0_IPIPE_MULT_DIV");
      strcpy(label[6], "BGP_PU0_IPIPE_ADD_SUB");
      strcpy(label[7], "BGP_PU0_IPIPE_LOGICAL_OPS");
      strcpy(label[8], "BGP_PU0_IPIPE_SHROTMK");
      strcpy(label[9], "BGP_PU0_IPIPE_BRANCHES");
      strcpy(label[10], "BGP_PU0_IPIPE_TLB_OPS");
      strcpy(label[11], "BGP_PU0_IPIPE_PROCESS_CONTROL");
      strcpy(label[12], "BGP_PU0_IPIPE_OTHER");
      strcpy(label[13], "BGP_PU0_DCACHE_LINEFILLINPROG");
      strcpy(label[14], "BGP_PU0_ICACHE_LINEFILLINPROG");
      strcpy(label[15], "BGP_PU0_DCACHE_MISS");
      strcpy(label[16], "BGP_PU0_DCACHE_HIT");
      strcpy(label[17], "BGP_PU0_DATA_LOADS");
      strcpy(label[18], "BGP_PU0_DATA_STORES");
      strcpy(label[19], "BGP_PU0_DCACHE_OPS");
      strcpy(label[20], "BGP_PU0_ICACHE_MISS");
      strcpy(label[21], "BGP_PU0_ICACHE_HIT");
      strcpy(label[22], "BGP_PU0_FPU_ADD_SUB_1");
      strcpy(label[23], "BGP_PU0_FPU_MULT_1");
      strcpy(label[24], "BGP_PU0_FPU_FMA_2");
      strcpy(label[25], "BGP_PU0_FPU_DIV_1");
      strcpy(label[26], "BGP_PU0_FPU_OTHER_NON_STORAGE_OPS");
      strcpy(label[27], "BGP_PU0_FPU_ADD_SUB_2");
      strcpy(label[28], "BGP_PU0_FPU_MULT_2");
      strcpy(label[29], "BGP_PU0_FPU_FMA_4");
      strcpy(label[30], "BGP_PU0_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS");
      strcpy(label[31], "BGP_PU0_FPU_QUADWORD_LOADS");
      strcpy(label[32], "BGP_PU0_FPU_OTHER_LOADS");
      strcpy(label[33], "BGP_PU0_FPU_QUADWORD_STORES");
      strcpy(label[34], "BGP_PU0_FPU_OTHER_STORES");
      strcpy(label[35], "BGP_PU1_JPIPE_INSTRUCTIONS");
      strcpy(label[36], "BGP_PU1_JPIPE_ADD_SUB");
      strcpy(label[37], "BGP_PU1_JPIPE_LOGICAL_OPS");
      strcpy(label[38], "BGP_PU1_JPIPE_SHROTMK");
      strcpy(label[39], "BGP_PU1_IPIPE_INSTRUCTIONS");
      strcpy(label[40], "BGP_PU1_IPIPE_MULT_DIV");
      strcpy(label[41], "BGP_PU1_IPIPE_ADD_SUB");
      strcpy(label[42], "BGP_PU1_IPIPE_LOGICAL_OPS");
      strcpy(label[43], "BGP_PU1_IPIPE_SHROTMK");
      strcpy(label[44], "BGP_PU1_IPIPE_BRANCHES");
      strcpy(label[45], "BGP_PU1_IPIPE_TLB_OPS");
      strcpy(label[46], "BGP_PU1_IPIPE_PROCESS_CONTROL");
      strcpy(label[47], "BGP_PU1_IPIPE_OTHER");
      strcpy(label[48], "BGP_PU1_DCACHE_LINEFILLINPROG");
      strcpy(label[49], "BGP_PU1_ICACHE_LINEFILLINPROG");
      strcpy(label[50], "BGP_PU1_DCACHE_MISS");
      strcpy(label[51], "BGP_PU1_DCACHE_HIT");
      strcpy(label[52], "BGP_PU1_DATA_LOADS");
      strcpy(label[53], "BGP_PU1_DATA_STORES");
      strcpy(label[54], "BGP_PU1_DCACHE_OPS");
      strcpy(label[55], "BGP_PU1_ICACHE_MISS");
      strcpy(label[56], "BGP_PU1_ICACHE_HIT");
      strcpy(label[57], "BGP_PU1_FPU_ADD_SUB_1");
      strcpy(label[58], "BGP_PU1_FPU_MULT_1");
      strcpy(label[59], "BGP_PU1_FPU_FMA_2");
      strcpy(label[60], "BGP_PU1_FPU_DIV_1");
      strcpy(label[61], "BGP_PU1_FPU_OTHER_NON_STORAGE_OPS");
      strcpy(label[62], "BGP_PU1_FPU_ADD_SUB_2");
      strcpy(label[63], "BGP_PU1_FPU_MULT_2");
      strcpy(label[64], "BGP_PU1_FPU_FMA_4");
      strcpy(label[65], "BGP_PU1_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS");
      strcpy(label[66], "BGP_PU1_FPU_QUADWORD_LOADS");
      strcpy(label[67], "BGP_PU1_FPU_OTHER_LOADS");
      strcpy(label[68], "BGP_PU1_FPU_QUADWORD_STORES");
      strcpy(label[69], "BGP_PU1_FPU_OTHER_STORES");
      strcpy(label[70], "BGP_PU0_L1_INVALIDATION_REQUESTS");
      strcpy(label[71], "BGP_PU1_L1_INVALIDATION_REQUESTS");
      strcpy(label[72], "BGP_PU0_L2_VALID_PREFETCH_REQUESTS");
      strcpy(label[73], "BGP_PU0_L2_PREFETCH_HITS_IN_FILTER");
      strcpy(label[74], "BGP_PU0_L2_PREFETCH_HITS_IN_STREAM");
      strcpy(label[75], "BGP_PU0_L2_CYCLES_PREFETCH_PENDING");
      strcpy(label[76], "BGP_PU0_L2_PAGE_ALREADY_IN_L2");
      strcpy(label[77], "BGP_PU0_L2_PREFETCH_SNOOP_HIT_SAME_CORE");
      strcpy(label[78], "BGP_PU0_L2_PREFETCH_SNOOP_HIT_OTHER_CORE");
      strcpy(label[79], "BGP_PU0_L2_PREFETCH_SNOOP_HIT_PLB");
      strcpy(label[80], "BGP_PU0_L2_CYCLES_READ_REQUEST_PENDING");
      strcpy(label[81], "BGP_PU0_L2_READ_REQUESTS");
      strcpy(label[82], "BGP_PU0_L2_DEVBUS_READ_REQUESTS");
      strcpy(label[83], "BGP_PU0_L2_L3_READ_REQUESTS");
      strcpy(label[84], "BGP_PU0_L2_NETBUS_READ_REQUESTS");
      strcpy(label[85], "BGP_PU0_L2_BLIND_DEV_READ_REQUESTS");
      strcpy(label[86], "BGP_PU0_L2_PREFETCHABLE_REQUESTS");
      strcpy(label[87], "BGP_PU0_L2_HIT");
      strcpy(label[88], "BGP_PU0_L2_SAME_CORE_SNOOPS");
      strcpy(label[89], "BGP_PU0_L2_OTHER_CORE_SNOOPS");
      strcpy(label[90], "BGP_PU0_L2_OTHER_DP_PU0_SNOOPS");
      strcpy(label[91], "BGP_PU0_L2_OTHER_DP_PU1_SNOOPS");
      strcpy(label[92], "BGP_RESERVED");
      strcpy(label[93], "BGP_RESERVED");
      strcpy(label[94], "BGP_RESERVED");
      strcpy(label[95], "BGP_RESERVED");
      strcpy(label[96], "BGP_RESERVED");
      strcpy(label[97], "BGP_RESERVED");
      strcpy(label[98], "BGP_RESERVED");
      strcpy(label[99], "BGP_RESERVED");
      strcpy(label[100], "BGP_RESERVED");
      strcpy(label[101], "BGP_PU0_L2_MEMORY_WRITES");
      strcpy(label[102], "BGP_PU0_L2_NETWORK_WRITES");
      strcpy(label[103], "BGP_PU0_L2_DEVBUS_WRITES");
      strcpy(label[104], "BGP_PU1_L2_VALID_PREFETCH_REQUESTS");
      strcpy(label[105], "BGP_PU1_L2_PREFETCH_HITS_IN_FILTER");
      strcpy(label[106], "BGP_PU1_L2_PREFETCH_HITS_IN_STREAM");
      strcpy(label[107], "BGP_PU1_L2_CYCLES_PREFETCH_PENDING");
      strcpy(label[108], "BGP_PU1_L2_PAGE_ALREADY_IN_L2");
      strcpy(label[109], "BGP_PU1_L2_PREFETCH_SNOOP_HIT_SAME_CORE");
      strcpy(label[110], "BGP_PU1_L2_PREFETCH_SNOOP_HIT_OTHER_CORE");
      strcpy(label[111], "BGP_PU1_L2_PREFETCH_SNOOP_HIT_PLB");
      strcpy(label[112], "BGP_PU1_L2_CYCLES_READ_REQUEST_PENDING");
      strcpy(label[113], "BGP_PU1_L2_READ_REQUESTS");
      strcpy(label[114], "BGP_PU1_L2_DEVBUS_READ_REQUESTS");
      strcpy(label[115], "BGP_PU1_L2_L3_READ_REQUESTS");
      strcpy(label[116], "BGP_PU1_L2_NETBUS_READ_REQUESTS");
      strcpy(label[117], "BGP_PU1_L2_BLIND_DEV_READ_REQUESTS");
      strcpy(label[118], "BGP_PU1_L2_PREFETCHABLE_REQUESTS");
      strcpy(label[119], "BGP_PU1_L2_HIT");
      strcpy(label[120], "BGP_PU1_L2_SAME_CORE_SNOOPS");
      strcpy(label[121], "BGP_PU1_L2_OTHER_CORE_SNOOPS");
      strcpy(label[122], "BGP_PU1_L2_OTHER_DP_PU0_SNOOPS");
      strcpy(label[123], "BGP_PU1_L2_OTHER_DP_PU1_SNOOPS");
      strcpy(label[124], "BGP_RESERVED");
      strcpy(label[125], "BGP_RESERVED");
      strcpy(label[126], "BGP_RESERVED");
      strcpy(label[127], "BGP_RESERVED");
      strcpy(label[128], "BGP_RESERVED");
      strcpy(label[129], "BGP_RESERVED");
      strcpy(label[130], "BGP_RESERVED");
      strcpy(label[131], "BGP_RESERVED");
      strcpy(label[132], "BGP_RESERVED");
      strcpy(label[133], "BGP_PU1_L2_MEMORY_WRITES");
      strcpy(label[134], "BGP_PU1_L2_NETWORK_WRITES");
      strcpy(label[135], "BGP_PU1_L2_DEVBUS_WRITES");
      strcpy(label[136], "BGP_L3_M0_RD0_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[137], "BGP_L3_M0_RD0_BURST_DELIVERED_L2");
      strcpy(label[138], "BGP_L3_M0_RD0_READ_RETURN_COLLISION");
      strcpy(label[139], "BGP_L3_M0_RD0_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[140], "BGP_L3_M0_RD0_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[141], "BGP_L3_M0_RD0_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[142], "BGP_L3_M0_RD0_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[143], "BGP_L3_M0_RD1_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[144], "BGP_L3_M0_RD1_BURST_DELIVERED_L2");
      strcpy(label[145], "BGP_L3_M0_RD1_READ_RETURN_COLLISION");
      strcpy(label[146], "BGP_L3_M0_RD1_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[147], "BGP_L3_M0_RD1_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[148], "BGP_L3_M0_RD1_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[149], "BGP_L3_M0_RD1_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[150], "BGP_L3_M0_DIR0_LOOKUPS");
      strcpy(label[151], "BGP_L3_M0_DIR0_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[152], "BGP_L3_M0_DIR1_LOOKUPS");
      strcpy(label[153], "BGP_L3_M0_DIR1_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[154], "BGP_L3_M0_MH_DDR_STORES");
      strcpy(label[155], "BGP_L3_M0_MH_DDR_FETCHES");
      strcpy(label[156], "BGP_L3_M1_RD0_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[157], "BGP_L3_M1_RD0_BURST_DELIVERED_L2");
      strcpy(label[158], "BGP_L3_M1_RD0_READ_RETURN_COLLISION");
      strcpy(label[159], "BGP_L3_M1_RD0_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[160], "BGP_L3_M1_RD0_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[161], "BGP_L3_M1_RD0_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[162], "BGP_L3_M1_RD0_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[163], "BGP_L3_M1_RD1_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[164], "BGP_L3_M1_RD1_BURST_DELIVERED_L2");
      strcpy(label[165], "BGP_L3_M1_RD1_READ_RETURN_COLLISION");
      strcpy(label[166], "BGP_L3_M1_RD1_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[167], "BGP_L3_M1_RD1_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[168], "BGP_L3_M1_RD1_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[169], "BGP_L3_M1_RD1_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[170], "BGP_L3_M1_DIR0_LOOKUPS");
      strcpy(label[171], "BGP_L3_M1_DIR0_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[172], "BGP_L3_M1_DIR1_LOOKUPS");
      strcpy(label[173], "BGP_L3_M1_DIR1_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[174], "BGP_L3_M1_MH_DDR_STORES");
      strcpy(label[175], "BGP_L3_M1_MH_DDR_FETCHES");
      strcpy(label[176], "BGP_PU0_SNOOP_PORT0_REMOTE_SOURCE_REQUESTS");
      strcpy(label[177], "BGP_PU0_SNOOP_PORT1_REMOTE_SOURCE_REQUESTS");
      strcpy(label[178], "BGP_PU0_SNOOP_PORT2_REMOTE_SOURCE_REQUESTS");
      strcpy(label[179], "BGP_PU0_SNOOP_PORT3_REMOTE_SOURCE_REQUESTS");
      strcpy(label[180], "BGP_PU0_SNOOP_PORT0_REJECTED_REQUESTS");
      strcpy(label[181], "BGP_PU0_SNOOP_PORT1_REJECTED_REQUESTS");
      strcpy(label[182], "BGP_PU0_SNOOP_PORT2_REJECTED_REQUESTS");
      strcpy(label[183], "BGP_PU0_SNOOP_PORT3_REJECTED_REQUESTS");
      strcpy(label[184], "BGP_PU0_SNOOP_L1_CACHE_WRAP");
      strcpy(label[185], "BGP_PU1_SNOOP_PORT0_REMOTE_SOURCE_REQUESTS");
      strcpy(label[186], "BGP_PU1_SNOOP_PORT1_REMOTE_SOURCE_REQUESTS");
      strcpy(label[187], "BGP_PU1_SNOOP_PORT2_REMOTE_SOURCE_REQUESTS");
      strcpy(label[188], "BGP_PU1_SNOOP_PORT3_REMOTE_SOURCE_REQUESTS");
      strcpy(label[189], "BGP_PU1_SNOOP_PORT0_REJECTED_REQUESTS");
      strcpy(label[190], "BGP_PU1_SNOOP_PORT1_REJECTED_REQUESTS");
      strcpy(label[191], "BGP_PU1_SNOOP_PORT2_REJECTED_REQUESTS");
      strcpy(label[192], "BGP_PU1_SNOOP_PORT3_REJECTED_REQUESTS");
      strcpy(label[193], "BGP_PU1_SNOOP_L1_CACHE_WRAP");
      strcpy(label[194], "BGP_TORUS_XP_PACKETS");
      strcpy(label[195], "BGP_TORUS_XP_32BCHUNKS");
      strcpy(label[196], "BGP_TORUS_XM_PACKETS");
      strcpy(label[197], "BGP_TORUS_XM_32BCHUNKS");
      strcpy(label[198], "BGP_TORUS_YP_PACKETS");
      strcpy(label[199], "BGP_TORUS_YP_32BCHUNKS");
      strcpy(label[200], "BGP_TORUS_YM_PACKETS");
      strcpy(label[201], "BGP_TORUS_YM_32BCHUNKS");
      strcpy(label[202], "BGP_TORUS_ZP_PACKETS");
      strcpy(label[203], "BGP_TORUS_ZP_32BCHUNKS");
      strcpy(label[204], "BGP_TORUS_ZM_PACKETS");
      strcpy(label[205], "BGP_TORUS_ZM_32BCHUNKS");
      strcpy(label[206], "BGP_DMA_PACKETS_INJECTED");
      strcpy(label[207], "BGP_DMA_DESCRIPTORS_READ_FROM_L3");
      strcpy(label[208], "BGP_DMA_FIFO_PACKETS_RECEIVED");
      strcpy(label[209], "BGP_DMA_COUNTER_PACKETS_RECEIVED");
      strcpy(label[210], "BGP_DMA_REMOTE_GET_PACKETS_RECEIVED");
      strcpy(label[211], "BGP_DMA_IDPU_READ_REQUESTS_TO_L3");
      strcpy(label[212], "BGP_DMA_READ_VALID_RETURNED");
      strcpy(label[213], "BGP_DMA_ACKED_READ_REQUESTS");
      strcpy(label[214], "BGP_DMA_CYCLES_RDPU_WRITE_ACTIVE");
      strcpy(label[215], "BGP_DMA_WRITE_REQUESTS_TO_L3");
      strcpy(label[216], "BGP_RESERVED");
      strcpy(label[217], "BGP_RESERVED");
      strcpy(label[218], "BGP_RESERVED");
      strcpy(label[219], "BGP_RESERVED");
      strcpy(label[220], "BGP_RESERVED");
      strcpy(label[221], "BGP_RESERVED");
      strcpy(label[222], "BGP_COL_AC_CH2_VC0_MATURE");
      strcpy(label[223], "BGP_COL_AC_CH1_VC0_MATURE");
      strcpy(label[224], "BGP_COL_AC_CH0_VC0_MATURE");
      strcpy(label[225], "BGP_COL_AC_INJECT_VC0_MATURE");
      strcpy(label[226], "BGP_COL_AC_CH2_VC1_MATURE");
      strcpy(label[227], "BGP_COL_AC_CH1_VC1_MATURE");
      strcpy(label[228], "BGP_COL_AC_CH0_VC1_MATURE");
      strcpy(label[229], "BGP_COL_AC_INJECT_VC1_MATURE");
      strcpy(label[230], "BGP_COL_AC_PENDING_REQUESTS");
      strcpy(label[231], "BGP_COL_AC_WAITING_REQUESTS");
      strcpy(label[232], "BGP_COL_AR2_PACKET_TAKEN");
      strcpy(label[233], "BGP_COL_AR1_PACKET_TAKEN");
      strcpy(label[234], "BGP_COL_AR0_PACKET_TAKEN");
      strcpy(label[235], "BGP_COL_ALC_PACKET_TAKEN");
      strcpy(label[236], "BGP_COL_AR0_VC0_DATA_PACKETS_RECEIVED");
      strcpy(label[237], "BGP_COL_AR0_VC1_DATA_PACKETS_RECEIVED");
      strcpy(label[238], "BGP_COL_AR1_VC0_DATA_PACKETS_RECEIVED");
      strcpy(label[239], "BGP_COL_AR1_VC1_DATA_PACKETS_RECEIVED");
      strcpy(label[240], "BGP_COL_AR2_VC0_DATA_PACKETS_RECEIVED");
      strcpy(label[241], "BGP_COL_AR2_VC1_DATA_PACKETS_RECEIVED");
      strcpy(label[242], "BGP_COL_AS0_VC0_DATA_PACKETS_SENT");
      strcpy(label[243], "BGP_COL_AS0_VC1_DATA_PACKETS_SENT");
      strcpy(label[244], "BGP_COL_AS1_VC0_DATA_PACKETS_SENT");
      strcpy(label[245], "BGP_COL_AS1_VC1_DATA_PACKETS_SENT");
      strcpy(label[246], "BGP_COL_AS2_VC0_DATA_PACKETS_SENT");
      strcpy(label[247], "BGP_COL_AS2_VC1_DATA_PACKETS_SENT");
      strcpy(label[248], "BGP_COL_INJECT_VC0_HEADER");
      strcpy(label[249], "BGP_COL_INJECT_VC1_HEADER");
      strcpy(label[250], "BGP_COL_RECEPTION_VC0_PACKET_ADDED");
      strcpy(label[251], "BGP_COL_RECEPTION_VC1_PACKET_ADDED");
      strcpy(label[252], "BGP_IC_TIMESTAMP");
      strcpy(label[253], "BGP_RESERVED");
      strcpy(label[254], "BGP_RESERVED");
      strcpy(label[255], "BGP_MISC_ELAPSED_TIME");
   }

   else if (mode == 1)
   {
      strcpy(label[0], "BGP_PU2_JPIPE_INSTRUCTIONS");
      strcpy(label[1], "BGP_PU2_JPIPE_ADD_SUB");
      strcpy(label[2], "BGP_PU2_JPIPE_LOGICAL_OPS");
      strcpy(label[3], "BGP_PU2_JPIPE_SHROTMK");
      strcpy(label[4], "BGP_PU2_IPIPE_INSTRUCTIONS");
      strcpy(label[5], "BGP_PU2_IPIPE_MULT_DIV");
      strcpy(label[6], "BGP_PU2_IPIPE_ADD_SUB");
      strcpy(label[7], "BGP_PU2_IPIPE_LOGICAL_OPS");
      strcpy(label[8], "BGP_PU2_IPIPE_SHROTMK");
      strcpy(label[9], "BGP_PU2_IPIPE_BRANCHES");
      strcpy(label[10], "BGP_PU2_IPIPE_TLB_OPS");
      strcpy(label[11], "BGP_PU2_IPIPE_PROCESS_CONTROL");
      strcpy(label[12], "BGP_PU2_IPIPE_OTHER");
      strcpy(label[13], "BGP_PU2_DCACHE_LINEFILLINPROG");
      strcpy(label[14], "BGP_PU2_ICACHE_LINEFILLINPROG");
      strcpy(label[15], "BGP_PU2_DCACHE_MISS");
      strcpy(label[16], "BGP_PU2_DCACHE_HIT");
      strcpy(label[17], "BGP_PU2_DATA_LOADS");
      strcpy(label[18], "BGP_PU2_DATA_STORES");
      strcpy(label[19], "BGP_PU2_DCACHE_OPS");
      strcpy(label[20], "BGP_PU2_ICACHE_MISS");
      strcpy(label[21], "BGP_PU2_ICACHE_HIT");
      strcpy(label[22], "BGP_PU2_FPU_ADD_SUB_1");
      strcpy(label[23], "BGP_PU2_FPU_MULT_1");
      strcpy(label[24], "BGP_PU2_FPU_FMA_2");
      strcpy(label[25], "BGP_PU2_FPU_DIV_1");
      strcpy(label[26], "BGP_PU2_FPU_OTHER_NON_STORAGE_OPS");
      strcpy(label[27], "BGP_PU2_FPU_ADD_SUB_2");
      strcpy(label[28], "BGP_PU2_FPU_MULT_2");
      strcpy(label[29], "BGP_PU2_FPU_FMA_4");
      strcpy(label[30], "BGP_PU2_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS");
      strcpy(label[31], "BGP_PU2_FPU_QUADWORD_LOADS");
      strcpy(label[32], "BGP_PU2_FPU_OTHER_LOADS");
      strcpy(label[33], "BGP_PU2_FPU_QUADWORD_STORES");
      strcpy(label[34], "BGP_PU2_FPU_OTHER_STORES");
      strcpy(label[35], "BGP_PU3_JPIPE_INSTRUCTIONS");
      strcpy(label[36], "BGP_PU3_JPIPE_ADD_SUB");
      strcpy(label[37], "BGP_PU3_JPIPE_LOGICAL_OPS");
      strcpy(label[38], "BGP_PU3_JPIPE_SHROTMK");
      strcpy(label[39], "BGP_PU3_IPIPE_INSTRUCTIONS");
      strcpy(label[40], "BGP_PU3_IPIPE_MULT_DIV");
      strcpy(label[41], "BGP_PU3_IPIPE_ADD_SUB");
      strcpy(label[42], "BGP_PU3_IPIPE_LOGICAL_OPS");
      strcpy(label[43], "BGP_PU3_IPIPE_SHROTMK");
      strcpy(label[44], "BGP_PU3_IPIPE_BRANCHES");
      strcpy(label[45], "BGP_PU3_IPIPE_TLB_OPS");
      strcpy(label[46], "BGP_PU3_IPIPE_PROCESS_CONTROL");
      strcpy(label[47], "BGP_PU3_IPIPE_OTHER");
      strcpy(label[48], "BGP_PU3_DCACHE_LINEFILLINPROG");
      strcpy(label[49], "BGP_PU3_ICACHE_LINEFILLINPROG");
      strcpy(label[50], "BGP_PU3_DCACHE_MISS");
      strcpy(label[51], "BGP_PU3_DCACHE_HIT");
      strcpy(label[52], "BGP_PU3_DATA_LOADS");
      strcpy(label[53], "BGP_PU3_DATA_STORES");
      strcpy(label[54], "BGP_PU3_DCACHE_OPS");
      strcpy(label[55], "BGP_PU3_ICACHE_MISS");
      strcpy(label[56], "BGP_PU3_ICACHE_HIT");
      strcpy(label[57], "BGP_PU3_FPU_ADD_SUB_1");
      strcpy(label[58], "BGP_PU3_FPU_MULT_1");
      strcpy(label[59], "BGP_PU3_FPU_FMA_2");
      strcpy(label[60], "BGP_PU3_FPU_DIV_1");
      strcpy(label[61], "BGP_PU3_FPU_OTHER_NON_STORAGE_OPS");
      strcpy(label[62], "BGP_PU3_FPU_ADD_SUB_2");
      strcpy(label[63], "BGP_PU3_FPU_MULT_2");
      strcpy(label[64], "BGP_PU3_FPU_FMA_4");
      strcpy(label[65], "BGP_PU3_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS");
      strcpy(label[66], "BGP_PU3_FPU_QUADWORD_LOADS");
      strcpy(label[67], "BGP_PU3_FPU_OTHER_LOADS");
      strcpy(label[68], "BGP_PU3_FPU_QUADWORD_STORES");
      strcpy(label[69], "BGP_PU3_FPU_OTHER_STORES");
      strcpy(label[70], "BGP_PU2_L1_INVALIDATION_REQUESTS");
      strcpy(label[71], "BGP_PU3_L1_INVALIDATION_REQUESTS");
      strcpy(label[72], "BGP_COL_AC_CH2_VC0_MATURE_UM1");
      strcpy(label[73], "BGP_COL_AC_CH1_VC0_MATURE_UM1");
      strcpy(label[74], "BGP_COL_AC_CH0_VC0_MATURE_UM1");
      strcpy(label[75], "BGP_COL_AC_INJECT_VC0_MATURE_UM1");
      strcpy(label[76], "BGP_COL_AC_CH2_VC1_MATURE_UM1");
      strcpy(label[77], "BGP_COL_AC_CH1_VC1_MATURE_UM1");
      strcpy(label[78], "BGP_COL_AC_CH0_VC1_MATURE_UM1");
      strcpy(label[79], "BGP_COL_AC_INJECT_VC1_MATURE_UM1");
      strcpy(label[80], "BGP_COL_AR0_VC0_EMPTY_PACKET");
      strcpy(label[81], "BGP_COL_AR0_VC1_EMPTY_PACKET");
      strcpy(label[82], "BGP_COL_AR0_IDLE_PACKET");
      strcpy(label[83], "BGP_COL_AR0_BAD_PACKET_MARKER");
      strcpy(label[84], "BGP_COL_AR0_VC0_CUT_THROUGH");
      strcpy(label[85], "BGP_COL_AR0_VC1_CUT_THROUGH");
      strcpy(label[86], "BGP_COL_AR0_HEADER_PARITY_ERROR");
      strcpy(label[87], "BGP_COL_AR0_UNEXPECTED_HEADER_ERROR");
      strcpy(label[88], "BGP_COL_AR0_RESYNC");
      strcpy(label[89], "BGP_COL_AR1_VC0_EMPTY_PACKET");
      strcpy(label[90], "BGP_COL_AR1_VC1_EMPTY_PACKET");
      strcpy(label[91], "BGP_COL_AR1_IDLE_PACKET");
      strcpy(label[92], "BGP_COL_AR1_BAD_PACKET_MARKER");
      strcpy(label[93], "BGP_COL_AR1_VC0_CUT_THROUGH");
      strcpy(label[94], "BGP_COL_AR1_VC1_CUT_THROUGH");
      strcpy(label[95], "BGP_COL_AR1_HEADER_PARITY_ERROR");
      strcpy(label[96], "BGP_COL_AR1_UNEXPECTED_HEADER_ERROR");
      strcpy(label[97], "BGP_COL_AR1_RESYNC");
      strcpy(label[98], "BGP_COL_AR2_VC0_EMPTY_PACKET");
      strcpy(label[99], "BGP_COL_AR2_VC1_EMPTY_PACKET");
      strcpy(label[100], "BGP_COL_AR2_IDLE_PACKET");
      strcpy(label[101], "BGP_COL_AR2_BAD_PACKET_MARKER");
      strcpy(label[102], "BGP_COL_AR2_VC0_CUT_THROUGH");
      strcpy(label[103], "BGP_COL_AR2_VC1_CUT_THROUGH");
      strcpy(label[104], "BGP_COL_AR2_HEADER_PARITY_ERROR");
      strcpy(label[105], "BGP_COL_AR2_UNEXPECTED_HEADER_ERROR");
      strcpy(label[106], "BGP_COL_AR2_RESYNC");
      strcpy(label[107], "BGP_COL_AS0_VC0_CUT_THROUGH");
      strcpy(label[108], "BGP_COL_AS0_VC1_CUT_THROUGH");
      strcpy(label[109], "BGP_COL_AS0_VC0_PACKETS_SENT");
      strcpy(label[110], "BGP_COL_AS0_VC1_PACKETS_SENT");
      strcpy(label[111], "BGP_COL_AS0_IDLE_PACKETS_SENT");
      strcpy(label[112], "BGP_COL_AS1_VC0_CUT_THROUGH");
      strcpy(label[113], "BGP_COL_AS1_VC1_CUT_THROUGH");
      strcpy(label[114], "BGP_COL_AS1_VC0_PACKETS_SENT");
      strcpy(label[115], "BGP_COL_AS1_VC1_PACKETS_SENT");
      strcpy(label[116], "BGP_COL_AS1_IDLE_PACKETS_SENT");
      strcpy(label[117], "BGP_COL_AS2_VC0_CUT_THROUGH");
      strcpy(label[118], "BGP_COL_AS2_VC1_CUT_THROUGH");
      strcpy(label[119], "BGP_COL_AS2_VC0_PACKETS_SENT");
      strcpy(label[120], "BGP_COL_AS2_VC1_PACKETS_SENT");
      strcpy(label[121], "BGP_COL_AS2_IDLE_PACKETS_SENT");
      strcpy(label[122], "BGP_COL_INJECT_VC0_PAYLOAD_ADDED");
      strcpy(label[123], "BGP_COL_INJECT_VC1_PAYLOAD_ADDED");
      strcpy(label[124], "BGP_COL_INJECT_VC0_PACKET_TAKEN");
      strcpy(label[125], "BGP_COL_INJECT_VC1_PACKET_TAKEN");
      strcpy(label[126], "BGP_COL_RECEPTION_VC0_HEADER_TAKEN");
      strcpy(label[127], "BGP_COL_RECEPTION_VC1_HEADER_TAKEN");
      strcpy(label[128], "BGP_COL_RECEPTION_VC0_PAYLOAD_TAKEN");
      strcpy(label[129], "BGP_COL_RECEPTION_VC1_PAYLOAD_TAKEN");
      strcpy(label[130], "BGP_COL_RECEPTION_VC0_PACKET_DISCARDED");
      strcpy(label[131], "BGP_COL_RECEPTION_VC1_PACKET_DISCARDED");
      strcpy(label[132], "BGP_PU2_L2_VALID_PREFETCH_REQUESTS");
      strcpy(label[133], "BGP_PU2_L2_PREFETCH_HITS_IN_FILTER");
      strcpy(label[134], "BGP_PU2_L2_PREFETCH_HITS_IN_STREAM");
      strcpy(label[135], "BGP_PU2_L2_CYCLES_PREFETCH_PENDING");
      strcpy(label[136], "BGP_PU2_L2_PAGE_ALREADY_IN_L2");
      strcpy(label[137], "BGP_PU2_L2_PREFETCH_SNOOP_HIT_SAME_CORE");
      strcpy(label[138], "BGP_PU2_L2_PREFETCH_SNOOP_HIT_OTHER_CORE");
      strcpy(label[139], "BGP_PU2_L2_PREFETCH_SNOOP_HIT_PLB");
      strcpy(label[140], "BGP_PU2_L2_CYCLES_READ_REQUEST_PENDING");
      strcpy(label[141], "BGP_PU2_L2_READ_REQUESTS");
      strcpy(label[142], "BGP_PU2_L2_DEVBUS_READ_REQUESTS");
      strcpy(label[143], "BGP_PU2_L2_L3_READ_REQUESTS");
      strcpy(label[144], "BGP_PU2_L2_NETBUS_READ_REQUESTS");
      strcpy(label[145], "BGP_PU2_L2_BLIND_DEV_READ_REQUESTS");
      strcpy(label[146], "BGP_PU2_L2_PREFETCHABLE_REQUESTS");
      strcpy(label[147], "BGP_PU2_L2_HIT");
      strcpy(label[148], "BGP_PU2_L2_SAME_CORE_SNOOPS");
      strcpy(label[149], "BGP_PU2_L2_OTHER_CORE_SNOOPS");
      strcpy(label[150], "BGP_PU2_L2_OTHER_DP_PU0_SNOOPS");
      strcpy(label[151], "BGP_PU2_L2_OTHER_DP_PU1_SNOOPS");
      strcpy(label[152], "BGP_RESERVED");
      strcpy(label[153], "BGP_RESERVED");
      strcpy(label[154], "BGP_RESERVED");
      strcpy(label[155], "BGP_RESERVED");
      strcpy(label[156], "BGP_RESERVED");
      strcpy(label[157], "BGP_RESERVED");
      strcpy(label[158], "BGP_RESERVED");
      strcpy(label[159], "BGP_RESERVED");
      strcpy(label[160], "BGP_RESERVED");
      strcpy(label[161], "BGP_PU2_L2_MEMORY_WRITES");
      strcpy(label[162], "BGP_PU2_L2_NETWORK_WRITES");
      strcpy(label[163], "BGP_PU2_L2_DEVBUS_WRITES");
      strcpy(label[164], "BGP_PU3_L2_VALID_PREFETCH_REQUESTS");
      strcpy(label[165], "BGP_PU3_L2_PREFETCH_HITS_IN_FILTER");
      strcpy(label[166], "BGP_PU3_L2_PREFETCH_HITS_IN_STREAM");
      strcpy(label[167], "BGP_PU3_L2_CYCLES_PREFETCH_PENDING");
      strcpy(label[168], "BGP_PU3_L2_PAGE_ALREADY_IN_L2");
      strcpy(label[169], "BGP_PU3_L2_PREFETCH_SNOOP_HIT_SAME_CORE");
      strcpy(label[170], "BGP_PU3_L2_PREFETCH_SNOOP_HIT_OTHER_CORE");
      strcpy(label[171], "BGP_PU3_L2_PREFETCH_SNOOP_HIT_PLB");
      strcpy(label[172], "BGP_PU3_L2_CYCLES_READ_REQUEST_PENDING");
      strcpy(label[173], "BGP_PU3_L2_READ_REQUESTS");
      strcpy(label[174], "BGP_PU3_L2_DEVBUS_READ_REQUESTS");
      strcpy(label[175], "BGP_PU3_L2_L3_READ_REQUESTS");
      strcpy(label[176], "BGP_PU3_L2_NETBUS_READ_REQUESTS");
      strcpy(label[177], "BGP_PU3_L2_BLIND_DEV_READ_REQUESTS");
      strcpy(label[178], "BGP_PU3_L2_PREFETCHABLE_REQUESTS");
      strcpy(label[179], "BGP_PU3_L2_HIT");
      strcpy(label[180], "BGP_PU3_L2_SAME_CORE_SNOOPS");
      strcpy(label[181], "BGP_PU3_L2_OTHER_CORE_SNOOPS");
      strcpy(label[182], "BGP_PU3_L2_OTHER_DP_PU0_SNOOPS");
      strcpy(label[183], "BGP_PU3_L2_OTHER_DP_PU1_SNOOPS");
      strcpy(label[184], "BGP_Reserved");
      strcpy(label[185], "BGP_Reserved");
      strcpy(label[186], "BGP_Reserved");
      strcpy(label[187], "BGP_Reserved");
      strcpy(label[188], "BGP_Reserved");
      strcpy(label[189], "BGP_Reserved");
      strcpy(label[190], "BGP_Reserved");
      strcpy(label[191], "BGP_Reserved");
      strcpy(label[192], "BGP_Reserved");
      strcpy(label[193], "BGP_PU3_L2_MEMORY_WRITES");
      strcpy(label[194], "BGP_PU3_L2_NETWORK_WRITES");
      strcpy(label[195], "BGP_PU3_L2_DEVBUS_WRITES");
      strcpy(label[196], "BGP_L3_M0_R2_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[197], "BGP_L3_M0_R2_BURST_DELIVERED_L2");
      strcpy(label[198], "BGP_L3_M0_R2_READ_RETURN_COLLISION");
      strcpy(label[199], "BGP_L3_M0_R2_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[200], "BGP_L3_M0_R2_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[201], "BGP_L3_M0_R2_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[202], "BGP_L3_M0_R2_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[203], "BGP_L3_M0_W0_DEPOSIT_REQUESTS");
      strcpy(label[204], "BGP_L3_M0_W0_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[205], "BGP_L3_M0_W1_DEPOSIT_REQUESTS");
      strcpy(label[206], "BGP_L3_M0_W1_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[207], "BGP_L3_M0_MH_ALLOCATION_REQUESTS");
      strcpy(label[208], "BGP_L3_M0_MH_CYCLES_ALLOCATION_REQUESTS_NOT_TAKEN");
      strcpy(label[209], "BGP_L3_M0_PF_PREFETCH_INTO_EDRAM");
      strcpy(label[210], "BGP_RESERVED");
      strcpy(label[211], "BGP_RESERVED");
      strcpy(label[212], "BGP_RESERVED");
      strcpy(label[213], "BGP_RESERVED");
      strcpy(label[214], "BGP_RESERVED");
      strcpy(label[215], "BGP_RESERVED");
      strcpy(label[216], "BGP_L3_M1_R2_SINGLE_LINE_DELIVERED_L2");
      strcpy(label[217], "BGP_L3_M1_R2_BURST_DELIVERED_L2");
      strcpy(label[218], "BGP_L3_M1_R2_READ_RETURN_COLLISION");
      strcpy(label[219], "BGP_L3_M1_R2_DIR0_HIT_OR_INFLIGHT");
      strcpy(label[220], "BGP_L3_M1_R2_DIR0_MISS_OR_LOCKDOWN");
      strcpy(label[221], "BGP_L3_M1_R2_DIR1_HIT_OR_INFLIGHT");
      strcpy(label[222], "BGP_L3_M1_R2_DIR1_MISS_OR_LOCKDOWN");
      strcpy(label[223], "BGP_L3_M1_W0_DEPOSIT_REQUESTS");
      strcpy(label[224], "BGP_L3_M1_W0_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[225], "BGP_L3_M1_W1_DEPOSIT_REQUESTS");
      strcpy(label[226], "BGP_L3_M1_W1_CYCLES_REQUESTS_NOT_TAKEN");
      strcpy(label[227], "BGP_L3_M1_MH_ALLOCATION_REQUESTS");
      strcpy(label[228], "BGP_L3_M1_MH_CYCLES_ALLOCATION_REQUESTS_NOT_TAKEN");
      strcpy(label[229], "BGP_L3_M1_PF_PREFETCH_INTO_EDRAM");
      strcpy(label[230], "BGP_Reserved");
      strcpy(label[231], "BGP_Reserved");
      strcpy(label[232], "BGP_Reserved");
      strcpy(label[233], "BGP_Reserved");
      strcpy(label[234], "BGP_Reserved");
      strcpy(label[235], "BGP_Reserved");
      strcpy(label[236], "BGP_PU2_SNOOP_PORT0_REMOTE_SOURCE_REQUESTS");
      strcpy(label[237], "BGP_PU2_SNOOP_PORT1_REMOTE_SOURCE_REQUESTS");
      strcpy(label[238], "BGP_PU2_SNOOP_PORT2_REMOTE_SOURCE_REQUESTS");
      strcpy(label[239], "BGP_PU2_SNOOP_PORT3_REMOTE_SOURCE_REQUESTS");
      strcpy(label[240], "BGP_PU2_SNOOP_PORT0_REJECTED_REQUESTS");
      strcpy(label[241], "BGP_PU2_SNOOP_PORT1_REJECTED_REQUESTS");
      strcpy(label[242], "BGP_PU2_SNOOP_PORT2_REJECTED_REQUESTS");
      strcpy(label[243], "BGP_PU2_SNOOP_PORT3_REJECTED_REQUESTS");
      strcpy(label[244], "BGP_PU2_SNOOP_L1_CACHE_WRAP");
      strcpy(label[245], "BGP_PU3_SNOOP_PORT0_REMOTE_SOURCE_REQUESTS");
      strcpy(label[246], "BGP_PU3_SNOOP_PORT1_REMOTE_SOURCE_REQUESTS");
      strcpy(label[247], "BGP_PU3_SNOOP_PORT2_REMOTE_SOURCE_REQUESTS");
      strcpy(label[248], "BGP_PU3_SNOOP_PORT3_REMOTE_SOURCE_REQUESTS");
      strcpy(label[249], "BGP_PU3_SNOOP_PORT0_REJECTED_REQUESTS");
      strcpy(label[250], "BGP_PU3_SNOOP_PORT1_REJECTED_REQUESTS");
      strcpy(label[251], "BGP_PU3_SNOOP_PORT2_REJECTED_REQUESTS");
      strcpy(label[252], "BGP_PU3_SNOOP_PORT3_REJECTED_REQUESTS");
      strcpy(label[253], "BGP_PU3_SNOOP_L1_CACHE_WRAP");
      strcpy(label[254], "BGP_Reserved");
      strcpy(label[255], "BGP_MISC_ELAPSED_TIME_UM1");
   }

   else if (mode == 2)
   {
      strcpy(label[0], "BGP_PU0_JPIPE_INSTRUCTIONS_UM2");
      strcpy(label[1], "BGP_PU0_JPIPE_ADD_SUB_UM2");
      strcpy(label[2], "BGP_PU0_JPIPE_LOGICAL_OPS_UM2");
      strcpy(label[3], "BGP_PU0_JPIPE_SHROTMK_UM2");
      strcpy(label[4], "BGP_PU0_IPIPE_INSTRUCTIONS_UM2");
      strcpy(label[5], "BGP_PU0_IPIPE_MULT_DIV_UM2");
      strcpy(label[6], "BGP_PU0_IPIPE_ADD_SUB_UM2");
      strcpy(label[7], "BGP_PU0_IPIPE_LOGICAL_OPS_UM2");
      strcpy(label[8], "BGP_PU0_IPIPE_SHROTMK_UM2");
      strcpy(label[9], "BGP_PU0_IPIPE_BRANCHES_UM2");
      strcpy(label[10], "BGP_PU0_IPIPE_TLB_OPS_UM2");
      strcpy(label[11], "BGP_PU0_IPIPE_PROCESS_CONTROL_UM2");
      strcpy(label[12], "BGP_PU0_IPIPE_OTHER_UM2");
      strcpy(label[13], "BGP_PU0_DCACHE_LINEFILLINPROG_UM2");
      strcpy(label[14], "BGP_PU0_ICACHE_LINEFILLINPROG_UM2");
      strcpy(label[15], "BGP_PU0_DCACHE_MISS_UM2");
      strcpy(label[16], "BGP_PU0_DCACHE_HIT_UM2");
      strcpy(label[17], "BGP_PU0_DATA_LOADS_UM2");
      strcpy(label[18], "BGP_PU0_DATA_STORES_UM2");
      strcpy(label[19], "BGP_PU0_DCACHE_OPS_UM2");
      strcpy(label[20], "BGP_PU0_ICACHE_MISS_UM2");
      strcpy(label[21], "BGP_PU0_ICACHE_HIT_UM2");
      strcpy(label[22], "BGP_PU0_FPU_ADD_SUB_1_UM2");
      strcpy(label[23], "BGP_PU0_FPU_MULT_1_UM2");
      strcpy(label[24], "BGP_PU0_FPU_FMA_2_UM2");
      strcpy(label[25], "BGP_PU0_FPU_DIV_1_UM2");
      strcpy(label[26], "BGP_PU0_FPU_OTHER_NON_STORAGE_OPS_UM2");
      strcpy(label[27], "BGP_PU0_FPU_ADD_SUB_2_UM2");
      strcpy(label[28], "BGP_PU0_FPU_MULT_2_UM2");
      strcpy(label[29], "BGP_PU0_FPU_FMA_4_UM2");
      strcpy(label[30], "BGP_PU0_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS_UM2");
      strcpy(label[31], "BGP_PU0_FPU_QUADWORD_LOADS_UM2");
      strcpy(label[32], "BGP_PU0_FPU_OTHER_LOADS_UM2");
      strcpy(label[33], "BGP_PU0_FPU_QUADWORD_STORES_UM2");
      strcpy(label[34], "BGP_PU0_FPU_OTHER_STORES_UM2");
      strcpy(label[35], "BGP_PU1_JPIPE_INSTRUCTIONS_UM2");
      strcpy(label[36], "BGP_PU1_JPIPE_ADD_SUB_UM2");
      strcpy(label[37], "BGP_PU1_JPIPE_LOGICAL_OPS_UM2");
      strcpy(label[38], "BGP_PU1_JPIPE_SHROTMK_UM2");
      strcpy(label[39], "BGP_PU1_IPIPE_INSTRUCTIONS_UM2");
      strcpy(label[40], "BGP_PU1_IPIPE_MULT_DIV_UM2");
      strcpy(label[41], "BGP_PU1_IPIPE_ADD_SUB_UM2");
      strcpy(label[42], "BGP_PU1_IPIPE_LOGICAL_OPS_UM2");
      strcpy(label[43], "BGP_PU1_IPIPE_SHROTMK_UM2");
      strcpy(label[44], "BGP_PU1_IPIPE_BRANCHES_UM2");
      strcpy(label[45], "BGP_PU1_IPIPE_TLB_OPS_UM2");
      strcpy(label[46], "BGP_PU1_IPIPE_PROCESS_CONTROL_UM2");
      strcpy(label[47], "BGP_PU1_IPIPE_OTHER_UM2");
      strcpy(label[48], "BGP_PU1_DCACHE_LINEFILLINPROG_UM2");
      strcpy(label[49], "BGP_PU1_ICACHE_LINEFILLINPROG_UM2");
      strcpy(label[50], "BGP_PU1_DCACHE_MISS_UM2");
      strcpy(label[51], "BGP_PU1_DCACHE_HIT_UM2");
      strcpy(label[52], "BGP_PU1_DATA_LOADS_UM2");
      strcpy(label[53], "BGP_PU1_DATA_STORES_UM2");
      strcpy(label[54], "BGP_PU1_DCACHE_OPS_UM2");
      strcpy(label[55], "BGP_PU1_ICACHE_MISS_UM2");
      strcpy(label[56], "BGP_PU1_ICACHE_HIT_UM2");
      strcpy(label[57], "BGP_PU1_FPU_ADD_SUB_1_UM2");
      strcpy(label[58], "BGP_PU1_FPU_MULT_1_UM2");
      strcpy(label[59], "BGP_PU1_FPU_FMA_2_UM2");
      strcpy(label[60], "BGP_PU1_FPU_DIV_1_UM2");
      strcpy(label[61], "BGP_PU1_FPU_OTHER_NON_STORAGE_OPS_UM2");
      strcpy(label[62], "BGP_PU1_FPU_ADD_SUB_2_UM2");
      strcpy(label[63], "BGP_PU1_FPU_MULT_2_UM2");
      strcpy(label[64], "BGP_PU1_FPU_FMA_4_UM2");
      strcpy(label[65], "BGP_PU1_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS_UM2");
      strcpy(label[66], "BGP_PU1_FPU_QUADWORD_LOADS_UM2");
      strcpy(label[67], "BGP_PU1_FPU_OTHER_LOADS_UM2");
      strcpy(label[68], "BGP_PU1_FPU_QUADWORD_STORES_UM2");
      strcpy(label[69], "BGP_PU1_FPU_OTHER_STORES_UM2");
      strcpy(label[70], "BGP_PU0_L1_INVALIDATION_UM2");
      strcpy(label[71], "BGP_PU1_L1_INVALIDATION_UM2");
      strcpy(label[72], "BGP_PU0_SNOOP_PORT0_CACHE_REJECTED_REQUEST");
      strcpy(label[73], "BGP_PU0_SNOOP_PORT1_CACHE_REJECTED_REQUEST");
      strcpy(label[74], "BGP_PU0_SNOOP_PORT2_CACHE_REJECTED_REQUEST");
      strcpy(label[75], "BGP_PU0_SNOOP_PORT3_CACHE_REJECTED_REQUEST");
      strcpy(label[76], "BGP_PU0_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[77], "BGP_PU0_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[78], "BGP_PU0_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[79], "BGP_PU0_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[80], "BGP_PU0_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[81], "BGP_PU0_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[82], "BGP_PU0_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[83], "BGP_PU0_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[84], "BGP_PU0_SNOOP_PORT0_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[85], "BGP_PU0_SNOOP_PORT1_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[86], "BGP_PU0_SNOOP_PORT2_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[87], "BGP_PU0_SNOOP_PORT3_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[88], "BGP_PU0_SNOOP_PORT0_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[89], "BGP_PU0_SNOOP_PORT1_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[90], "BGP_PU0_SNOOP_PORT2_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[91], "BGP_PU0_SNOOP_PORT3_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[92], "BGP_PU0_SNOOP_PORT0_UPDATED_CACHE_LINE");
      strcpy(label[93], "BGP_PU0_SNOOP_PORT1_UPDATED_CACHE_LINE");
      strcpy(label[94], "BGP_PU0_SNOOP_PORT2_UPDATED_CACHE_LINE");
      strcpy(label[95], "BGP_PU0_SNOOP_PORT3_UPDATED_CACHE_LINE");
      strcpy(label[96], "BGP_PU0_SNOOP_PORT0_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[97], "BGP_PU0_SNOOP_PORT1_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[98], "BGP_PU0_SNOOP_PORT2_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[99], "BGP_PU0_SNOOP_PORT3_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[100], "BGP_PU1_SNOOP_PORT0_CACHE_REJECTED_REQUEST");
      strcpy(label[101], "BGP_PU1_SNOOP_PORT1_CACHE_REJECTED_REQUEST");
      strcpy(label[102], "BGP_PU1_SNOOP_PORT2_CACHE_REJECTED_REQUEST");
      strcpy(label[103], "BGP_PU1_SNOOP_PORT3_CACHE_REJECTED_REQUEST");
      strcpy(label[104], "BGP_PU1_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[105], "BGP_PU1_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[106], "BGP_PU1_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[107], "BGP_PU1_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[108], "BGP_PU1_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[109], "BGP_PU1_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[110], "BGP_PU1_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[111], "BGP_PU1_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[112], "BGP_PU1_SNOOP_PORT0_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[113], "BGP_PU1_SNOOP_PORT1_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[114], "BGP_PU1_SNOOP_PORT2_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[115], "BGP_PU1_SNOOP_PORT3_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[116], "BGP_PU1_SNOOP_PORT0_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[117], "BGP_PU1_SNOOP_PORT1_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[118], "BGP_PU1_SNOOP_PORT2_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[119], "BGP_PU1_SNOOP_PORT3_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[120], "BGP_PU1_SNOOP_PORT0_UPDATED_CACHE_LINE");
      strcpy(label[121], "BGP_PU1_SNOOP_PORT1_UPDATED_CACHE_LINE");
      strcpy(label[122], "BGP_PU1_SNOOP_PORT2_UPDATED_CACHE_LINE");
      strcpy(label[123], "BGP_PU1_SNOOP_PORT3_UPDATED_CACHE_LINE");
      strcpy(label[124], "BGP_PU1_SNOOP_PORT0_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[125], "BGP_PU1_SNOOP_PORT1_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[126], "BGP_PU1_SNOOP_PORT2_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[127], "BGP_PU1_SNOOP_PORT3_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[128], "BGP_TORUS_XP_TOKEN_ACK_PACKETS");
      strcpy(label[129], "BGP_TORUS_XP_ACKS");
      strcpy(label[130], "BGP_TORUS_XP_VCD0_32BCHUNKS");
      strcpy(label[131], "BGP_TORUS_XP_VCD1_32BCHUNKS");
      strcpy(label[132], "BGP_TORUS_XP_VCBN_32BCHUNKS");
      strcpy(label[133], "BGP_TORUS_XP_VCBP_32BCHUNKS");
      strcpy(label[134], "BGP_TORUS_XP_NO_TOKENS");
      strcpy(label[135], "BGP_TORUS_XP_NO_VCD0_TOKENS");
      strcpy(label[136], "BGP_TORUS_XP_NO_VCBN_TOKENS");
      strcpy(label[137], "BGP_TORUS_XP_NO_VCBP_TOKENS");
      strcpy(label[138], "BGP_TORUS_XM_TOKEN_ACK_PACKETS");
      strcpy(label[139], "BGP_TORUS_XM_ACKS");
      strcpy(label[140], "BGP_TORUS_XM_VCD0_32BCHUNKS");
      strcpy(label[141], "BGP_TORUS_XM_VCD1_32BCHUNKS");
      strcpy(label[142], "BGP_TORUS_XM_VCBN_32BCHUNKS");
      strcpy(label[143], "BGP_TORUS_XM_VCBP_32BCHUNKS");
      strcpy(label[144], "BGP_TORUS_XM_NO_TOKENS");
      strcpy(label[145], "BGP_TORUS_XM_NO_VCD0_TOKENS");
      strcpy(label[146], "BGP_TORUS_XM_NO_VCBN_TOKENS");
      strcpy(label[147], "BGP_TORUS_XM_NO_VCBP_TOKENS");
      strcpy(label[148], "BGP_TORUS_YP_TOKEN_ACK_PACKETS");
      strcpy(label[149], "BGP_TORUS_YP_ACKS");
      strcpy(label[150], "BGP_TORUS_YP_VCD0_32BCHUNKS");
      strcpy(label[151], "BGP_TORUS_YP_VCD1_32BCHUNKS");
      strcpy(label[152], "BGP_TORUS_YP_VCBN_32BCHUNKS");
      strcpy(label[153], "BGP_TORUS_YP_VCBP_32BCHUNKS");
      strcpy(label[154], "BGP_TORUS_YP_NO_TOKENS");
      strcpy(label[155], "BGP_TORUS_YP_NO_VCD0_TOKENS");
      strcpy(label[156], "BGP_TORUS_YP_NO_VCBN_TOKENS");
      strcpy(label[157], "BGP_TORUS_YP_NO_VCBP_TOKENS");
      strcpy(label[158], "BGP_TORUS_YM_TOKEN_ACK_PACKETS");
      strcpy(label[159], "BGP_TORUS_YM_ACKS");
      strcpy(label[160], "BGP_TORUS_YM_VCD0_32BCHUNKS");
      strcpy(label[161], "BGP_TORUS_YM_VCD1_32BCHUNKS");
      strcpy(label[162], "BGP_TORUS_YM_VCBN_32BCHUNKS");
      strcpy(label[163], "BGP_TORUS_YM_VCBP_32BCHUNKS");
      strcpy(label[164], "BGP_TORUS_YM_NO_TOKENS");
      strcpy(label[165], "BGP_TORUS_YM_NO_VCD0_TOKENS");
      strcpy(label[166], "BGP_TORUS_YM_NO_VCBN_TOKENS");
      strcpy(label[167], "BGP_TORUS_YM_NO_VCBP_TOKENS");
      strcpy(label[168], "BGP_TORUS_ZP_TOKEN_ACK_PACKETS");
      strcpy(label[169], "BGP_TORUS_ZP_ACKS");
      strcpy(label[170], "BGP_TORUS_ZP_VCD0_32BCHUNKS");
      strcpy(label[171], "BGP_TORUS_ZP_VCD1_32BCHUNKS");
      strcpy(label[172], "BGP_TORUS_ZP_VCBN_32BCHUNKS");
      strcpy(label[173], "BGP_TORUS_ZP_VCBP_32BCHUNKS");
      strcpy(label[174], "BGP_TORUS_ZP_NO_TOKENS");
      strcpy(label[175], "BGP_TORUS_ZP_NO_VCD0_TOKENS");
      strcpy(label[176], "BGP_TORUS_ZP_NO_VCBN_TOKENS");
      strcpy(label[177], "BGP_TORUS_ZP_NO_VCBP_TOKENS");
      strcpy(label[178], "BGP_TORUS_ZM_TOKEN_ACK_PACKETS");
      strcpy(label[179], "BGP_TORUS_ZM_ACKS");
      strcpy(label[180], "BGP_TORUS_ZM_VCD0_32BCHUNKS");
      strcpy(label[181], "BGP_TORUS_ZM_VCD1_32BCHUNKS");
      strcpy(label[182], "BGP_TORUS_ZM_VCBN_32BCHUNKS");
      strcpy(label[183], "BGP_TORUS_ZM_VCBP_32BCHUNKS");
      strcpy(label[184], "BGP_TORUS_ZM_NO_TOKENS");
      strcpy(label[185], "BGP_TORUS_ZM_NO_VCD0_TOKENS");
      strcpy(label[186], "BGP_RESERVED");
      strcpy(label[187], "BGP_RESERVED");
      strcpy(label[188], "BGP_RESERVED");
      strcpy(label[189], "BGP_RESERVED");
      strcpy(label[190], "BGP_RESERVED");
      strcpy(label[191], "BGP_RESERVED");
      strcpy(label[192], "BGP_RESERVED");
      strcpy(label[193], "BGP_RESERVED");
      strcpy(label[194], "BGP_RESERVED");
      strcpy(label[195], "BGP_RESERVED");
      strcpy(label[196], "BGP_RESERVED");
      strcpy(label[197], "BGP_RESERVED");
      strcpy(label[198], "BGP_RESERVED");
      strcpy(label[199], "BGP_RESERVED");
      strcpy(label[200], "BGP_RESERVED");
      strcpy(label[201], "BGP_RESERVED");
      strcpy(label[202], "BGP_RESERVED");
      strcpy(label[203], "BGP_RESERVED");
      strcpy(label[204], "BGP_COL_AR2_ABORT_UM2");
      strcpy(label[205], "BGP_COL_AR1_ABORT_UM2");
      strcpy(label[206], "BGP_COL_AR0_ABORT_UM2");
      strcpy(label[207], "BGP_COL_A_LOCAL_CLIENT_ABORT");
      strcpy(label[208], "BGP_COL_AR0_VC0_FULL");
      strcpy(label[209], "BGP_COL_AR0_VC1_FULL");
      strcpy(label[210], "BGP_COL_AR1_VC0_FULL");
      strcpy(label[211], "BGP_COL_AR1_VC1_FULL");
      strcpy(label[212], "BGP_COL_AR2_VC0_FULL");
      strcpy(label[213], "BGP_COL_AR2_VC1_FULL");
      strcpy(label[214], "BGP_COL_AS0_VC0_EMPTY");
      strcpy(label[215], "BGP_COL_AS0_VC1_EMPTY");
      strcpy(label[216], "BGP_COL_AS0_RESENDS");
      strcpy(label[217], "BGP_COL_AS1_VC0_EMPTY");
      strcpy(label[218], "BGP_COL_AS1_VC1_EMPTY");
      strcpy(label[219], "BGP_COL_AS1_RESENDS");
      strcpy(label[220], "BGP_COL_AS2_VC0_EMPTY");
      strcpy(label[221], "BGP_COL_AS2_VC1_EMPTY");
      strcpy(label[222], "BGP_COL_AS2_RESENDS");
      strcpy(label[223], "BGP_RESERVED");
      strcpy(label[224], "BGP_RESERVED");
      strcpy(label[225], "BGP_RESERVED");
      strcpy(label[226], "BGP_RESERVED");
      strcpy(label[227], "BGP_RESERVED");
      strcpy(label[228], "BGP_RESERVED");
      strcpy(label[229], "BGP_RESERVED");
      strcpy(label[230], "BGP_RESERVED");
      strcpy(label[231], "BGP_RESERVED");
      strcpy(label[232], "BGP_RESERVED");
      strcpy(label[233], "BGP_RESERVED");
      strcpy(label[234], "BGP_RESERVED");
      strcpy(label[235], "BGP_RESERVED");
      strcpy(label[236], "BGP_RESERVED");
      strcpy(label[237], "BGP_RESERVED");
      strcpy(label[238], "BGP_RESERVED");
      strcpy(label[239], "BGP_RESERVED");
      strcpy(label[240], "BGP_RESERVED");
      strcpy(label[241], "BGP_RESERVED");
      strcpy(label[242], "BGP_RESERVED");
      strcpy(label[243], "BGP_RESERVED");
      strcpy(label[244], "BGP_RESERVED");
      strcpy(label[245], "BGP_RESERVED");
      strcpy(label[246], "BGP_RESERVED");
      strcpy(label[247], "BGP_RESERVED");
      strcpy(label[248], "BGP_RESERVED");
      strcpy(label[249], "BGP_RESERVED");
      strcpy(label[250], "BGP_RESERVED");
      strcpy(label[251], "BGP_RESERVED");
      strcpy(label[252], "BGP_RESERVED");
      strcpy(label[253], "BGP_RESERVED");
      strcpy(label[254], "BGP_RESERVED");
      strcpy(label[255], "BGP_MISC_ELAPSED_TIME_UM2");
   }

   else if (mode == 3)
   {
      strcpy(label[0], "BGP_PU2_JPIPE_INSTRUCTIONS_UM3");
      strcpy(label[1], "BGP_PU2_JPIPE_ADD_SUB_UM3");
      strcpy(label[2], "BGP_PU2_JPIPE_LOGICAL_OPS_UM3");
      strcpy(label[3], "BGP_PU2_JPIPE_SHROTMK_UM3");
      strcpy(label[4], "BGP_PU2_IPIPE_INSTRUCTIONS_UM3");
      strcpy(label[5], "BGP_PU2_IPIPE_MULT_DIV_UM3");
      strcpy(label[6], "BGP_PU2_IPIPE_ADD_SUB_UM3");
      strcpy(label[7], "BGP_PU2_IPIPE_LOGICAL_OPS_UM3");
      strcpy(label[8], "BGP_PU2_IPIPE_SHROTMK_UM3");
      strcpy(label[9], "BGP_PU2_IPIPE_BRANCHES_UM3");
      strcpy(label[10], "BGP_PU2_IPIPE_TLB_OPS_UM3");
      strcpy(label[11], "BGP_PU2_IPIPE_PROCESS_CONTROL_UM3");
      strcpy(label[12], "BGP_PU2_IPIPE_OTHER_UM3");
      strcpy(label[13], "BGP_PU2_DCACHE_LINEFILLINPROG_UM3");
      strcpy(label[14], "BGP_PU2_ICACHE_LINEFILLINPROG_UM3");
      strcpy(label[15], "BGP_PU2_DCACHE_MISS_UM3");
      strcpy(label[16], "BGP_PU2_DCACHE_HIT_UM3");
      strcpy(label[17], "BGP_PU2_DATA_LOADS_UM3");
      strcpy(label[18], "BGP_PU2_DATA_STORES_UM3");
      strcpy(label[19], "BGP_PU2_DCACHE_OPS_UM3");
      strcpy(label[20], "BGP_PU2_ICACHE_MISS_UM3");
      strcpy(label[21], "BGP_PU2_ICACHE_HIT_UM3");
      strcpy(label[22], "BGP_PU2_FPU_ADD_SUB_1_UM3");
      strcpy(label[23], "BGP_PU2_FPU_MULT_1_UM3");
      strcpy(label[24], "BGP_PU2_FPU_FMA_2_UM3");
      strcpy(label[25], "BGP_PU2_FPU_DIV_1_UM3");
      strcpy(label[26], "BGP_PU2_FPU_OTHER_NON_STORAGE_OPS_UM3");
      strcpy(label[27], "BGP_PU2_FPU_ADD_SUB_2_UM3");
      strcpy(label[28], "BGP_PU2_FPU_MULT_2_UM3");
      strcpy(label[29], "BGP_PU2_FPU_FMA_4_UM3");
      strcpy(label[30], "BGP_PU2_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS_UM3");
      strcpy(label[31], "BGP_PU2_FPU_QUADWORD_LOADS_UM3");
      strcpy(label[32], "BGP_PU2_FPU_OTHER_LOADS_UM3");
      strcpy(label[33], "BGP_PU2_FPU_QUADWORD_STORES_UM3");
      strcpy(label[34], "BGP_PU2_FPU_OTHER_STORES_UM3");
      strcpy(label[35], "BGP_PU3_JPIPE_INSTRUCTIONS_UM3");
      strcpy(label[36], "BGP_PU3_JPIPE_ADD_SUB_UM3");
      strcpy(label[37], "BGP_PU3_JPIPE_LOGICAL_OPS_UM3");
      strcpy(label[38], "BGP_PU3_JPIPE_SHROTMK_UM3");
      strcpy(label[39], "BGP_PU3_IPIPE_INSTRUCTIONS_UM3");
      strcpy(label[40], "BGP_PU3_IPIPE_MULT_DIV_UM3");
      strcpy(label[41], "BGP_PU3_IPIPE_ADD_SUB_UM3");
      strcpy(label[42], "BGP_PU3_IPIPE_LOGICAL_OPS_UM3");
      strcpy(label[43], "BGP_PU3_IPIPE_SHROTMK_UM3");
      strcpy(label[44], "BGP_PU3_IPIPE_BRANCHES_UM3");
      strcpy(label[45], "BGP_PU3_IPIPE_TLB_OPS_UM3");
      strcpy(label[46], "BGP_PU3_IPIPE_PROCESS_CONTROL_UM3");
      strcpy(label[47], "BGP_PU3_IPIPE_OTHER_UM3");
      strcpy(label[48], "BGP_PU3_DCACHE_LINEFILLINPROG_UM3");
      strcpy(label[49], "BGP_PU3_ICACHE_LINEFILLINPROG_UM3");
      strcpy(label[50], "BGP_PU3_DCACHE_MISS_UM3");
      strcpy(label[51], "BGP_PU3_DCACHE_HIT_UM3");
      strcpy(label[52], "BGP_PU3_DATA_LOADS_UM3");
      strcpy(label[53], "BGP_PU3_DATA_STORES_UM3");
      strcpy(label[54], "BGP_PU3_DCACHE_OPS_UM3");
      strcpy(label[55], "BGP_PU3_ICACHE_MISS_UM3");
      strcpy(label[56], "BGP_PU3_ICACHE_HIT_UM3");
      strcpy(label[57], "BGP_PU3_FPU_ADD_SUB_1_UM3");
      strcpy(label[58], "BGP_PU3_FPU_MULT_1_UM3");
      strcpy(label[59], "BGP_PU3_FPU_FMA_2_UM3");
      strcpy(label[60], "BGP_PU3_FPU_DIV_1_UM3");
      strcpy(label[61], "BGP_PU3_FPU_OTHER_NON_STORAGE_OPS_UM3");
      strcpy(label[62], "BGP_PU3_FPU_ADD_SUB_2_UM3");
      strcpy(label[63], "BGP_PU3_FPU_MULT_2_UM3");
      strcpy(label[64], "BGP_PU3_FPU_FMA_4_UM3");
      strcpy(label[65], "BGP_PU3_FPU_DUAL_PIPE_OTHER_NON_STORAGE_OPS_UM3");
      strcpy(label[66], "BGP_PU3_FPU_QUADWORD_LOADS_UM3");
      strcpy(label[67], "BGP_PU3_FPU_OTHER_LOADS_UM3");
      strcpy(label[68], "BGP_PU3_FPU_QUADWORD_STORES_UM3");
      strcpy(label[69], "BGP_PU3_FPU_OTHER_STORES_UM3");
      strcpy(label[70], "BGP_PU2_L1_INVALIDATION_UM3");
      strcpy(label[71], "BGP_PU3_L1_INVALIDATION_UM3");
      strcpy(label[72], "BGP_COL_A_CH2_VC0_HAVE");
      strcpy(label[73], "BGP_COL_A_CH1_VC0_HAVE");
      strcpy(label[74], "BGP_COL_A_CH0_VC0_HAVE");
      strcpy(label[75], "BGP_COL_A_INJECT_VC0_HAVE");
      strcpy(label[76], "BGP_COL_A_CH2_VC1_HAVE");
      strcpy(label[77], "BGP_COL_A_CH1_VC1_HAVE");
      strcpy(label[78], "BGP_COL_A_CH0_VC1_HAVE");
      strcpy(label[79], "BGP_COL_A_INJECT_VC1_HAVE");
      strcpy(label[80], "BGP_COL_AC_GREEDY_MODE");
      strcpy(label[81], "BGP_COL_AC_PENDING_REQUESTS_UM3");
      strcpy(label[82], "BGP_COL_AC_WAITING_REQUESTS_UM3");
      strcpy(label[83], "BGP_COL_ACLS0_WINS");
      strcpy(label[84], "BGP_COL_ACLS1_WINS");
      strcpy(label[85], "BGP_COL_ACLS2_WINS");
      strcpy(label[86], "BGP_COL_ACLS3_WINS");
      strcpy(label[87], "BGP_COL_ACLS4_WINS");
      strcpy(label[88], "BGP_COL_ACLS5_WINS");
      strcpy(label[89], "BGP_COL_ACLS6_WINS");
      strcpy(label[90], "BGP_COL_ACLS7_WINS");
      strcpy(label[91], "BGP_COL_ACLS8_WINS");
      strcpy(label[92], "BGP_COL_ACLS9_WINS");
      strcpy(label[93], "BGP_COL_ACLS10_WINS");
      strcpy(label[94], "BGP_COL_ACLS11_WINS");
      strcpy(label[95], "BGP_COL_ACLS12_WINS");
      strcpy(label[96], "BGP_COL_ACLS13_WINS");
      strcpy(label[97], "BGP_COL_ACLS14_WINS");
      strcpy(label[98], "BGP_COL_ACLS15_WINS");
      strcpy(label[99], "BGP_COL_AS2_BUSY");
      strcpy(label[100], "BGP_COL_AS1_BUSY");
      strcpy(label[101], "BGP_COL_AS1_BUSY_RECEPTION");
      strcpy(label[102], "BGP_COL_ALC_BUSY");
      strcpy(label[103], "BGP_COL_AR2_BUSY");
      strcpy(label[104], "BGP_COL_AR1_BUSY");
      strcpy(label[105], "BGP_COL_AR0_BUSY");
      strcpy(label[106], "BGP_COL_ALC_BUSY_INJECT");
      strcpy(label[107], "BGP_COL_ALU_BUSY");
      strcpy(label[108], "BGP_COL_AR2_ABORT_UM3");
      strcpy(label[109], "BGP_COL_AR1_ABORT_UM3");
      strcpy(label[110], "BGP_COL_AR0_ABORT_UM3");
      strcpy(label[111], "BGP_COL_ALC_ABORT");
      strcpy(label[112], "BGP_COL_AR2_PACKET_TAKEN_UM3");
      strcpy(label[113], "BGP_COL_AR1_PACKET_TAKEN_UM3");
      strcpy(label[114], "BGP_COL_AR0_PACKET_TAKEN_UM3");
      strcpy(label[115], "BGP_COL_ALC_PACKET_TAKEN_UM3");
      strcpy(label[116], "BGP_COL_AR0_VC0_DATA_PACKET_RECEIVED");
      strcpy(label[117], "BGP_COL_AR0_VC1_DATA_PACKET_RECEIVED");
      strcpy(label[118], "BGP_COL_AR0_VC1_FULL_UM3");
      strcpy(label[119], "BGP_COL_AR0_HEADER_PARITY_ERROR_UM3");
      strcpy(label[120], "BGP_COL_AR1_VC0_DATA_PACKET_RECEIVED");
      strcpy(label[121], "BGP_COL_AR1_VC1_DATA_PACKET_RECEIVED");
      strcpy(label[122], "BGP_COL_AR1_VC0_FULL_UM3");
      strcpy(label[123], "BGP_COL_AR1_VC1_FULL_UM3");
      strcpy(label[124], "BGP_COL_AR2_VC0_DATA_PACKET_RECEIVED");
      strcpy(label[125], "BGP_COL_AR2_VC1_DATA_PACKET_RECEIVED");
      strcpy(label[126], "BGP_COL_AR2_VC0_FULL_UM3");
      strcpy(label[127], "BGP_COL_AR2_VC1_FULL_UM3");
      strcpy(label[128], "BGP_COL_AS0_VC0_EMPTY_UM3");
      strcpy(label[129], "BGP_COL_AS0_VC1_EMPTY_UM3");
      strcpy(label[130], "BGP_COL_AS0_VC0_DATA_PACKETS_SENT_UM3");
      strcpy(label[131], "BGP_COL_AS0_VC1_DATA_PACKETS_SENT_UM3");
      strcpy(label[132], "BGP_COL_AS0_RESENDS_UM3");
      strcpy(label[133], "BGP_COL_AS1_VC0_EMPTY_UM3");
      strcpy(label[134], "BGP_COL_AS1_VC1_EMPTY_UM3");
      strcpy(label[135], "BGP_COL_AS1_VC0_DATA_PACKETS_SENT_UM3");
      strcpy(label[136], "BGP_COL_AS1_VC1_DATA_PACKETS_SENT_UM3");
      strcpy(label[137], "BGP_COL_AS1_RESENDS_UM3");
      strcpy(label[138], "BGP_COL_AS2_VC0_EMPTY_UM3");
      strcpy(label[139], "BGP_COL_AS2_VC1_EMPTY_UM3");
      strcpy(label[140], "BGP_COL_AS2_VC0_DATA_PACKETS_SENT_UM3");
      strcpy(label[141], "BGP_COL_AS2_VC1_DATA_PACKETS_SENT_UM3");
      strcpy(label[142], "BGP_COL_AS2_RESENDS_UM3");
      strcpy(label[143], "BGP_COL_INJECT_VC0_HEADER_ADDED");
      strcpy(label[144], "BGP_COL_INJECT_VC1_HEADER_ADDED");
      strcpy(label[145], "BGP_COL_RECEPTION_VC0_PACKED_ADDED");
      strcpy(label[146], "BGP_COL_RECEPTION_VC1_PACKED_ADDED");
      strcpy(label[147], "BGP_PU2_SNOOP_PORT0_CACHE_REJECTED_REQUEST");
      strcpy(label[148], "BGP_PU2_SNOOP_PORT1_CACHE_REJECTED_REQUEST");
      strcpy(label[149], "BGP_PU2_SNOOP_PORT2_CACHE_REJECTED_REQUEST");
      strcpy(label[150], "BGP_PU2_SNOOP_PORT3_CACHE_REJECTED_REQUEST");
      strcpy(label[151], "BGP_PU2_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[152], "BGP_PU2_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[153], "BGP_PU2_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[154], "BGP_PU2_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[155], "BGP_PU2_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[156], "BGP_PU2_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[157], "BGP_PU2_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[158], "BGP_PU2_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[159], "BGP_PU2_SNOOP_PORT0_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[160], "BGP_PU2_SNOOP_PORT1_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[161], "BGP_PU2_SNOOP_PORT2_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[162], "BGP_PU2_SNOOP_PORT3_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[163], "BGP_PU2_SNOOP_PORT0_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[164], "BGP_PU2_SNOOP_PORT1_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[165], "BGP_PU2_SNOOP_PORT2_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[166], "BGP_PU2_SNOOP_PORT3_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[167], "BGP_PU2_SNOOP_PORT0_UPDATED_CACHE_LINE");
      strcpy(label[168], "BGP_PU2_SNOOP_PORT1_UPDATED_CACHE_LINE");
      strcpy(label[169], "BGP_PU2_SNOOP_PORT2_UPDATED_CACHE_LINE");
      strcpy(label[170], "BGP_PU2_SNOOP_PORT3_UPDATED_CACHE_LINE");
      strcpy(label[171], "BGP_PU2_SNOOP_PORT0_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[172], "BGP_PU2_SNOOP_PORT1_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[173], "BGP_PU2_SNOOP_PORT2_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[174], "BGP_PU2_SNOOP_PORT3_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[175], "BGP_PU3_SNOOP_PORT0_CACHE_REJECTED_REQUEST");
      strcpy(label[176], "BGP_PU3_SNOOP_PORT1_CACHE_REJECTED_REQUEST");
      strcpy(label[177], "BGP_PU3_SNOOP_PORT2_CACHE_REJECTED_REQUEST");
      strcpy(label[178], "BGP_PU3_SNOOP_PORT3_CACHE_REJECTED_REQUEST");
      strcpy(label[179], "BGP_PU3_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[180], "BGP_PU3_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[181], "BGP_PU3_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[182], "BGP_PU3_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_ACTIVE_SET");
      strcpy(label[183], "BGP_PU3_SNOOP_PORT0_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[184], "BGP_PU3_SNOOP_PORT1_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[185], "BGP_PU3_SNOOP_PORT2_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[186], "BGP_PU3_SNOOP_PORT3_HIT_STREAM_REGISTER_IN_HISTORY_SET");
      strcpy(label[187], "BGP_PU3_SNOOP_PORT0_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[188], "BGP_PU3_SNOOP_PORT1_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[189], "BGP_PU3_SNOOP_PORT2_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[190], "BGP_PU3_SNOOP_PORT3_STREAM_REGISTER_REJECTED_REQUEST");
      strcpy(label[191], "BGP_PU3_SNOOP_PORT0_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[192], "BGP_PU3_SNOOP_PORT1_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[193], "BGP_PU3_SNOOP_PORT2_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[194], "BGP_PU3_SNOOP_PORT3_RANGE_FILTER_REJECTED_REQUEST");
      strcpy(label[195], "BGP_PU3_SNOOP_PORT0_UPDATED_CACHE_LINE");
      strcpy(label[196], "BGP_PU3_SNOOP_PORT1_UPDATED_CACHE_LINE");
      strcpy(label[197], "BGP_PU3_SNOOP_PORT2_UPDATED_CACHE_LINE");
      strcpy(label[198], "BGP_PU3_SNOOP_PORT3_UPDATED_CACHE_LINE");
      strcpy(label[199], "BGP_PU3_SNOOP_PORT0_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[200], "BGP_PU3_SNOOP_PORT1_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[201], "BGP_PU3_SNOOP_PORT2_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[202], "BGP_PU3_SNOOP_PORT3_FILTERED_BY_CACHE_AND_REGISTERS");
      strcpy(label[203], "BGP_RESERVED");
      strcpy(label[204], "BGP_RESERVED");
      strcpy(label[205], "BGP_RESERVED");
      strcpy(label[206], "BGP_RESERVED");
      strcpy(label[207], "BGP_RESERVED");
      strcpy(label[208], "BGP_RESERVED");
      strcpy(label[209], "BGP_RESERVED");
      strcpy(label[210], "BGP_RESERVED");
      strcpy(label[211], "BGP_RESERVED");
      strcpy(label[212], "BGP_RESERVED");
      strcpy(label[213], "BGP_RESERVED");
      strcpy(label[214], "BGP_RESERVED");
      strcpy(label[215], "BGP_RESERVED");
      strcpy(label[216], "BGP_RESERVED");
      strcpy(label[217], "BGP_RESERVED");
      strcpy(label[218], "BGP_RESERVED");
      strcpy(label[219], "BGP_RESERVED");
      strcpy(label[220], "BGP_RESERVED");
      strcpy(label[221], "BGP_RESERVED");
      strcpy(label[222], "BGP_RESERVED");
      strcpy(label[223], "BGP_RESERVED");
      strcpy(label[224], "BGP_RESERVED");
      strcpy(label[225], "BGP_RESERVED");
      strcpy(label[226], "BGP_RESERVED");
      strcpy(label[227], "BGP_RESERVED");
      strcpy(label[228], "BGP_RESERVED");
      strcpy(label[229], "BGP_RESERVED");
      strcpy(label[230], "BGP_RESERVED");
      strcpy(label[231], "BGP_RESERVED");
      strcpy(label[232], "BGP_RESERVED");
      strcpy(label[233], "BGP_RESERVED");
      strcpy(label[234], "BGP_RESERVED");
      strcpy(label[235], "BGP_RESERVED");
      strcpy(label[236], "BGP_RESERVED");
      strcpy(label[237], "BGP_RESERVED");
      strcpy(label[238], "BGP_RESERVED");
      strcpy(label[239], "BGP_RESERVED");
      strcpy(label[240], "BGP_RESERVED");
      strcpy(label[241], "BGP_RESERVED");
      strcpy(label[242], "BGP_RESERVED");
      strcpy(label[243], "BGP_RESERVED");
      strcpy(label[244], "BGP_RESERVED");
      strcpy(label[245], "BGP_RESERVED");
      strcpy(label[246], "BGP_RESERVED");
      strcpy(label[247], "BGP_RESERVED");
      strcpy(label[248], "BGP_RESERVED");
      strcpy(label[249], "BGP_RESERVED");
      strcpy(label[250], "BGP_RESERVED");
      strcpy(label[251], "BGP_RESERVED");
      strcpy(label[252], "BGP_RESERVED");
      strcpy(label[253], "BGP_RESERVED");
      strcpy(label[254], "BGP_RESERVED");
      strcpy(label[255], "BGP_MISC_ELAPSED_TIME_UM3");
   }

}


