/*
 *	Author: Victor Junqiu WEI
 *	Email: wjqjsnj@gmail.com
 */

#ifndef DATA_UTILITY_H
#define	DATA_UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_struct.h"
#include "bst.h"
//#include "costenum.h"

#ifndef WIN32
#include<sys/resource.h>
#endif

#define TEXT_COL_MAX			1075
#define MAX_FILENAME_LENG		256			//maximum filename length.
#define	MAX_DESCRIPTION_LENG	256

#define	KEY_RANGE_SIZE			1
#define K_NN_PARAMETER			10

#define	CONFIG_FILE				"IRTree_config.txt"
#define	COSKQ_CONFIG_FILE		"config.txt"
#define	COSKQ_STAT_FILE			"nskq_stat.txt"
#define	COSKQ_RESULT_FILE		"nskq_result.txt"

#define GEN_IRTREE_CONFIG		"gen_irtree_config.txt"

#define SYN_CONFIG_FILE			"syn_config_key.txt"


//The structure for storing the frequency of the keywords.
//
typedef struct key_freq
{
	float*		freq;
	KEY_TYPE	key;
}	key_freq_t;


//The structure for storing the objects.
typedef struct 
{
	int			dim;
	int			obj_n;
	obj_t*		obj_v;

	int			key_n;
	//key_freq_t*	key_freq_v;
	bst_t*		key_freq_v;
}	data_t;


//The structure for storing the configuration info.
typedef	struct
{
//	int		M;
//	int		m;
	int		dim;
	//int	height_max;
	int		obj_n;
	int		key_n;

	int		split_opt;

	//char	data_file[ MAX_FILENAME_LENG];
	char	loc_file[ MAX_FILENAME_LENG];
	char	doc_file[ MAX_FILENAME_LENG];

	char	tree_file[ MAX_FILENAME_LENG];

	//IR-tree augmentation.

}	IRTree_config_t;

//The structure for storing the configuration information for the CoSKQ problem.
typedef struct coskq_config
{
	//Dataset specific.
	int		dim;
	int		obj_n;
	int		key_n;

	char	loc_file[ MAX_FILENAME_LENG];
	char	doc_file[ MAX_FILENAME_LENG];
	char	tree_file[ MAX_FILENAME_LENG];

	//IR-tree specific.
	//int		split_opt;

	//Cost measurement.
	int		cost_measure;

	//CoSKQ specific.
	int		alg_opt;
	int		prune_opt;

	int		q_key_n;
	int		q_set_size;

	//range	key_freq_range;
	int		low;
	int		high;

	//Simp-Appro.
	float	diag;

	int		tree_tag;

	//ratio.
	float	ratio_thr;

	//CostEnum specific.
}	coskq_config_t;

//The structure of the config file for generating synthetic data.
typedef struct syn_config
{

	int		dim;
	
	//for base data.
	int		obj_n;
	int		key_n;
	char	loc_file[ MAX_FILENAME_LENG];
	char	doc_file[ MAX_FILENAME_LENG];

	//for new data.
	int		new_obj_n;
	char	new_loc_file[ MAX_FILENAME_LENG];
	char	new_doc_file[ MAX_FILENAME_LENG];
	
	float	dev;

}	syn_config_t;

//
typedef struct key_freq_pair
{
	int		key;
	int		freq;
}	key_freq_pair_t;


extern coskq_stat_t stat_v;

//extern IRTree_t IRTree_v;

#ifndef WIN32

void GetCurTime( rusage* curTime);

void GetTime( struct rusage* timeStart, struct rusage* timeEnd, float* userTime, float* sysTime);

#endif

IRTree_config_t* read_config_irtree( );

coskq_config_t* read_config_coskq( );

void add_keyword_entry( k_node_t* &k_node_v, KEY_TYPE key);

void copy_k_list( k_node_t* k_head1, k_node_t* k_head2);

void print_k_list( k_node_t* k_head, FILE* o_fp);

void print_coskq_stat( coskq_config_t* cfg, int cnt);

void alloc_obj( obj_t* obj_v, int dim);

data_t*	read_data_irtree( IRTree_config_t* cfg);

data_t* alloc_data( int num);

data_t*	read_data_coskq( coskq_config_t* cfg);

void release_k_list( k_node_t* k_node_v);

void release_data( data_t* data_v);


//Query & dataset generation utilities.
int rand_i( int min, int max);

float rand_f( float min_v, float max_v);

bool is_old( int* rand_v, int cur, int var);


//Generate synthetic datasets.
float gaussian_f( float mean, float s_var);

range* collect_data_range( data_t* data_v);

void print_data( data_t* data_v, syn_config_t* s_cfg);

void gen_syn_data( syn_config_t* s_cfg);

void gen_syn_data_keyword( syn_config_t* s_cfg);

void batch_gen_syn_data( );



#endif




