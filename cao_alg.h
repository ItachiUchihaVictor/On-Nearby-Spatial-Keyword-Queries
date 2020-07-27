
/*
 *	Author: Victor Junqiu WEI
 *	Author: wjqjsnj@gmail.com
 */

/*
 *	The implementations of the algorithms in "Collective Spatial Keyword Querying"
 *	by Xin Cao et al. 
 */

#ifndef CAO_ALG_H
#define CAO_ALG_H

#include "data_struct2.h"
#include "costenum.h"
#include "b_heap.h"
#include "bst.h"

extern	int	cost_tag;

cns_t* alloc_cns( );

cns_t* copy_cns( cns_t* cns_v);

bool has_same_content_cns( cns_t* cns_v1, cns_t* cns_v2); 

void add_cns_entry( cns_t* cns_v, node_t* node_v);

void print_cns( cns_t* cns_v);

void release_cns( cns_t* cns_v);

cns_list_t* alloc_cns_list( );

cns_list_t* copy_cns_list( cns_list_t* cns_list_v);

void add_cns_list_entry( cns_list_t* cns_list_v, cns_t* cns_v);

void print_cns_list( cns_list_t* cns_list_v);

void release_cns_list( cns_list_t* cns_list_v);

cns_list_set_t* alloc_cns_list_set( );

void release_cns_list_set( cns_list_set_t* cns_list_set_v); 

obj_set_list_t* alloc_obj_set_list( );

void add_obj_set_list_entry( obj_set_list_t* obj_set_list_v, obj_set_t* obj_set_v);

void release_obj_set_list( obj_set_list_t* obj_set_list_v);

B_KEY_TYPE calc_minDist_node( node_t* node_v1, node_t* node_v2);

B_KEY_TYPE calc_low_cost( cns_t* cns_v, query_t* q);

obj_set_t* Cao_Appro1( query_t* q);

bool is_contained_query( query_t* q, KEY_TYPE key);

bool is_contained_obj_set( obj_set_t* obj_set_v, KEY_TYPE key);

KEY_TYPE find_critical_keyword( obj_set_t* obj_set_v, query_t* q);

obj_set_t* Cao_Appro2( query_t* q);

obj_t* join_check_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2);

//obj_set_t* merge_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2);

void reformat_obj_set_list( obj_set_list_t* obj_set_list_v);

bool is_existing_obj_set( obj_set_list_t* obj_set_list_v, obj_set_t* obj_set_v);

obj_set_list_t* get_obj_set_list_apriori( node_t* leaf_node_v, B_KEY_TYPE cost, query_t* q, int len_bound);

obj_set_list_t* copy_obj_set_list( obj_set_list_t* obj_set_list_v);

void append_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2); 

obj_set_list_t* exhaustive_node_set_search_sub( obj_set_list_t** oList, int sta, int end, query_t* q);

obj_set_t* exhaustive_node_set_search( cns_t* cns_v, B_KEY_TYPE cost, query_t* q);

cns_list_t* const_L1_apriori( node_t* node_v, query_t* q, B_KEY_TYPE cost);

node_t* join_check( cns_t* cns_v1, cns_t* cns_v2);

//bool is_contained_node( cns_t* cns_v, node_t* node_v);

//cns_t* merge_cns( cns_t* cns_v1, cns_t* cns_v2);

cns_list_t* reformat_cns_list_set( cns_list_set_t* cns_list_set_v);

void append_cns( cns_t* cns_v1, cns_t* cns_v2);

bool is_existing_cns( cns_list_t* cns_list_v, cns_t* cns_v);

cns_list_t* enumerate_node_set_sub( cns_list_t** cList, int sta, int end, query_t* q);

cns_list_t* enumerate_node_set( cns_t* cns_v, B_KEY_TYPE cost, query_t* q);

bool is_covered_cns( cns_t* cns_v, query_t* q);

obj_set_t* Cao_Exact( query_t* q);


#endif