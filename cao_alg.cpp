
/*
 *	Author: Victor Junqiu WEI
 *	Author: wjqjsnj@gmail.com
 */

#include "cao_alg.h"

/*
 *	Allocate a cns_t structure.
 */
cns_t* alloc_cns( )
{
	cns_t* cns_v;

	cns_v = ( cns_t*)malloc( sizeof( cns_t));
	memset( cns_v, 0, sizeof( cns_t));

	cns_v->list_head = ( node_list_t*)malloc( sizeof( node_list_t));
	memset( cns_v->list_head, 0, sizeof( node_list_t));

	/*s*/
	stat_v.memory_v += sizeof( cns_t) + sizeof( node_list_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return cns_v;
}

/*
 *	Copy a cns_t structure @cns_v.
 */
cns_t* copy_cns( cns_t* cns_v)
{
	cns_t* rtn;
	node_list_t* node_list_iter;

	rtn = alloc_cns( );

	node_list_iter = cns_v->list_head->next;
	while( node_list_iter != NULL)
	{
		add_cns_entry( rtn, node_list_iter->node_v);

		node_list_iter = node_list_iter->next;
	}
	
	return rtn;
}

/*
 *	Check whether two cns_t structures @cns_v1 and @cns_v2 have the same content.
 */
bool has_same_content_cns( cns_t* cns_v1, cns_t* cns_v2)
{
	node_list_t* iter1, *iter2;

	if( cns_v1->node_n != cns_v2->node_n)
		return false;

	iter1 = cns_v1->list_head->next;
	while( iter1 != NULL)
	{
		iter2 = cns_v2->list_head->next;
		while( iter2 != NULL)
		{
			if( iter2->node_v == iter1->node_v)
				break;

			iter2 = iter2->next;
		}

		if( iter2 == NULL)
			return false;

		iter1 = iter1->next;
	}

	return true;
}

/*
 *	Add an node @node_v into the cns_t structure @cns_v.
 */
void add_cns_entry( cns_t* cns_v, node_t* node_v)
{
	node_list_t* node_list_v;

	node_list_v = ( node_list_t*)malloc( sizeof( node_list_t));
	memset( node_list_v, 0, sizeof( node_list_t));

	node_list_v->node_v =  node_v;
	node_list_v->next =  cns_v->list_head->next;
	cns_v->list_head->next = node_list_v;

	cns_v->node_n ++;

	/*s*/
	stat_v.memory_v += sizeof( node_list_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 *	Print a cns_t structure.
 */
void print_cns( cns_t* cns_v)
{
	node_list_t* iter;
	
	iter = cns_v->list_head->next;
	while( iter != NULL)
	{
		printf( "%Np  ", iter->node_v);

		iter = iter->next;
	}
	printf( "\n");
}

/*
 *	Release a cns_t structure.
 */
void release_cns( cns_t* cns_v)
{
	node_list_t* node_list_v1, *node_list_v2;

	node_list_v1 = cns_v->list_head;
	while( node_list_v1->next != NULL)
	{
		node_list_v2 = node_list_v1->next;
		free( node_list_v1);
		
		node_list_v1 = node_list_v2;
	}
	free( node_list_v1);

	/*s*/
	stat_v.memory_v -= ( cns_v->node_n + 1) * sizeof( node_list_t) + sizeof( cns_t);
	/*s*/

	free( cns_v);
}

/*
 *	Alloc a cns_list_t structure.
 */
cns_list_t* alloc_cns_list( )
{
	cns_list_t* cns_list_v;

	cns_list_v = ( cns_list_t*)malloc( sizeof( cns_list_t));
	memset( cns_list_v, 0, sizeof( cns_list_t));

	cns_list_v->head = alloc_cns( );

	/*s*/
	stat_v.memory_v += sizeof( cns_list_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return cns_list_v;
}

/*
 *	Copy a cns_list_v structure @cns_list_v.
 */
cns_list_t* copy_cns_list( cns_list_t* cns_list_v)
{
	cns_list_t* rtn;
	cns_t* cns_v, *cns_iter;

	rtn = alloc_cns_list( );
	cns_v = rtn->head;

	cns_iter = cns_list_v->head->next;
	while( cns_iter != NULL)
	{
		cns_v->next = copy_cns( cns_iter);
		cns_v = cns_v->next;

		cns_iter = cns_iter->next;
	}

	return rtn;
}

/*
 *	Print a cns_list_t structure.
 */
void print_cns_list( cns_list_t* cns_list_v)
{
	cns_t* iter;

	iter = cns_list_v->head->next;
	while( iter != NULL)
	{
		print_cns( iter);

		iter = iter->next;
	}
	printf( "\n");
}

/*
 *	Add an cns_t entry @cns_v into the cns_list_t structure @cns_list_v.
 */
void add_cns_list_entry( cns_list_t* cns_list_v, cns_t* cns_v)
{
	cns_v->next = cns_list_v->head->next;
	cns_list_v->head->next = cns_v;

	cns_list_v->cns_n ++;
}

/*
 *	Release a cns_list_t structure.
 */
void release_cns_list( cns_list_t* cns_list_v)
{
	cns_t* cns_v1, *cns_v2;

	cns_v1 = cns_list_v->head;
	while( cns_v1->next != NULL)
	{
		cns_v2 = cns_v1->next;
		release_cns( cns_v1);
		
		cns_v1 = cns_v2;
	}
	release_cns( cns_v1);

	free( cns_list_v);

	/*s*/
	stat_v.memory_v -= sizeof( cns_list_t);
	/*s*/
}

/*
 *	Alloc a cns_list_set_t structure.
 */
cns_list_set_t* alloc_cns_list_set( )
{
	cns_list_set_t* cns_list_set_v;

	cns_list_set_v = ( cns_list_set_t*)malloc( sizeof( cns_list_set_t));
	memset( cns_list_set_v, 0, sizeof( cns_list_set_t));

	cns_list_set_v->head = alloc_cns_list( );

	/*s*/
	stat_v.memory_v += sizeof( cns_list_set_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return cns_list_set_v;
}

/*
 *	Release a cns_list_set_t structure.
 */
void release_cns_list_set( cns_list_set_t* cns_list_set_v)
{
	cns_list_t* cns_list_v1, *cns_list_v2;

	cns_list_v1 = cns_list_set_v->head;
	while( cns_list_v1->next != NULL)
	{
		cns_list_v2 = cns_list_v1->next;
		release_cns_list( cns_list_v1);

		cns_list_v1 = cns_list_v2;
	}
	release_cns_list( cns_list_v1);

	free( cns_list_set_v);

	/*s*/
	stat_v.memory_v -= sizeof( cns_list_set_t);
	/*s*/
}

/*
 *	Allocate an obj_set_list_t structure.
 */
obj_set_list_t* alloc_obj_set_list( )
{
	obj_set_list_t* obj_set_list_v;

	obj_set_list_v = ( obj_set_list_t*)malloc( sizeof( obj_set_list_t));
	memset( obj_set_list_v, 0, sizeof( obj_set_list_t));

	/*s*/
	stat_v.memory_v += sizeof( obj_set_list_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return obj_set_list_v;
}

/*
 *	Add an entry @obj_set_v into @obj_set_list_v.
 */
void add_obj_set_list_entry( obj_set_list_t* obj_set_list_v, obj_set_t* obj_set_v)
{
	obj_set_list_t* tmp;

	tmp = ( obj_set_list_t*)malloc( sizeof( obj_set_list_t));
	memset( tmp, 0, sizeof( obj_set_list_t));

	tmp->obj_set_v = obj_set_v;
	tmp->next = obj_set_list_v->next;
	obj_set_list_v->next = tmp;

	obj_set_list_v->obj_set_n ++;

	/*s*/
	stat_v.memory_v += sizeof( obj_set_list_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 * Print an obj_set_list_t structure @obj_set_list_v.
 */
void print_obj_set_list( obj_set_list_t* obj_set_list_v)
{
	obj_set_list_t* iter1, *iter2;

	iter1 = obj_set_list_v;
	while( iter1 != NULL)
	{
		iter2 = iter1->next;
		while( iter2 != NULL)
		{
			print_obj_set( iter2->obj_set_v, stdout);
			printf( "\n");

			iter2 = iter2->next;
		}

		iter1 = iter1->down;
	}
}

/*
 * Release an obj_set_list_t structure.
 */
void release_obj_set_list( obj_set_list_t* obj_set_list_v)
{
	obj_set_list_t* iter1, *iter2, *tmp;

	iter1 = obj_set_list_v;
	while( iter1 != NULL)
	{
		iter2 = iter1;
		iter1 = iter1->down;
		while( iter2->next != NULL)
		{
			tmp = iter2->next;
			
			if( iter2->obj_set_v != NULL)
				release_obj_set( iter2->obj_set_v);
			free( iter2);

			iter2 = tmp;

			/*s*/
			stat_v.memory_v -= sizeof( obj_set_list_t);
			/*s*/
		}
		if( iter2->obj_set_v != NULL)
			release_obj_set( iter2->obj_set_v);
		free( iter2);

		/*s*/
		stat_v.memory_v -= sizeof( obj_set_list_t);
		/*s*/
	}
}

/*
 *	Calculate the minimum distance between the MBRs of two nodes @node_v1 and @node_v2.
 */
B_KEY_TYPE calc_minDist_node( node_t* node_v1, node_t* node_v2)
{
	int i;
	float diff;
	B_KEY_TYPE minDist;
	range* MBR1, *MBR2;

	if( node_v1->parent == NULL)
		MBR1 = get_MBR_node( node_v1, IRTree_v.dim);
	else
		MBR1 = node_v1->parent->MBRs[ node_v1->loc];

	if( node_v2->parent == NULL)
		MBR2 = get_MBR_node( node_v2, IRTree_v.dim);
	else
		MBR2 = node_v2->parent->MBRs[ node_v2->loc];

	minDist = 0;
	for( i=0; i<IRTree_v.dim; i++)
	{
		if( MBR1[ i].max < MBR2[ i].min)
			diff = MBR2[ i].min - MBR1[ i].max;
		else if( MBR1[ i].min > MBR2[ i].max)
			diff = MBR1[ i].min - MBR2[ i].max;
		else
			diff = 0;

		minDist += pow( diff, 2);
	}

	if( node_v1->parent == NULL)
	{
		free( MBR1);
		
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof( range);
		/*s*/
	}
	if( node_v2->parent == NULL)
	{
		free( MBR2);
	
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof( range);
		/*s*/
	}

	return sqrt( minDist);		
}

/*
 *	Calculate the lower bound of the cost of a covering node set.
 */
B_KEY_TYPE calc_low_cost( cns_t* cns_v, query_t* q)
{
	B_KEY_TYPE comp1, comp2, tmp;
	node_list_t* node_list_iter, *node_list_iter2;
	node_t* node_v;
	range* MBR;

	//Compute the first component.
	comp1 = 0;
	node_list_iter = cns_v->list_head->next;
	while( node_list_iter != NULL)
	{
		node_v = node_list_iter->node_v;
		if( node_v->parent == NULL)
			MBR = get_MBR_node( node_v, IRTree_v.dim);
		else
			MBR = node_v->parent->MBRs[ node_v->loc];

		tmp = calc_minDist( MBR, q->loc_v);
		if( tmp > comp1)
			comp1 = tmp;

		if( node_v->parent == NULL)
		{		
			free( MBR);

			/*s*/
			stat_v.memory_v -= IRTree_v.dim * sizeof( range);
			/*s*/
		}

		node_list_iter = node_list_iter->next;
	}

	//Compute the second component.
	comp2 = 0;
	node_list_iter = cns_v->list_head->next;
	while( node_list_iter->next != NULL)
	{
		node_list_iter2 = node_list_iter->next;
		while( node_list_iter2 != NULL)
		{
			tmp = calc_minDist_node( node_list_iter->node_v, node_list_iter2->node_v);

			if( tmp > comp2)
				comp2 = tmp;

			node_list_iter2 = node_list_iter2->next;
		}

		node_list_iter = node_list_iter->next;
	}

	if( cost_tag == 1)
		return comp1 + comp2;
	else
		return max( comp1, comp2);
}

/*
 *	The implementation of the "Cao-Appro1" algrithm.
 *
 *	return NULL if no feasible set is possible, 
 *	i.e., the keywords in @q could not be covered,
 *	otherwise, return the approximate solution.
 */
obj_set_t* Cao_Appro1( query_t* q)
{
	obj_set_t* S;
	obj_t* obj_v;
	k_node_t* k_node_v;

	//printf( "Cao-Appro1:\n");

	S = alloc_obj_set( );

	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		obj_v = const_NN_key( q->loc_v, k_node_v->key, NULL);

		if( obj_v == NULL)
		{
			printf( "No solution exists!\n");
			release_obj_set( S);
			return NULL;
		}

		add_obj_set_entry( obj_v, S);

		k_node_v = k_node_v->next;
	}

	return S;
}

/*
 *	Check whether a keyword @key is contained in the query @q.
 */
bool is_contained_query( query_t* q, KEY_TYPE key)
{
	k_node_t* k_node_iter;

	k_node_iter = q->psi_v->k_head->next;
	while( k_node_iter != NULL)
	{
		if( k_node_iter->key == key)
			return true;

		k_node_iter = k_node_iter->next;
	}
	
	return false;
}

/*
 *	Check wheher the objects in @obj_set_v contain a keyword @key.
 */
bool is_contained_obj_set( obj_set_t* obj_set_v, KEY_TYPE key)
{
	obj_node_t* obj_node_iter;

	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		if( has_key_obj( obj_node_iter->obj_v, key))
			return true;

		obj_node_iter = obj_node_iter->next;
	}

	return false;
}

/*
 *	Find the keyword that is only contained in the fartheset obj
 *	in @obj_set_v from @q.
 */
KEY_TYPE find_critical_keyword( obj_set_t* obj_set_v, query_t* q)
{
	KEY_TYPE c_key;
	B_KEY_TYPE max_dist, dist;
	obj_t* obj_v;
	obj_node_t* obj_node_v, *obj_node_pre, *obj_node_target_pre, *obj_node_target;
	k_node_t* k_node_v;

	max_dist = 0;
	obj_v = NULL;
	obj_node_target_pre = NULL;

	//bug.
	//Remove identical obj_v in obj_set_v.
	remove_identical_obj( obj_set_v);

	//Figure out the obj that is the farthest from @q.
	obj_node_pre = obj_set_v->head;
	obj_node_v = obj_node_pre->next;
	while( obj_node_v != NULL)
	{
		dist = calc_minDist( obj_node_v->obj_v->MBR, q->loc_v);
		if( dist > max_dist)
		{
			max_dist = dist;
			obj_v = obj_node_v->obj_v;
			obj_node_target_pre = obj_node_pre;
		}

		obj_node_pre = obj_node_v;
		obj_node_v = obj_node_pre->next;
	}

	obj_node_target = obj_node_target_pre->next;

	//Remove obj_node_target from obj_set_v temporarily [for checking use].
	obj_node_target_pre->next = obj_node_target->next;
	obj_set_v->obj_n --;

	//Traverse all the keywords of obj_v to find the critical keyword c_key.
	k_node_v = obj_v->k_head->next;
	while( k_node_v != NULL)
	{
		c_key = k_node_v->key;

		if( is_contained_query( q, c_key) &&
			!is_contained_obj_set( obj_set_v, c_key))
			break;

		k_node_v = k_node_v->next;
	}

/*t*/
	if( k_node_v == NULL)
	{
		fprintf( stderr, "Inconsistency in find_critical_keyword.\n");
		exit( 0);		
	}
/*t*/

	//Include back the obj_node_target.
	obj_node_target->next = obj_node_target_pre->next;
	obj_node_target_pre->next = obj_node_target;
	obj_set_v->obj_n ++;

	return c_key;
}

/*
 *	The implementation of the "Cao-Appro2" algorithm.
 *
 *	Online "k-NN" process.
 */
obj_set_t* Cao_Appro2( query_t* q)
{
	int i, rear, top;
	KEY_TYPE c_key;
	B_KEY_TYPE costV, costV_1, dist;
	b_heap_t* U;
	obj_set_t* V, *V_1, *V_tmp;
	void* e;
	node_t* node_v;
	obj_t* obj_v;
	bst_node_t* bst_node_v;
	BIT_TYPE p_list;
	loc_t* loc_v;
	query_t* q_new;



	//printf( "Cao-Appro2:\n");

	U = alloc_b_heap( INI_HEAP_SIZE);

	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);

	b_h_insert( U, rear++);
	
	V = Cao_Appro1( q);	

	/*s*/
	stat_v.n_k ++;
	/*s*/

	if( V == NULL)
	{
		release_b_heap( U);
		return NULL;
	}

	costV = comp_cost( V, q);

	//Find the keyword only contained by the object that is farthest from @q in V.
	c_key = find_critical_keyword( V, q);

	//Best-first search process.
	while( !b_h_is_empty( U))
	{
		top = b_h_get_top( U);
		e = U->obj_arr[ top].element;

		/*t/
			if( e == NULL)
			printf( "");
		/*t*/

		if( U->obj_arr[ top].e_tag == 1)
		{
			//e is an node_t*.
			node_v = ( node_t*)e;

			//Check the distance.
			dist = calc_minDist_node( node_v, q->loc_v);

			bst_node_v = bst_search( node_v->bst_v, c_key);
			if( bst_node_v == NULL)
				continue;
		
			p_list = bst_node_v->p_list;
			for( i=0; i<node_v->num; i++)
			{
				//Check the c_key keyword.
				if( !get_k_bit( p_list, i))
					continue;

				/*t/
				if( node_v->child[ i] == NULL)
					printf( "");
				/*t*/

				/*t/
				if( rear == 9999)
					printf( "");
				/*t*/

				U->obj_arr[ rear].element = node_v->child[ i];
				if( node_v ->level > 0)
				{
					//node_v is an inner-node.
					U->obj_arr[ rear].e_tag = 1;
					U->obj_arr[ rear].key = calc_minDist( node_v->MBRs[ i], q->loc_v);
				}
				else
				{
					//node_v is a leaf-node.
					U->obj_arr[ rear].e_tag = 2;
					U->obj_arr[ rear].key = calc_minDist( node_v->MBRs[ i], q->loc_v);
				}

				//Enqueue.
				b_h_insert( U, rear++);
			}//
		}
		else
		{
			//e is an obj_t*.
			obj_v = ( obj_t*)e;
			loc_v = get_obj_loc( obj_v);

			//Construct a new query instance.
			q_new = alloc_query( );
			q_new->loc_v = loc_v;

			q_new->psi_v = alloc_psi( );
			copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			V_1 = Cao_Appro1( q_new);
			//add_obj_set_entry( obj_v, V_1);

			if( V_1 == NULL)
			{
				release_query( q_new);
				continue;
			}

			costV_1 = comp_cost( V_1, q);

			if( costV_1 < costV)
			{
				costV = costV_1;
				V_tmp = V;
				V = V_1;
				
				release_obj_set( V_tmp);
			}
			else
				release_obj_set( V_1);

			release_query( q_new);		
			
			/*s*/
			stat_v.n_k ++;
			/*s*/
		}//else
	}//while

	release_b_heap( U);

	return V;
}

/*
 *	[Apriori] The join check operation over two sets of objects @obj_set_v1 and @obj_set_v2.
 */
obj_t* join_check_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	int cnt;
	obj_t* obj_v;
	obj_node_t* obj_node_iter1, *obj_node_iter2;

	obj_v = NULL;
	cnt = 0;
	obj_node_iter1 = obj_set_v1->head->next;
	while( obj_node_iter1 != NULL)
	{
		obj_node_iter2 = obj_set_v2->head->next;
		while( obj_node_iter2 != NULL)
		{
			if( obj_node_iter2->obj_v == obj_node_iter1->obj_v)
				break;

			obj_node_iter2 = obj_node_iter2->next;
		}

		if( obj_node_iter2 == NULL)
		{
			obj_v = obj_node_iter1->obj_v;
			cnt ++;
			if( cnt >= 2)
				return NULL;
		}

		obj_node_iter1 = obj_node_iter1->next;
	}

	return obj_v;
}

/*
 *	Linerize the cross link @obj_set_list_v.
 *	From "cross list" to a "single list".
 */
void reformat_obj_set_list( obj_set_list_t* obj_set_list_v)
{
	obj_set_list_t* iter1, *iter2, *tmp;

	iter1 = obj_set_list_v->down;
	obj_set_list_v->down = NULL;
	while( iter1 != NULL)
	{
		iter2 = iter1;
		while( iter2->next != NULL)
			iter2 = iter2->next;

		iter2->next = obj_set_list_v->next;
		obj_set_list_v->next = iter1->next;
		obj_set_list_v->obj_set_n += iter1->obj_set_n;

		tmp = iter1;
		iter1 = iter1->down;
		free( tmp);

		/*s*/
		stat_v.memory_v -= sizeof( obj_set_list_t);
		/*s*/
	}
}

/*
 *	Check whether an obj_set_t structure obj_set_v exists
 *	in a obj_set_list_t structure @obj_set_list_v.
 */
bool is_existing_obj_set( obj_set_list_t* obj_set_list_v, obj_set_t* obj_set_v)
{
	obj_set_list_t* iter;

	iter = obj_set_list_v->next;
	while( iter != NULL)
	{
		if( has_same_content_obj_set( iter->obj_set_v, obj_set_v))
			return true;

		iter = iter->next;
	}

	return false;
}

/*
 *	[Apriori] Compute the subsets of objects under a leaf node @leaf_node_v.
 *
 *	@cost is the best-known cost.
 *	@q is the query.
 *	@len_bound bounds the number of objects in the subsets.
 */
obj_set_list_t* get_obj_set_list_apriori( node_t* leaf_node_v, B_KEY_TYPE cost, query_t* q, int len_bound)
{
	int i, cnt;
	B_KEY_TYPE low_cost;
	obj_set_t* obj_set_v, *obj_set_v1, *obj_set_v2;
	obj_set_list_t* obj_set_list_v, *L_1, *L_pre, *L_m, *obj_set_list_v1, *obj_set_list_v2;;
	obj_t* obj_v;

	obj_set_list_v = alloc_obj_set_list( );

	//Construct L1.
	L_1 = obj_set_list_v;
	for( i=0; i<leaf_node_v->num; i++)
	{
		obj_v = ( obj_t*)( leaf_node_v->child[ i]);
		
		low_cost = calc_minDist( obj_v->MBR, q->loc_v);
		if( low_cost < cost && is_relevant_obj( obj_v, q))
		{
			//Add obj_v into L1.
			obj_set_v = alloc_obj_set( );
			add_obj_set_entry( obj_v, obj_set_v);

			add_obj_set_list_entry( L_1, obj_set_v);
		}
	}

	//Construct L2, L3, ..., Lm.
	cnt = 1;
	L_pre = L_1;
	while( L_pre->obj_set_n >= 2 && cnt < len_bound)
	{
		L_m = alloc_obj_set_list( );

		obj_set_list_v1 = L_pre->next;
		while( obj_set_list_v1->next != NULL)
		{
			obj_set_v1 = obj_set_list_v1->obj_set_v;

			obj_set_list_v2 = obj_set_list_v1->next;
			while( obj_set_list_v2 != NULL)
			{
				obj_set_v2 = obj_set_list_v2->obj_set_v;

				//Join check on obj_set_v1 and obj_set_v2.
				if( ( obj_v = join_check_obj_set( obj_set_v1, obj_set_v2)))
				{
					//Merge.
					//obj_set_v = merge_obj_set( obj_set_v1, obj_set_v2);
					obj_set_v = copy_obj_set( obj_set_v2);

					add_obj_set_entry( obj_v, obj_set_v);
					
					//
					if( comp_cost( obj_set_v, q) < cost &&
						!is_existing_obj_set( L_m, obj_set_v))
						add_obj_set_list_entry( L_m, obj_set_v);
					else
						release_obj_set( obj_set_v);
				}

				obj_set_list_v2 = obj_set_list_v2->next;
			}

			obj_set_list_v1 = obj_set_list_v1->next;
		}

		L_pre->down = L_m;
		L_pre = L_m;
		cnt ++;
	}

	reformat_obj_set_list( obj_set_list_v);

	return obj_set_list_v;
}

/*
 *	Copy an obj_set_list_t structure @obj_set_list_v.
 */
obj_set_list_t* copy_obj_set_list( obj_set_list_t* obj_set_list_v)
{	
	obj_set_list_t* rtn, *tmp, *obj_set_list_iter;

	rtn = alloc_obj_set_list( );

	obj_set_list_iter = obj_set_list_v->next;
	while( obj_set_list_iter != NULL)
	{
		tmp = ( obj_set_list_t*)malloc( sizeof( obj_set_list_t));
		memset( tmp, 0, sizeof( obj_set_list_t));

		tmp->obj_set_v = obj_set_list_iter->obj_set_v;
		tmp->next = rtn->next;
		rtn->next = tmp;

		obj_set_list_iter = obj_set_list_iter->next;

		/*s*/
		stat_v.memory_v += sizeof( obj_set_list_t);
		if( stat_v.memory_v > stat_v.memory_max)
			stat_v.memory_max = stat_v.memory_v;
		/*s*/
	}

	rtn->obj_set_n = obj_set_list_v->obj_set_n;
	
	return rtn;
}

/*
 *	Append the objs in @obj_set_v2 to @obj_set_v1.
 *
 *	In this implementation, @obj_set_v2 becomes empty after this process.
 */
void append_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	obj_node_t* obj_node_iter;

	obj_node_iter = obj_set_v2->head->next;
	while( obj_node_iter != NULL)
	{
		add_obj_set_entry( obj_node_iter->obj_v, obj_set_v1);

		obj_node_iter = obj_node_iter->next;
	}

	/*
	obj_node_iter = obj_set_v1->head;
	while( obj_node_iter->next != NULL)
		obj_node_iter = obj_node_iter->next;

	obj_node_iter->next = obj_set_v2->head->next;
	obj_set_v1->obj_n += obj_set_v2->obj_n;

	obj_set_v2->head->next = NULL;
	obj_set_v2->obj_n = 0;
	*/
}

/*
 *	[Recursive] The sub-procedure of the "exhaustive_node_set_search" function.
 */
obj_set_list_t* exhaustive_node_set_search_sub( obj_set_list_t** oList, int sta, int end, query_t* q)
{
	obj_set_list_t *obj_set_list_v, *obj_set_list_iter1, *obj_set_list_iter2, *obj_set_list_v1;
	obj_set_t *obj_set_v, *obj_set_v1, *obj_set_v2;

	//Base step.
	if( sta == end)
		return copy_obj_set_list( oList[ sta]);

	//Recursive step.
	obj_set_list_v = alloc_obj_set_list( );

	obj_set_list_v1 = exhaustive_node_set_search_sub( oList, sta+1, end, q);

	obj_set_list_iter1 = oList[ sta]->next;
	while( obj_set_list_iter1 != NULL)
	{
		obj_set_v1 = obj_set_list_iter1->obj_set_v;

		/*t/
		print_obj_set( obj_set_v1, stdout);
		/*t*/

		obj_set_list_iter2 = obj_set_list_v1->next;
		while( obj_set_list_iter2 != NULL)
		{
			obj_set_v2 = obj_set_list_iter2->obj_set_v;

			/*t/
			print_obj_set( obj_set_v2, stdout);
			/*t*/

			//Alloc a new obj_set_t structure.
			obj_set_v = alloc_obj_set( );

			//Append the objs in obj_set_v1 to obj_set_v.
			append_obj_set( obj_set_v, obj_set_v1);

			//Append the objs in obj_set_v2 to obj_set_v.	
			append_obj_set( obj_set_v, obj_set_v2);

			/*t/
			print_obj_set( obj_set_v, stdout);
			/*t*/
		
			//Add obj_set_v into obj_set_list_v.
			//if( is_covered_obj_set( obj_set_v, q))
			add_obj_set_list_entry( obj_set_list_v, obj_set_v);

			obj_set_list_iter2 = obj_set_list_iter2->next;
		}

		obj_set_list_iter1 = obj_set_list_iter1->next;
	}

	//Release the resources.
	release_obj_set_list( obj_set_list_v1);

	return obj_set_list_v;
}

/*
 *	Construct the feasible set with the smallest cost wrt @q
 *	from a node set @cns_v.
 */
obj_set_t* exhaustive_node_set_search( cns_t* cns_v, B_KEY_TYPE cost, query_t* q)
{
	int i, node_n, len_bound;
	B_KEY_TYPE min_cost, cost_tmp;
	obj_set_t* rtn;
	node_list_t* node_list_iter;
	node_t* leaf_node_v;
	obj_set_list_t* obj_set_list_v, *obj_set_list_iter1, *obj_set_list_iter2;
	obj_set_list_t* obj_set_list_best, *obj_set_list_best_pre;
	obj_set_list_t** oList;

	node_n = cns_v->node_n;

	oList = ( obj_set_list_t**)malloc( node_n * sizeof( obj_set_list_t*));
	memset( oList, 0, node_n * sizeof( obj_set_list_t*));

	/*s*/
	stat_v.memory_v += node_n * sizeof( obj_set_list_t*);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	len_bound = q->psi_v->key_n - node_n + 1;
	node_list_iter = cns_v->list_head->next;
	for( i=0; i<node_n; i++)
	{
		leaf_node_v = node_list_iter->node_v;

		oList[ i] = get_obj_set_list_apriori( leaf_node_v, cost, q, len_bound); 

		/*t/
		print_obj_set_list( oList[ i]);
		/*t*/

		node_list_iter = node_list_iter->next;
	}

	//Construct the set of obj_set_t structure by picking one obj_set_t structure
	//from oList[ i] for i=0, 1, ..., node_n-1.	
	obj_set_list_v = exhaustive_node_set_search_sub( oList, 0, node_n-1, q);

	/*t/
	print_obj_set_list( obj_set_list_v);
	/*t*/

	//Find the obj_set_t structure with the smallest cost wrt @q within obj_set_list_v.
	min_cost = DBL_MAX;
	obj_set_list_best = NULL;
	obj_set_list_best_pre = NULL;
	
	obj_set_list_iter1 = obj_set_list_v;
	obj_set_list_iter2 = obj_set_list_iter1->next;
	while( obj_set_list_iter2 != NULL)
	{
		/*s*/
		stat_v.feasible_set_n ++;
		/*s*/

		cost_tmp = comp_cost( obj_set_list_iter2->obj_set_v, q);
		if( cost_tmp < min_cost && 
			is_covered_obj_set( obj_set_list_iter2->obj_set_v, q))
		{
			/*t/
			print_obj_set( obj_set_list_iter2->obj_set_v, stdout);
			/*t*/

			min_cost = cost_tmp;
			obj_set_list_best = obj_set_list_iter2;
			obj_set_list_best_pre = obj_set_list_iter1;
		}

		obj_set_list_iter1 = obj_set_list_iter2;
		obj_set_list_iter2 = obj_set_list_iter2->next;
	}

	if( obj_set_list_best_pre)
	{
		rtn = obj_set_list_best->obj_set_v;

		obj_set_list_best_pre->next = obj_set_list_best->next;
		obj_set_list_v->obj_set_n --;
	}
	else
		rtn = NULL;

	//Release the resources.
	for( i=0; i<node_n-1; i++)
		release_obj_set_list( oList[ i]);
	free( oList);

	release_obj_set_list( obj_set_list_v);

	if( rtn != NULL)
	{
		free( obj_set_list_best);

		/*s*/
		stat_v.memory_v -=  sizeof( obj_set_list_t);
		/*s*/
	}

	/*s*/
	stat_v.memory_v -= node_n * sizeof( obj_set_list_t*);
	/*s*/

	return rtn;
}

/*
 *	[Apriori] Construct the L1 set.
 *
 *	@node_v is an inner-node.
 */
cns_list_t* const_L1_apriori( node_t* node_v, query_t* q, B_KEY_TYPE cost)
{
	int i;
	B_KEY_TYPE dist;
	cns_list_t* cns_list_v;
	cns_t* cns_v;

	cns_list_v = alloc_cns_list( );

	for( i=0; i<node_v->num; i++)
	{
		dist = calc_minDist_node( ( node_t*)( node_v->child[ i]), q->loc_v);
		if( dist >= cost || !is_relevant_node( ( node_t*)( node_v->child[ i]), q))
			continue;
		
		cns_v = alloc_cns( );
		add_cns_entry( cns_v, ( node_t*)( node_v->child[ i]));

		cns_v->next = cns_list_v->head->next;
		cns_list_v->head->next = cns_v;
		cns_list_v->cns_n ++;
	}

	return cns_list_v;
}

/*
 *	The checking for the join step.
 *
 *	Check whether two node sets @cns_v1 and @cns_v2
 *	share s-1 nodes, where s is the number of nodes in @cns_v1 or @cns_v2.
 *
 *	return the NULL if the join check fails,
 *	otherwise, return the (only) node_t in cns_v1 that is not contained in cns_v2.
 *
 */
node_t* join_check( cns_t* cns_v1, cns_t* cns_v2)
{
	int cnt;
	node_list_t* node_list_iter1, *node_list_iter2;
	node_t* node_v, *node_v1;

	node_v = NULL;		//Coding strategy.
	cnt = 0;
	node_list_iter1 = cns_v1->list_head->next;
	while( node_list_iter1 != NULL)
	{
		node_v1 = node_list_iter1->node_v;

		node_list_iter2 = cns_v2->list_head->next;
		while( node_list_iter2 != NULL)
		{
			if( node_list_iter2->node_v == node_v1)
				break;

			node_list_iter2 = node_list_iter2->next;
		}

		if( node_list_iter2 == NULL)
		{
			node_v = node_v1;
			cnt ++;

			if( cnt >= 2)
				return NULL;			
		}

		node_list_iter1 = node_list_iter1->next;
	}

	return node_v;
}

/*
 *	Check whether @cns_v contains the node @node_v.
 *
bool is_contained_node( cns_t* cns_v, node_t* node_v)
{
	node_list_t* node_list_iter;

	node_list_iter = cns_v->list_head->next;
	while( node_list_iter != NULL)
	{
		if( node_list_iter->node_v == node_v)
			return true;

		node_list_iter = node_list_iter->next;
	}

	return false;
}
*/

/*
 *	Merge two sets of nodes @cns_v1 and @cns_v2 into one.
 *
cns_t* merge_cns( cns_t* cns_v1, cns_t* cns_v2)
{
	cns_t* cns_v;
	node_list_t* node_list_v;

	cns_v = alloc_cns( );

	//Process cns_v1.
	node_list_v = cns_v1->list_head->next;
	while( node_list_v != NULL)
	{
		add_cns_entry( cns_v, node_list_v->node_v);

		node_list_v = node_list_v->next;
	}

	//Process cns_v2.
	node_list_v = cns_v2->list_head->next;
	while( node_list_v != NULL)
	{
		if( !is_contained_node( cns_v1, node_list_v->node_v))
		{
			add_cns_entry( cns_v, node_list_v->node_v);
			break;
		}
		
		node_list_v = node_list_v->next;
	}
	
	return cns_v;
}*/

/*
 *	Construct a cns_list_t structure based on a cns_list_set_t structure @cns_list_set_v.
 *
 *	Method: Store all cns_list_t structures in the cns_list_set_t structure @cns_list_set_v
 *			into a single cns_list_t structure.
 *
 *	@cns_list_set_v is destroyed.
 */
cns_list_t* reformat_cns_list_set( cns_list_set_t* cns_list_set_v)
{
	cns_list_t* cns_list_v, *cns_list_iter, *tmp;
	cns_t* cns_iter;

	cns_list_v = alloc_cns_list( );

	cns_list_iter = cns_list_set_v->head->next;
	while( cns_list_iter != NULL)
	{
		cns_iter = cns_list_iter->head;
		while( cns_iter->next != NULL)
			cns_iter = cns_iter->next;

		cns_iter->next = cns_list_v->head->next;
		cns_list_v->head->next = cns_list_iter->head->next;

		cns_list_v->cns_n += cns_list_iter->cns_n;

		tmp = cns_list_iter;
		cns_list_iter = cns_list_iter->next;
		
		free( tmp);

		/*s*/
		stat_v.memory_v -= sizeof( cns_list_t);
		/*s*/
	}

	cns_list_set_v->head->next = NULL;
	cns_list_set_v->list_n = 0;

	release_cns_list_set( cns_list_set_v);

	return cns_list_v;
}

/*
 *	Append the cns_t structure @cns_v2 to @cns_v1.
 */
void append_cns( cns_t* cns_v1, cns_t* cns_v2)
{
	node_list_t* node_list_iter;

	node_list_iter = cns_v2->list_head->next;
	while( node_list_iter != NULL)
	{
		add_cns_entry( cns_v1, node_list_iter->node_v);

		node_list_iter = node_list_iter->next;
	}
}

/*
 *	Check whether a cns_t structure @cns_v exists in a cns_list_t structure @cns_list_v.
 */
bool is_existing_cns( cns_list_t* cns_list_v, cns_t* cns_v)
{
	cns_t* iter;

	iter = cns_list_v->head->next;
	while( iter != NULL)
	{
		if( has_same_content_cns( iter, cns_v))
			return true;

		iter = iter->next;
	}

	return false;	
}

/*
 *	For a cns_list_t structure whose entries are formed by
 *	picking one cns_t structure from @cList[ i] for each i in [@sta, @end].
 */
cns_list_t* enumerate_node_set_sub( cns_list_t** cList, int sta, int end, query_t* q)
{
	cns_t *cns_v, *cns_iter1, *cns_iter2;
	cns_list_t *cns_list_v, *cns_list_v1;

	//Base step.
	if( sta == end)
		return copy_cns_list( cList[ sta]);

	//Recursive step.
	cns_list_v = alloc_cns_list( );

	cns_list_v1 = enumerate_node_set_sub( cList, sta+1, end, q);

	cns_iter1 = cList[ sta]->head->next;
	while( cns_iter1 != NULL)
	{
		cns_iter2 = cns_list_v1->head->next;
		while( cns_iter2 != NULL)
		{
			cns_v = alloc_cns( );

			//Append cns_iter2 to cns_v.
			append_cns( cns_v, cns_iter2);

			//Append cns_iter1 to cns_v.
			append_cns( cns_v, cns_iter1);

			//Add cns_v into cns_list_v;
			//if( is_covered_cns( cns_v, q))
			add_cns_list_entry( cns_list_v, cns_v);

			cns_iter2 = cns_iter2->next;
		}

		cns_iter1 = cns_iter1->next;
	}
	
	return cns_list_v;
}

/*
 *	The implementation of the "EnumeratNodeSets" algorithm (Algorithm 7 of Cao's paper)
 */
cns_list_t* enumerate_node_set( cns_t* cns_v, B_KEY_TYPE cost, query_t* q)
{
	int i, s, list_n;
	node_t* node_v, *node_v1;
	node_list_t* n_list_v;
	cns_t* cns_v1, *cns_v2, *cns_v3;
	cns_list_t* cns_list_v, *cns_list_cur, *rtn;
	cns_list_t** cList;
	cns_list_set_t* cns_list_set_v;

	//Initialize the structure for storing cList_1, cList_2, ....
	list_n = cns_v->node_n;
	cList = ( cns_list_t**)malloc( list_n * sizeof( cns_list_t*));
	memset( cList, 0, list_n * sizeof( cns_list_t*));


	/*s*/
	stat_v.memory_v += list_n * sizeof( cns_list_t*);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	n_list_v = cns_v->list_head->next;
	for( i=0; i<list_n; i++)
	{
		node_v = n_list_v->node_v;

		//Initailize the structure @cns_list_set_v for storing L_i^1, L_i^2, ....
		cns_list_set_v = alloc_cns_list_set( );
		cns_list_cur = cns_list_set_v->head;
		
		//Apriori process.
		//Construct L_i^1.
		cns_list_v = const_L1_apriori( node_v, q, cost);
		
		cns_list_cur->next = cns_list_v;
		cns_list_set_v->list_n ++;

		cns_list_cur = cns_list_v;

		//Construct L_i^2, L_i^3, ....
		for( s=2; s <= q->psi_v->key_n - list_n + 1; s++)
		{
			cns_list_v = alloc_cns_list( );

			//cns_list_cur: the base set.
			cns_v1 = cns_list_cur->head->next;
			if( cns_v1 == NULL)
			{
				//cns_v1 could be NULL.
				break;				
			}

			while( cns_v1->next != NULL)
			{
				cns_v2 = cns_v1->next;
				while( cns_v2 != NULL)
				{
					//If cns_v1 and cns_v2 share the first s-1 nodes.
					if( ( node_v1 = join_check( cns_v1, cns_v2)))
					{
						//cns_v3 = merge_cns( cns_v1, cns_v2);
						cns_v3 = copy_cns( cns_v2);
						add_cns_entry( cns_v3, node_v1);

						if( calc_low_cost( cns_v3, q) < cost && 
							!is_existing_cns( cns_list_v, cns_v3))
							add_cns_list_entry( cns_list_v, cns_v3);
						else
							release_cns( cns_v3);
					}

					cns_v2 = cns_v2->next;
				}

				cns_v1 = cns_v1->next;
			}

			/*t*/
			//print_cns_list( cns_list_v);
			/*t*/

			cns_list_cur->next = cns_list_v;
			cns_list_set_v->list_n ++;

			cns_list_cur = cns_list_v;

			if( cns_list_cur->cns_n <= 1)
				break;
		}//for(s)

		//Maintain cList_i, i.e., cList[ i].
		cList[ i] = reformat_cns_list_set( cns_list_set_v);

		/*t*/
		//print_cns_list( cList[ i]);
		/*t*/

		//Release cns_list_set_v.
		//release_cns_list_set( cns_list_set_v);

		n_list_v = n_list_v->next;
	}

	//Construct setList;
	rtn = enumerate_node_set_sub( cList, 0, list_n-1, q);	


	/*t*/
	//print_cns_list( rtn);
	/*t*/

	//Release the resources.
	for( i=0; i<list_n; i++)
		release_cns_list( cList[ i]);
	free( cList);
	
	/*s*/
	stat_v.memory_v -= list_n * sizeof( cns_list_t*);
	/*s*/

	return rtn;
}

/*
 *	Check whether the keywords of a query @q are covered by a node set @cns_v.
 */
bool is_covered_cns( cns_t* cns_v, query_t* q)
{
	KEY_TYPE key;
	k_node_t* k_node_iter;
	node_list_t* node_list_iter;

	k_node_iter = q->psi_v->k_head->next;
	while( k_node_iter != NULL)
	{
		key = k_node_iter->key;

		node_list_iter = cns_v->list_head->next;
		while( node_list_iter != NULL)
		{
			if( has_key_node( node_list_iter->node_v, key))
				break;

			node_list_iter = node_list_iter->next;
		}

		if( node_list_iter == NULL)
			return false;

		k_node_iter = k_node_iter->next;
	}

	return true;
}

/*
 *	The implementation of the "Cao-Exact" algorithm (Algorithm 6 of Cao's Paper).
 *
 *	Heap option 1: b_heap_t structure.
 *	Heap option 2: bst_t structure.
 *	
 *	The option 2 is impossible since the dead lock: 
 *		bst_node_t-->cns_t-->node_t-->bst_node_t.
 *
 *	Since Option 1 requires an input size, which is difficult to estimate,
 *	we adopt option 2.
 *	We use the "void*" type in bst_node_t in order to break the above dead lock.
 */
obj_set_t* Cao_Exact( query_t* q)
{
	//printf("Cao_Exact\n");
	int cnt;
	bst_t* heap_v;
	cns_t* cns_v, *cns_list_iter, *cns_tmp; 
	obj_set_t* V, *V_1, *V_tmp;
	B_KEY_TYPE costV, costV_1;
	bst_node_t* bst_node_v, *top;
	cns_list_t* cns_list_v;
	
	//printf( "Cao-Exact:\n");

	//Compute an approximate solution.
	V = Cao_Appro1( q);
	if( V == NULL)
	{
		return NULL;
	}

	costV = comp_cost( V, q);
	
	//Initialize a queue.
	heap_v = bst_ini( );

	//Enqueue the root.
	cns_v = alloc_cns( );
	add_cns_entry( cns_v, IRTree_v.root);

	bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
	memset( bst_node_v, 0, sizeof( bst_node_t));

	bst_node_v->key = calc_low_cost( cns_v, q);
	bst_node_v->cns_v = ( void*)cns_v;
	bst_insert( heap_v, bst_node_v);

	/*s*/
	stat_v.memory_v += sizeof( bst_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	top = bst_get_min( heap_v->root);

	cnt = 0;
	//Best-first search process.
	while( top != NULL)
	{
		cns_v = ( cns_t*)top->cns_v;

		/*s*/
		cnt ++;
		if( cnt % 100 == 0)
			printf( "#cns: %i\n", cnt);
		/*t*/

		/*t*/
		//print_cns( cns_v);
		/*t*/

		if( calc_low_cost( cns_v, q) >= costV)
			break;

		if( cns_v->list_head->next->node_v->level == 0)
		{
			//The nodes in cns_v are leaf-nodes.
			V_1 = exhaustive_node_set_search( cns_v, costV, q);
			if( V_1)
			{
				costV_1 = comp_cost( V_1, q);
				if( costV_1 < costV)	// && is_covered_cns( cns_v, q)
				{
					V_tmp = V;
					V = V_1;
					release_obj_set( V_tmp);
					
					costV = costV_1;
				}
				else
					release_obj_set( V_1);
			}
		}
		else
		{
			//The nodes in cns_v are inner-nodes.
			cns_list_v = enumerate_node_set( cns_v, costV, q);

			cns_tmp = cns_list_v->head;
			cns_list_iter = cns_tmp->next;
			while( cns_list_iter != NULL)
			{
				if( is_covered_cns( cns_list_iter, q))
				{
					//Un-link the corresponding cns_t structure from the list.
					cns_tmp->next = cns_list_iter->next;
					cns_list_iter->next = NULL;

					//
					bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
					memset( bst_node_v, 0, sizeof( bst_node_t));

					bst_node_v->cns_v = cns_list_iter;
					bst_node_v->key = calc_low_cost( cns_list_iter, q);
					bst_insert( heap_v, bst_node_v);

					/*s*/
					stat_v.memory_v += sizeof( bst_node_t);
					if( stat_v.memory_v > stat_v.memory_max)
						stat_v.memory_max = stat_v.memory_v;
					/*s*/					

					cns_list_iter = cns_tmp->next;
					continue;
				}
				
				cns_tmp = cns_list_iter;
				cns_list_iter = cns_tmp->next;
			}//while

			release_cns_list( cns_list_v);
		}

		release_cns( cns_v);
		bst_delete( heap_v, top);
		free( top);			//Memory usage consideration.

		/*s*/
		stat_v.node_set_n ++;
		stat_v.memory_v -= sizeof( bst_node_t);
		/*s*/

		top = bst_get_min( heap_v->root);
	}

	bst_release( heap_v);

	return V;
}
