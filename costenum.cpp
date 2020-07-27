/*
 *	Author: Victor Junqiu WEI
 *	Email: wjqjsnj@gmail.com
 */


#include "Voronoic.h"
#include "cao_alg.h"
#include "toplevel_tree.h"



double distance_sqr(struct point_t const* a, struct point_t const* b)
{
    return (a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y);
}

double distance(struct point_t const* a, struct point_t const* b)
{
    return sqrt(distance_sqr(a, b));
}

int intersect(struct circle_t const circles[], struct point_t intersections[])
{
    double a, b, c, p, q, s; 
    double cos_value[2], sin_value[2]; 
    double d 
             = distance(&circles[0].center, &circles[1].center);

    if (d > circles[0].r + circles[1].r
        || d < fabs(circles[0].r - circles[1].r))
    {
        return 0;
    }

    a = 2.0 * circles[0].r * (circles[0].center.x - circles[1].center.x);
    b = 2.0 * circles[0].r * (circles[0].center.y - circles[1].center.y);
    c = circles[1].r * circles[1].r - circles[0].r * circles[0].r
        - distance_sqr(&circles[0].center, &circles[1].center);
    p = a * a + b * b;
    q = -2.0 * a * c;

   
    if (d == circles[0].r + circles[1].r
     || d == fabs(circles[0].r - circles[1].r))
    {
        cos_value[0] = -q / p / 2.0;
        sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);

        intersections[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
        intersections[0].y = circles[0].r * sin_value[0] + circles[0].center.y;

        
        if (distance_sqr(&intersections[0], &circles[1].center)
            != circles[1].r * circles[1].r)
        {
            intersections[0].y = circles[0].center.y
                                 - circles[0].r * sin_value[0];
        }
        return 1;
    }

    s = c * c - b * b;
    cos_value[0] = (sqrt(q * q - 4.0 * p * s) - q) / p / 2.0;
    cos_value[1] = (-sqrt(q * q - 4.0 * p * s) - q) / p / 2.0;
    sin_value[0] = sqrt(1 - cos_value[0] * cos_value[0]);
    sin_value[1] = sqrt(1 - cos_value[1] * cos_value[1]);

    intersections[0].x = circles[0].r * cos_value[0] + circles[0].center.x;
    intersections[1].x = circles[0].r * cos_value[1] + circles[0].center.x;
    intersections[0].y = circles[0].r * sin_value[0] + circles[0].center.y;
    intersections[1].y = circles[0].r * sin_value[1] + circles[0].center.y;

    
    if (distance_sqr(&intersections[0], &circles[1].center)
        != circles[1].r * circles[1].r)
    {
        intersections[0].y = circles[0].center.y - circles[0].r * sin_value[0];
    }
    if (distance_sqr(&intersections[1], &circles[1].center)
        != circles[1].r * circles[1].r)
    {
        intersections[1].y = circles[0].center.y - circles[0].r * sin_value[1];
    }

   
    if (intersections[0].y == intersections[1].y
     && intersections[0].x == intersections[1].x)
    {
        intersections[1].y = -intersections[1].y;
    }
    return 2;
}
/*
 *	Allocate a loc_t structure.
 */
loc_t* alloc_loc( int dim)
{
	loc_t* loc_v;

	loc_v = ( loc_t*)malloc( sizeof( loc_t));
	memset( loc_v, 0, sizeof( loc_t));

	loc_v->coord = ( float*)malloc( dim * sizeof( float));
	memset( loc_v->coord, 0, dim * sizeof( float));

	loc_v->dim = dim;

	return loc_v;
}

/*
 *	Retrieve the loc_t information of an obj_t structure @obj_v.
 */
loc_t* get_obj_loc( obj_t* obj_v)
{
	int i;
	loc_t* loc_v;

	loc_v = alloc_loc( IRTree_v.dim);
	for( i=0; i<IRTree_v.dim; i++)
		loc_v->coord[ i] = obj_v->MBR[ i].min;
	
	return loc_v;
}

/*
 *	Copy a loc_t structure.
 */
loc_t* copy_loc( loc_t* loc_v)
{
	int j;
	loc_t* loc_v1;

	loc_v1 = alloc_loc( loc_v->dim);

	for( j=0; j<loc_v->dim; j++)
		loc_v1->coord[ j] = loc_v->coord[ j];

	return loc_v1;
}

/*
 *	Release a loc_t structure.
 */
void release_loc( loc_t* loc_v)
{
	free( loc_v->coord);
	free( loc_v);
}

/*
 *	Allocate a psi_t structure.
 */
psi_t* alloc_psi( )
{
	psi_t* psi_v;

	psi_v = ( psi_t*)malloc( sizeof( psi_t));
	memset( psi_v, 0, sizeof( psi_t));

	psi_v->k_head = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( psi_v->k_head, 0, sizeof( k_node_t));



	return psi_v;
}

/*
 *	Add an keyword into the psi_t structure.
 */
void add_psi_entry( psi_t* psi_v, KEY_TYPE key)
{
	k_node_t* k_node_v;

	k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( k_node_v, 0, sizeof( k_node_t));

	k_node_v->key = key;
	k_node_v->next =  psi_v->k_head->next;
	psi_v->k_head->next = k_node_v;
	psi_v->key_n ++;
}

/*
 *	Release a psi_t structure.
 */
void release_psi( psi_t* psi_v)
{
	k_node_t* k_node_v1, *k_node_v2;

	k_node_v1 = psi_v->k_head;
	while( k_node_v1->next != NULL)
	{
		k_node_v2 = k_node_v1->next;
		free( k_node_v1);
		k_node_v1 = k_node_v2;
	}
	free( k_node_v1);
}

/*
 *	Allocate a query_t structure.
 */
query_t* alloc_query( )
{
	query_t* q;

	q = ( query_t*)malloc( sizeof( query_t));
	memset( q, 0, sizeof( query_t));

	//q->loc_v = alloc_loc( dim);

	//q->psi_v = alloc_psi( );

	return q;
}

/*
 *	Print a query_t structure.
 */
void print_query( query_t* q, FILE* o_fp)
{
	int i;
	k_node_t* k_node_v;

	//Location.
	fprintf( o_fp, "%f", q->loc_v->coord[ 0]);
	for( i=1; i<q->loc_v->dim; i++)
	{
		fprintf( o_fp, ",%f", q->loc_v->coord[ i]);
	}
	fprintf( o_fp, "\n");

	//Keywords.
	fprintf( o_fp, "%i", q->psi_v->key_n);
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		fprintf( o_fp, ",%i", ( int)( k_node_v->key));
		k_node_v = k_node_v->next;
	}
	fprintf( o_fp, "\n\n");

	return;
}

/*
 *	Read a query instance from a file stream @i_fp.
 */
query_t* read_query( FILE* i_fp)
{
	int i;
	char des;
	KEY_TYPE key;
	query_t* q;
	k_node_t* k_node_v;

	q = alloc_query( );

	q->loc_v = alloc_loc( IRTree_v.dim);
	q->psi_v = alloc_psi( );

	//Location.
	if( fscanf( i_fp, "%f", &q->loc_v->coord[ 0]) == EOF)
	{
		release_query( q);
		return NULL;
	}

	for( i=1; i<IRTree_v.dim; i++)
		fscanf( i_fp, "%c%f", &des, &q->loc_v->coord[ i]);

	//Keywords.
	k_node_v = q->psi_v->k_head;
	fscanf( i_fp, "%i", &q->psi_v->key_n);
	for( i=0; i<q->psi_v->key_n; i++)
	{
		fscanf( i_fp, "%c%lf", &des, &key);
		add_keyword_entry( k_node_v, key);
	}

	return q;		
}

/*
 *	Release a query_t structure.
 */
void release_query( query_t* q)
{
	release_loc( q->loc_v);
	release_psi( q->psi_v);

	free( q);
}

/*
 *	Allocate the memory for a disk_t structure.
 */
disk_t* alloc_disk( int dim)
{
	disk_t* disk_v;

	disk_v = ( disk_t*)malloc( sizeof( disk_t));
	memset( disk_v, 0, sizeof( disk_t));

	disk_v->center = alloc_loc( dim);

	/*s*/
	stat_v.memory_v += sizeof( disk_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
	

	return disk_v;
}

/*
 *	Initialize a disk_t structure.
 */
void set_disk( disk_t* disk_v, loc_t* loc_v, B_KEY_TYPE radius)
{
	int i;

	disk_v->radius = radius;
	
	for( i=0; i<loc_v->dim; i++)
		disk_v->center->coord[ i] = loc_v->coord[ i];
}

/*
 *	Construct a disk_t structure with its center of @loc_v and
 *	its radius of @radius.
 */
disk_t* const_disk( loc_t* loc_v, B_KEY_TYPE radius)
{
	int i;
	disk_t* disk_v;


	disk_v = alloc_disk( IRTree_v.dim);
	for( i=0; i<IRTree_v.dim; i++)
		disk_v->center->coord[ i] = loc_v->coord[ i];
	disk_v->center->dim = IRTree_v.dim;

	disk_v->radius = radius;

	return disk_v;
}

/*
 *	Release a disk_t structure.
 */
void release_disk( disk_t* disk_v)
{
	release_loc( disk_v->center);
	free( disk_v);

	/*s*/
	stat_v.memory_v -= sizeof( disk_t);
	/*s*/
}

/*
 *	Allocate the memory for an obj_set_t structure.
 */
obj_set_t* alloc_obj_set( )
{
	obj_set_t* obj_set_v;

	obj_set_v = ( obj_set_t*)malloc( sizeof( obj_set_t));
	memset( obj_set_v, 0, sizeof( obj_set_t));

	obj_set_v->head = ( obj_node_t*)malloc( sizeof( obj_node_t));
	memset( obj_set_v->head, 0, sizeof( obj_node_t));

	/*s*/
	stat_v.memory_v += sizeof( obj_set_t) + sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return obj_set_v;
}

/*
 *	Copy an obj_set_t structure.
 */
obj_set_t* copy_obj_set( obj_set_t* obj_set_v)
{
	obj_set_t* rtn;
	obj_node_t* obj_node_v, *obj_node_iter;

	rtn= alloc_obj_set( );
	obj_node_v = rtn->head;

	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		obj_node_v->next = ( obj_node_t*)malloc( sizeof( obj_node_t));
		memset( obj_node_v->next, 0, sizeof( obj_node_t));

		obj_node_v = obj_node_v->next;
		obj_node_v->obj_v = obj_node_iter->obj_v;
		
		obj_node_iter = obj_node_iter->next;
	}

	rtn->obj_n = obj_set_v->obj_n;

	/*s*/
	stat_v.memory_v += rtn->obj_n * sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return rtn;		
}

/*
 *	Check whether two obj_set_t structures obj_set_v1 and obj_set_v2 have the same content.
 */
bool has_same_content_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	obj_node_t* iter1, *iter2;

	if( obj_set_v1->obj_n != obj_set_v2->obj_n)
		return false;

	iter1 = obj_set_v1->head->next;
	while( iter1 != NULL)
	{
		iter2 = obj_set_v2->head->next;
		while( iter2 != NULL)
		{
			if( iter2->obj_v == iter1->obj_v)
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
 *	Remove the identical objs from @obj_set_v.
 */
void remove_identical_obj( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_iter1, *obj_node_iter2, *tmp;

	obj_node_iter1 = obj_set_v->head->next;
	while( obj_node_iter1->next != NULL)
	{
		tmp = obj_node_iter1;	
		obj_node_iter2 = obj_node_iter1->next;
		while( obj_node_iter2 != NULL)
		{
			if( obj_node_iter2->obj_v == obj_node_iter1->obj_v)
			{
				//remove.
				tmp->next = obj_node_iter2->next;
				obj_set_v->obj_n --;
				
				free( obj_node_iter2);
				
				/*s*/
				stat_v.memory_v -= sizeof( obj_node_t);
				/*s*/

				obj_node_iter2 = tmp->next;
				continue;
			}

			tmp = obj_node_iter2;
			obj_node_iter2 = obj_node_iter2->next;
		}

		obj_node_iter1 = obj_node_iter1->next;
		if( obj_node_iter1 == NULL)
			break;
	}

	return;
}

/*
 *	Print an obj_set_t structure.
 */
void print_obj_set( obj_set_t* obj_set_v, FILE* o_fp)
{
	int i;
	obj_node_t* obj_node_iter;
	obj_t* obj_v;

	if( !obj_set_v)
		return;

	for( i=0; i<20; i++)
		fprintf( o_fp, "==");
	fprintf( o_fp, "\n");

	fprintf( o_fp, "%i\n\n", obj_set_v->obj_n);
	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		obj_v = obj_node_iter->obj_v;
		fprintf( o_fp, "%i:\t", obj_v->id);
		print_k_list( obj_v->k_head, o_fp);

		obj_node_iter = obj_node_iter->next;
	}
	fprintf( o_fp, "\n");
}

/*
 *	Release the memory of an obj_set_t structure.
 */
void release_obj_set( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v1, *obj_node_v2;

	if( !obj_set_v)
		return;

	/*s*/
	stat_v.memory_v -= sizeof( obj_set_t) + 
						( obj_set_v->obj_n + 1) * sizeof( obj_node_t);
	/*s*/

	obj_node_v1 = obj_set_v->head;
	while( obj_node_v1->next != NULL)
	{
		obj_node_v2 = obj_node_v1->next;
		free( obj_node_v1);

		obj_node_v1 = obj_node_v2;
	}
	free( obj_node_v1);

	free( obj_set_v);
}

/*
 *	Allocate a IF_entry_t structure.
 *
IF_entry_t* alloc_IF_entry( )
{
	IF_entry_t* IF_entry_v;

	IF_entry_v = ( IF_entry_t*)malloc( sizeof( IF_entry_t));
	memset( IF_entry_v, 0, sizeof( IF_entry_t));

	IF_entry_v->p_list = alloc_obj_set( );

	return IF_entry_v;
}*/

/*
 *	Release a IF_entry_t structure.
 *
void release_IF_entry( IF_entry_t* IF_entry_v)
{
	release_obj_set( IF_entry_v->p_list);
	free( IF_entry_v);
}*/

/*
 *	Allocate a IF_t structure.
 *
IF_t* alloc_IF( )
{
	int i;
	IF_t* IF_v;

	IF_v = ( IF_t*)malloc( sizeof( IF_t));
	memset( IF_v, 0, sizeof( IF_t));

	IF_v->bst_v = bst_ini( );

	//IF_v->key_n = key_n;

	//IF_v->entry_v = ( IF_entry_t*)malloc( key_n * sizeof( IF_entry_t));
	//memset( IF_v->entry_v, 0, key_n * sizeof( IF_entry_t));

	//for( i=0; i<key_n; i++)
		//IF_v->entry_v[ i].p_list = alloc_obj_set( );
		//IF_v->entry_v[ i] = alloc_IF_entry( );

	return IF_v;
}*/

/*
 *	Release a IF_t structure.
 *
void release_IF( IF_t* IF_v)
{
	//int i;
	
	//for( i=0; i<IF_v->key_n; i++)
		//release_obj_set( IF_v->entry_v[ i].p_list);
		//release_IF_entry( IF_v->entry_v[ i]);

	//free( IF_v->entry_v);
	bst_release( IF_v->bst_v);
	free( IF_v);
}*/

/*
 *	Calculate the distance between two locations @loc_v1 and @loc_v2.
 */
B_KEY_TYPE calc_dist_loc( loc_t* loc_v1, loc_t* loc_v2)
{
	int i;
	B_KEY_TYPE dist;

	dist = 0;
	for( i=0; i<loc_v1->dim; i++)
		dist += pow( loc_v1->coord[ i] - loc_v2->coord[ i], 2);

	return sqrt( dist);
}

/*
 *	Calculate the distance between two objects @obj_v1 and @obj_v2.
 */
B_KEY_TYPE calc_dist_obj( obj_t* obj_v1, obj_t* obj_v2)
{
	B_KEY_TYPE dist;
	loc_t* loc_v1, *loc_v2;

	loc_v1 = get_obj_loc( obj_v1);
	loc_v2 = get_obj_loc( obj_v2);

	dist = calc_dist_loc( loc_v1, loc_v2);

	release_loc( loc_v1);
	release_loc( loc_v2);

	return dist;	
}

B_KEY_TYPE calc_dist_int(query_ssq * q, obj* obj1,obj* obj2){
	B_KEY_TYPE dist,radius,q1,q2;
	loc_t* loc1,*loc2,*loci1,*loci2;
	loc1=get_obj_loc(obj1);
	loc2=get_obj_loc(obj2);
	radius=calc_dist_loc(loc1,loc2);
	q1=calc_dist_loc(q->loc_v,loc1);
	q2=calc_dist_loc(q->loc_v,loc2);
	if(q1<=radius && q2<=radius) return 0.0;
	loci1->coord[0]=loc1->coord[0]-(radius)/(q1)*(loc1->coord[0]-q->loc_v->coord[0]);
	loci1->coord[1]=loc1->coord[1]-(radius)/(q1)*(loc1->coord[1]-q->loc_v->coord[1]);
	loci2->coord[0]=loc2->coord[0]-(radius)/(q2)*(loc2->coord[0]-q->loc_v->coord[0]);
	loci2->coord[1]=loc2->coord[1]-(radius)/(q2)*(loc2->coord[1]-q->loc_v->coord[1]);
	loci1->dim=2;
	loci2->dim=2;

	if(calc_dist_loc(loci1,loc2)<=radius && calc_dist_loc(loci2,loc1)<=radius) return min(q1-radius,q2-radius);
	if(calc_dist_loc(loci1,loc2)>radius && calc_dist_loc(loci2,loc1)<=radius) return q2-radius;
	if(calc_dist_loc(loci1,loc2)<=radius && calc_dist_loc(loci2,loc1)>radius) return q1-radius;

	loc_t* loc3,*loc4;
	struct circle_t circles[2];
    struct point_t points[2];
	circles[0].center.x=loc1->coord[0];
	circles[0].center.y=loc1->coord[1];
	circles[1].center.x=loc2->coord[0];
	circles[1].center.y=loc2->coord[1];
	circles[0].r=radius;
	circles[1].r=radius;

	    switch (intersect(circles, points)) {
        case 0:
            puts("No intersection.");
            break;
        case 1:
            printf("error");
            break;
    }
		loc3->dim=2;
		loc4->dim=2;
		loc3->coord[0]=points[0].x;
		loc3->coord[1]=points[0].y;
		loc4->coord[0]=points[1].x;
		loc4->coord[1]=points[1].y;
		return min(calc_dist_loc(q->loc_v,loc3),calc_dist_loc(q->loc_v,loc4));

}
/*
 *	Check whether an object @obj_v contains the keyword @key.
 */
bool has_key_obj( obj_t* obj_v, KEY_TYPE key)
{
	k_node_t* k_node_v;

	k_node_v = obj_v->k_head->next;
	while( k_node_v != NULL)
	{
		if( k_node_v->key == key)
			return true;

		k_node_v = k_node_v->next;
	}

	return false;
}

/*
 *	Check whether an object @obj_v is "relevant" to the query @q.
 *	That is, whether @obj_v contains a keyword in the query @q.
 */
bool is_relevant_obj( obj_t* obj_v, query_t* q)
{
	KEY_TYPE key;
	k_node_t* k_node_v;

	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		key = k_node_v->key;
		if( has_key_obj( obj_v, key))
			return true;

		k_node_v = k_node_v->next;
	}

	return false;		
}

/*
 *	Check whether the keywords in the query @q are covered by a set of objs in @obj_set_v.
 */
bool is_covered_obj_set( obj_set_t* obj_set_v, query_t* q)
{
	KEY_TYPE key;
	k_node_t* k_node_iter;
	obj_node_t* obj_node_iter;

	k_node_iter = q->psi_v->k_head->next;
	while( k_node_iter != NULL)
	{
		key = k_node_iter->key;

		obj_node_iter = obj_set_v->head->next;
		while( obj_node_iter != NULL)
		{
			if( has_key_obj( obj_node_iter->obj_v, key))
				break;

			obj_node_iter = obj_node_iter->next;
		}

		if( obj_node_iter == NULL)
			return false;

		k_node_iter = k_node_iter->next;
	}
	
	return true;
}

/*
 *	Check whether the sub-tree rooted at a node @node_v
 *	contains a keyword @key.
 *
 */
BIT_TYPE has_key_node( node_t* node_v, KEY_TYPE key)
{
	bst_node_t* bst_node_v;

	if( ( bst_node_v = bst_search( node_v->bst_v, key)))
		return bst_node_v->p_list;
	
	return 0;
}

/*
 *	Check whether a node @node_v is "relevant" or not.
 */
BIT_TYPE is_relevant_node( node_t* node_v, query_t* q)
{
	BIT_TYPE res, res_t;

	k_node_t* k_node_v;

	res = 0;
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		res_t = has_key_node( node_v, k_node_v->key);
		union_bit( res, res_t);

		if( res == UINT_MAX)
			return res;

		k_node_v = k_node_v->next;
	}

	return res;
}

/*
 *	Calculate the minimum distance between the MBR of a node @node_v 
 *	and a location @loc_v.
 */
B_KEY_TYPE calc_minDist( range* MBR, loc_t* loc_v)
{
	int i;
	B_KEY_TYPE m_dist;
	
	//Calculate the minimum distance between the MBR and the location.
	m_dist = 0;
	for( i=0; i<loc_v->dim; i++)
	{
		if( loc_v->coord[ i] < MBR[ i].min)
			m_dist += pow( MBR[ i].min - loc_v->coord[ i], 2);
		else if( loc_v->coord[ i] > MBR[ i].max)
			m_dist += pow( loc_v->coord[ i] - MBR[ i].max, 2);
	}

	return sqrt( m_dist);
}

/*
 *	Calculate the maxDist between a MBR @MBR and 
 */
B_KEY_TYPE calc_maxDist( range* MBR, loc_t* loc_v)
{
	int i;
	B_KEY_TYPE maxDist;

	maxDist = 0;
	for( i=0; i<loc_v->dim; i++)
	{
		if( loc_v->coord[ i] <= ( MBR[ i].min + MBR[ i].max) / 2)
			maxDist += pow( MBR[ i].max - loc_v->coord[ i], 2);
		else
			maxDist += pow( loc_v->coord[ i] - MBR[ i].min, 2);
	}

	return sqrt( maxDist);
}

/*
 *	Calculate the minimum distance between the MBR of a node @node_v and a location @loc_v.
 */
B_KEY_TYPE calc_minDist_node( node_t* node_v, loc_t* loc_v)
{
	B_KEY_TYPE dist;
	range* MBR;

	if( node_v->parent == NULL)
		MBR = get_MBR_node( node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[ node_v->loc];

	dist = calc_minDist( MBR, loc_v);

	if( node_v->parent == NULL)
	{
		free( MBR);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}

	return dist;
}

B_KEY_TYPE calc_minDist_nodep( node_t* node_v1, node_t* node_v2){
	B_KEY_TYPE dist[8],dist_v;
	range* MBR1;
	range* MBR2;
	loc_t * loc_v;
	//dist=100000000000000.0;
	if( node_v1->parent == NULL)
		MBR1 = get_MBR_node( node_v1, IRTree_v.dim);
	else
		MBR1 = node_v1->parent->MBRs[ node_v1->loc];

	if( node_v2->parent == NULL)
		MBR2 = get_MBR_node( node_v2, IRTree_v.dim);
	else
		MBR2 = node_v2->parent->MBRs[ node_v2->loc];
	loc_v=alloc_loc(2);
	loc_v->coord[0]=MBR1[0].min;
	loc_v->coord[1]=MBR1[1].min;
	dist[0]=calc_minDist(MBR2, loc_v);
	
	loc_v->coord[0]=MBR1[0].max;
	loc_v->coord[1]=MBR1[1].min;
	dist[1]=calc_minDist(MBR2, loc_v);

	loc_v->coord[0]=MBR1[0].max;
	loc_v->coord[1]=MBR1[1].max;
	dist[2]=calc_minDist(MBR2, loc_v);

	loc_v->coord[0]=MBR1[0].min;
	loc_v->coord[1]=MBR1[1].max;
	dist[3]=calc_minDist(MBR2, loc_v);

	loc_v->coord[0]=MBR2[0].min;
	loc_v->coord[1]=MBR2[1].min;
	dist[4]=calc_minDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].max;
	loc_v->coord[1]=MBR2[1].min;
	dist[5]=calc_minDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].max;
	loc_v->coord[1]=MBR2[1].max;
	dist[6]=calc_minDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].min;
	loc_v->coord[1]=MBR2[1].max;
	dist[7]=calc_minDist(MBR1, loc_v);

	dist_v=dist[7];
	for(int i=0;i<7;i++){
		dist_v=min(dist[i],dist_v);
	}
	if( node_v1->parent == NULL)
	{
		free( MBR1);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}

	if( node_v2->parent == NULL)
	{
		free( MBR2);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}
	return dist_v;
}

B_KEY_TYPE calc_maxDist_nodep( node_t* node_v1, node_t* node_v2){
	B_KEY_TYPE dist[8],dist_v;
	range* MBR1;
	range* MBR2;
	loc_t * loc_v=alloc_loc(2);
	//dist=100000000000000.0;
	if( node_v1->parent == NULL)
		MBR1 = get_MBR_node( node_v1, IRTree_v.dim);
	else
		MBR1 = node_v1->parent->MBRs[ node_v1->loc];

	if( node_v2->parent == NULL)
		MBR2 = get_MBR_node( node_v2, IRTree_v.dim);
	else
		MBR2 = node_v2->parent->MBRs[ node_v2->loc];
//	loc_v->dim=2;
	loc_v->coord[0]=MBR1[0].min;
	loc_v->coord[1]=MBR1[1].min;
	dist[0]=calc_maxDist(MBR2, loc_v);
	
	loc_v->coord[0]=MBR1[0].max;
	loc_v->coord[1]=MBR1[1].min;
	dist[1]=calc_maxDist(MBR2, loc_v);

	loc_v->coord[0]=MBR1[0].max;
	loc_v->coord[1]=MBR1[1].max;
	dist[2]=calc_maxDist(MBR2, loc_v);

	loc_v->coord[0]=MBR1[0].min;
	loc_v->coord[1]=MBR1[1].max;
	dist[3]=calc_maxDist(MBR2, loc_v);

	loc_v->coord[0]=MBR2[0].min;
	loc_v->coord[1]=MBR2[1].min;
	dist[4]=calc_maxDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].max;
	loc_v->coord[1]=MBR2[1].min;
	dist[5]=calc_maxDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].max;
	loc_v->coord[1]=MBR2[1].max;
	dist[6]=calc_maxDist(MBR1, loc_v);

	loc_v->coord[0]=MBR2[0].min;
	loc_v->coord[1]=MBR2[1].max;
	dist[7]=calc_maxDist(MBR1, loc_v);

	dist_v=dist[7];
	for(int i=0;i<7;i++){
		dist_v=max(dist[i],dist_v);
	}
	if( node_v1->parent == NULL)
	{
		free( MBR1);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}

	if( node_v2->parent == NULL)
	{
		free( MBR2);

		/*s*/
		stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
		/*s*/
	}
	release_loc(loc_v);
	return dist_v;
	
}

B_KEY_TYPE calc_maxDist_nodep( node_t* node_v1, node_t* node_v2);
/*
 *	Check whether a MBR @MBR overlaps with a disk @disk_v.
 */
bool is_overlap( range* MBR, disk_t* disk_v)
{
	B_KEY_TYPE min_dist;

	min_dist = calc_minDist( MBR, disk_v->center);
	if(  min_dist <= disk_v->radius)
		return true;
	else
		return false;
}

/*
 *	Constrained NN search.
 *	Find the NN that contains the keyword @key and
 *	is located in the dist @disk_v.
 *
 *	disk_v == NULL indicates that the disk is the whole space,
 *	thus, in this case, this constraint is invalid (interface consideration).
 *
 *	Method: best-first search.
 *	
 */
obj_t* const_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v)
{	
	int size, top, i, rear;
	BIT_TYPE p_list;

	B_KEY_TYPE min_dist, c_dist;
	range* MBR;
	obj_t* obj_v;
	node_t* node_v;
	b_heap_t* b_h;

	//Checking 1.
	if( IRTree_v.root->num == 0)
		return NULL;

	//Checking 2: involved in the search process implicitly.
	//if( !( bst_node_v = bst_search( IRTree.root->bst_v, key)))
		//return NULL;

	//Checking 3:
	MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
	if( disk_v != NULL && !is_overlap( MBR, disk_v))
		return NULL;

	size = IRTree_v.obj_n / M + 1;
	b_h = alloc_b_heap( size);

	rear = 1;
	b_h->obj_arr[ rear].node_v = IRTree_v.root;
	b_h->obj_arr[ rear].key = calc_minDist( MBR, loc_v);
	free( MBR);

	/*s*/
	stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
	/*s*/

	b_h_insert( b_h, rear++);

	min_dist = DBL_MAX;
	obj_v = NULL;

	while( !b_h_is_empty( b_h))
	{
		top = b_h_get_top( b_h);
		if( b_h->obj_arr[ top].key >= min_dist)
			break;

		node_v = b_h->obj_arr[ top].node_v;

		//Keyword constraint (for all its entries).
		//if( !( bst_node_v = bst_search( node_v->bst_v, key)))
			//continue;
		if( ( p_list = has_key_node( node_v, key)) == 0)
			continue;
		
		//bst_node_v != NULL.
		for( i=0; i<node_v->num; i++)
		{
			//Keyword constraint (for a specific entry).
			if( !get_k_bit( p_list, i))
				continue;

			//the i'th entry contains the keyword.
			if( node_v->level == 0)
			{
				//node_v is a leaf-node.
				MBR = ( ( obj_t*)( node_v->child[ i]))->MBR;
				
				//Disk constraint (if any, i.e., disk_v != NULL).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Update the current best solution.
					min_dist = c_dist;
					obj_v = ( obj_t*)( node_v->child[ i]);
				}
			}
			else
			{
				//node_v is an inner-node.
				MBR = node_v->MBRs[ ( ( node_t*)( node_v->child[ i]))->loc];

				//Disk constraint (if any).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Enqueue the child.
					b_h->obj_arr[ rear].node_v = ( node_t*)( node_v->child[ i]);
					b_h->obj_arr[ rear].key = c_dist;

					b_h_insert( b_h, rear++);
				}					
			}
		}//for
	}//while

	release_b_heap( b_h);
		
	return obj_v;
}

/*
 *	The kNN version of "const_NN_key".
 *
 *	Simply use a list for storing the kNN results.
 */
obj_set_t* const_k_NN_key( loc_t* loc_v, KEY_TYPE key, disk_t* disk_v, int k)
{	
	int size, top, i, rear, cnt;
	BIT_TYPE p_list;
	B_KEY_TYPE min_dist, c_dist;
	range* MBR;
	//obj_t* obj_v;
	node_t* node_v;
	b_heap_t* b_h;
	obj_node_t* obj_node_v, *obj_node_v1, *obj_node_rear;
	obj_set_t* obj_set_v;

	//Checking 1.
	if( IRTree_v.root->num == 0)
		return NULL;

	//Checking 2: involved in the search process implicitly.
	//if( !( bst_node_v = bst_search( IRTree.root->bst_v, key)))
		//return NULL;

	//Checking 3:
	MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
	if( disk_v != NULL && !is_overlap( MBR, disk_v))
		return NULL;

	//Initialize a heap.
	size = IRTree_v.obj_n / M + 1;
	b_h = alloc_b_heap( size);

	rear = 1;
	b_h->obj_arr[ rear].node_v = IRTree_v.root;
	b_h->obj_arr[ rear].key = calc_minDist( MBR, loc_v);
	free( MBR);

	/*s*/
	stat_v.memory_v -= sizeof( IRTree_v.dim * sizeof( range));
	/*s*/

	b_h_insert( b_h, rear++);


	obj_set_v = alloc_obj_set( );
	min_dist = DBL_MAX;
	//obj_v = NULL;

	while( !b_h_is_empty( b_h))
	{
		top = b_h_get_top( b_h);
		if( b_h->obj_arr[ top].key >= min_dist)
			break;

		node_v = b_h->obj_arr[ top].node_v;

		//Keyword constraint (for all its entries).
		//if( !( bst_node_v = bst_search( node_v->bst_v, key)))
			//continue;
		if( ( p_list = has_key_node( node_v, key)) == 0)
			continue;
		
		//bst_node_v != NULL.
		for( i=0; i<node_v->num; i++)
		{
			//Keyword constraint (for a specific entry).
			if( !get_k_bit( p_list, i))
				continue;

			//the i'th entry contains the keyword.
			if( node_v->level == 0)
			{
				//node_v is a leaf-node.
				MBR = ( ( obj_t*)( node_v->child[ i]))->MBR;
				
				//Disk constraint (if any, i.e., disk_v != NULL).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Update the current best solution.
					/*
					min_dist = c_dist;
					obj_v = ( obj_t*)( node_v->child[ i]);
					*/
					obj_node_v = ( obj_node_t*)malloc( sizeof( obj_node_t));
					memset( obj_node_v, 0, sizeof( obj_node_t));

					obj_node_v->obj_v = ( obj_t*)( node_v->child[ i]);
					obj_node_v->dist = c_dist;

					obj_node_v1 = obj_set_v->head;
					while( obj_node_v1->next != NULL && obj_node_v1->next->dist < c_dist)
						obj_node_v1 = obj_node_v1->next;

					obj_node_v->next = obj_node_v1->next;
					obj_node_v1->next = obj_node_v;
					obj_set_v->obj_n ++;

					//Update min_dist.
					cnt = 0;
					obj_node_rear = obj_set_v->head;
					while( obj_node_rear->next != NULL && cnt < k)
					{
						obj_node_rear = obj_node_rear->next;

						cnt++;
					}
					min_dist = obj_node_rear->dist;

					if( obj_node_rear->next != NULL)
					{
						free( obj_node_rear->next);
						obj_node_rear->next = NULL;
						obj_set_v->obj_n --;
					}
				}
			}
			else
			{
				//node_v is an inner-node.
				MBR = node_v->MBRs[ ( ( node_t*)( node_v->child[ i]))->loc];

				//Disk constraint (if any).
				if( disk_v != NULL && !is_overlap( MBR, disk_v))
					continue;

				c_dist = calc_minDist( MBR, loc_v);
				if( c_dist < min_dist)
				{
					//Enqueue the child.
					b_h->obj_arr[ rear].node_v = ( node_t*)( node_v->child[ i]);
					b_h->obj_arr[ rear].key = c_dist;

					b_h_insert( b_h, rear++);
				}					
			}
		}//for
	}//while

	release_b_heap( b_h);
		
	return obj_set_v;
}

/*
 *	Check whether a MBR @MBR is enclosed entirely by a disk @disk_v.
 */
bool is_enclosed( range* MBR, disk_t* disk_v)
{
	if( calc_maxDist( MBR, disk_v->center) <= disk_v->radius)
		return true;
	else
		return false;
}

/*
 *	Add an object entry @obj_v to @obj_set_v.
 */
void add_obj_set_entry( obj_t* obj_v, obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v;

	obj_node_v = ( obj_node_t*)malloc( sizeof( obj_node_t));
	memset( obj_node_v, 0, sizeof( obj_node_t));

	obj_node_v->obj_v = obj_v;
	obj_node_v->next =  obj_set_v->head->next;
	obj_set_v->head->next = obj_node_v;
	obj_set_v->obj_n ++;

	/*s*/
	stat_v.memory_v += sizeof( obj_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 *	Remove the first entry from the list of objects @obj_set_v.
 */	
void remove_obj_set_entry( obj_set_t* obj_set_v)
{
	obj_node_t* obj_node_v;

	obj_node_v = obj_set_v->head->next;
	obj_set_v->head->next = obj_node_v->next;
	obj_set_v->obj_n --;

	free( obj_node_v);

	/*s*/
	stat_v.memory_v -= sizeof( obj_node_t);
	/*s*/
}

/*
 *	Retrieve all the objects located at the sub-tree rooted at @node_v.
 *	The retrieved objects are stored in obj_set_v.
 */
void retrieve_sub_tree( node_t* node_v, obj_set_t* &obj_set_v, query_t* q)
{
	int i;
	BIT_TYPE p_list;

	if( node_v->level == 0)
	{
		//node_v is a leaf-node.
		//Retrieve all its objects.
		for( i=0; i<node_v->num; i++)
		{
			if( is_relevant_obj( ( obj_t*)( node_v->child[ i]), q))
				add_obj_set_entry( ( obj_t*)( node_v->child[ i]), obj_set_v);			
		}
	}
	else
	{
		//node_v is an inner-node.
		//Invoke the function recursively.
		p_list = is_relevant_node( node_v, q);
		for( i=0; i<node_v->num; i++)
		{
			if( get_k_bit( p_list, i))
				retrieve_sub_tree( ( node_t*)( node_v->child[ i]), obj_set_v, q);
		}
	}
}

/*
 *	Range query on the sub-tree rooted at @node_v.
 *	@disk_v indicates the range which is a circle.
 *
 *	The results are stored in @obj_set_v.
 */
void range_query_sub( node_t* node_v, disk_t* disk_v, obj_set_t* &obj_set_v, query_t* q)
{
	int i;
	BIT_TYPE p_list;
	range* MBR;

	if( node_v->parent == NULL)
		MBR = get_MBR_node( node_v, IRTree_v.dim);
	else
		MBR = node_v->parent->MBRs[ node_v->loc];

	//No overlapping.
	if( !is_overlap( MBR, disk_v))
		return;

	//Enclosed entrely.
	if( is_enclosed( MBR, disk_v))
	{
		retrieve_sub_tree( node_v, obj_set_v, q);
		if( node_v->parent == NULL)
		{
			free( MBR);

			/*s*/
			stat_v.memory_v -= IRTree_v.dim * sizeof( range);
			/*s*/
		}

		return;
	}

	//The remaining cases.
	if( node_v->level == 0)
	{
		//node_v is a leaf-node.
		for( i=0; i<node_v->num; i++)
		{
			if( is_enclosed( ( ( obj_t*)( node_v->child[ i]))->MBR, disk_v) &&
				is_relevant_obj( ( obj_t*)( node_v->child[ i]), q))
				add_obj_set_entry( ( obj_t*)( node_v->child[ i]), obj_set_v);
		}
	}
	else
	{
		//node_v is an inner-node.
		p_list = is_relevant_node( node_v, q);
		
		for( i=0; i<node_v->num; i++)
		{
			if( get_k_bit( p_list, i))
				range_query_sub( ( node_t*)( node_v->child[ i]), disk_v, obj_set_v, q);
		}
	}

	if( node_v->parent == NULL)
	{
		free( MBR);
		/*s*/
		stat_v.memory_v -= IRTree_v.dim * sizeof( range);
		/*s*/
	}
}

/*
 *	Circle range query.
 *
 *	DFS: recursive implementation.
 */
obj_set_t* range_query( disk_t* disk_v, query_t* q)
{
	obj_set_t* obj_set_v;

	obj_set_v = alloc_obj_set( );
	range_query_sub( IRTree_v.root, disk_v, obj_set_v, q);

	return obj_set_v;
}

/*
 *	Exclusion the objects on the boundary of disk @disk_v from @obj_set_v.
 */
void refine_region( obj_set_t* obj_set_v, disk_t* disk_v)
{
	B_KEY_TYPE dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = obj_set_v->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		loc_v = get_obj_loc( obj_node_v2->obj_v);
		dist = calc_dist_loc( loc_v, disk_v->center);
		
		release_loc( loc_v);
		
		if( dist == disk_v->radius)
		{
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			obj_set_v->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/
			
			continue;
		}

		obj_node_v1 = obj_node_v2;
		obj_node_v2 = obj_node_v1->next;
	}
}

/*
 *	construct a psi_t structure based on an object @obj_v.
 */
psi_t* get_psi_obj( obj_t* obj_v)
{
	psi_t* psi_v;
	k_node_t* k_node_v, *k_node_v1;

	psi_v = alloc_psi( );
	k_node_v = psi_v->k_head;

	k_node_v1 = obj_v->k_head->next;
	while( k_node_v1 != NULL)
	{	
		add_keyword_entry( k_node_v, k_node_v1->key);

		k_node_v1 = k_node_v1->next;
	}

	return psi_v;
}

/*
 *	Construct a psi_t structure based on a list of keywords @k_head.
 */
psi_t* const_psi( k_node_t* k_head)
{
	psi_t* psi_v;
	k_node_t* k_node_v;

	psi_v = ( psi_t*)malloc( sizeof( psi_t));
	memset( psi_v, 0, sizeof( psi_t));

	/*s*/
	stat_v.memory_v += sizeof( psi_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	psi_v->k_head = k_head;
	
	k_node_v = k_head->next;
	while( k_node_v != NULL)
	{
		psi_v->key_n ++;
		
		k_node_v = k_node_v->next;
	}

	return psi_v;
}

/*
 *	Exclude the keywords that occur in @k_head2 from @k_head1.
 *
 *	Return the resulting keywords in @k_head1.
 */
k_node_t* key_exclusion( k_node_t* k_head1, k_node_t* k_head2)
{
	int tag;
	k_node_t* k_node_v, *k_node_v1, *k_node_v2, *k_node_v3;

	k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( k_node_v, 0, sizeof( k_node_t));
	k_node_v3 = k_node_v;

	/*s*/
	stat_v.memory_v += sizeof( k_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	k_node_v1 = k_head1->next;
	while( k_node_v1 != NULL)
	{
		tag = 0;
		k_node_v2 = k_head2->next;
		while( k_node_v2 != NULL)
		{
			if( k_node_v2->key == k_node_v1->key)
			{
				tag = 1;
				break;
			}

			k_node_v2 = k_node_v2->next;
		}

		if( tag == 0)
		{
			//The current keyword should not excluded.
			add_keyword_entry( k_node_v3, k_node_v1->key);
		}

		k_node_v1 = k_node_v1->next;
	}

	return k_node_v;
}

/*
 *	Exclude the keywords that occur in @psi_v2 from @psi_v1.
 *
 *	Return the resulting keywords in @psi_v1.
 */
psi_t* psi_exclusion( psi_t* psi_v1, psi_t* psi_v2)
{
	int tag;
	k_node_t* k_node_v, *k_node_v1, *k_node_v2;
	psi_t* psi_v;

	psi_v = alloc_psi( );

	k_node_v = psi_v->k_head;
	k_node_v1 = psi_v1->k_head->next;
	while( k_node_v1 != NULL)
	{
		tag = 0;
		k_node_v2 = psi_v2->k_head->next;
		while( k_node_v2 != NULL)
		{
			if( k_node_v2->key == k_node_v1->key)
			{
				tag = 1;
				break;
			}

			k_node_v2 = k_node_v2->next;
		}

		if( tag == 0)
		{
			//The current keyword should not excluded.
			add_keyword_entry( k_node_v, k_node_v1->key);
			psi_v->key_n ++;
		}

		k_node_v1 = k_node_v1->next;
	}

	return psi_v;
}

/*
 *	Check whether the object set @O_t coveres the keyword set @psi.
 */
bool FeasibilityCheck( obj_set_t* O_t, psi_t* psi)
{
	int tag;
	KEY_TYPE key;
	k_node_t* k_node_v, *k_node_v1;
	obj_node_t* obj_node_v;

	k_node_v = psi->k_head->next;
	while( k_node_v != NULL)
	{
		key = k_node_v->key;

		tag = 0;
		obj_node_v =  O_t->head->next;
		while( obj_node_v != NULL)
		{
			k_node_v1 = obj_node_v->obj_v->k_head->next;
			while( k_node_v1 != 0)
			{
				if( k_node_v1->key == key)
				{
					tag = 1;
					break;
				}

				k_node_v1 = k_node_v1->next;
			}

			if( tag == 1)
				break;

			obj_node_v = obj_node_v->next;
		}

		if( tag == 0)
			return false;
		
		k_node_v = k_node_v->next;
	}

	return true;
}

/*
 *	Filter out the objects that are not located in range @disk_v from @O_t.
 */
void obj_filter_range( obj_set_t* &O_t, disk_t* disk_v)
{
	B_KEY_TYPE dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = O_t->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		loc_v = get_obj_loc( obj_node_v2->obj_v);
		dist = calc_dist_loc( loc_v, disk_v->center);
		if( dist > disk_v->radius)
		{
			//the obj is not located in the disk.
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			O_t->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/
		}
		else
		{
			obj_node_v1 = obj_node_v2;
			obj_node_v2 = obj_node_v2->next;
		}

		release_loc( loc_v);
	}		
}

/*
 *	Construct the IF on a set of objects @obj_set_v for the keywords in @psi_v.
 *
 *	1. The IF structure is indexed by a binary search tree.
 *	2. No ordering is imposed in IF.
 */
bst_t* const_IF( obj_set_t* obj_set_v, psi_t* psi_v)
{
	int i;
	bst_t* IF_v;
	k_node_t* k_node_v;
	obj_node_t* obj_node_v;
	bst_node_t* bst_node_v;

	//IF_v = alloc_IF( psi_v->key_n);
	IF_v = bst_ini( );

	k_node_v = psi_v->k_head->next;
	for( i=0; i<psi_v->key_n; i++)
	{
		//IF_v->entry_v[ i].key = k_node_v->key;
		bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
		memset( bst_node_v, 0, sizeof( bst_node_t));

		/*s*/
		stat_v.memory_v += sizeof( bst_node_t);
		/*s*/

		bst_node_v->key = k_node_v->key;
		bst_node_v->p_list_obj = alloc_obj_set( );
		
		obj_node_v = obj_set_v->head->next;
		while( obj_node_v != NULL)
		{
			if( has_key_obj( obj_node_v->obj_v, k_node_v->key))
				add_obj_set_entry( obj_node_v->obj_v, bst_node_v->p_list_obj);

			obj_node_v = obj_node_v->next;
		}

		bst_insert( IF_v, bst_node_v);

		k_node_v = k_node_v->next;
	}
	
	/*s*/
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return IF_v;
}

/*
 *	Construct a obj_set_t structure containing 3 objects
 *	from a tri_t structure @triplet_v.
 */
obj_set_t* const_obj_set( tri_t* triplet_v)
{
	obj_set_t* obj_set_v;

	obj_set_v = alloc_obj_set( );
	obj_set_v->obj_n = 3;

	//Including the objects.
	add_obj_set_entry( triplet_v->o, obj_set_v);
	add_obj_set_entry( triplet_v->o_1, obj_set_v);
	add_obj_set_entry( triplet_v->o_2, obj_set_v);

	return obj_set_v;
}

/*
 *	Check the distance constraint.
 */
bool check_dist_constraint( obj_set_t* obj_set_v, obj_t* obj_v, obj_t* o, B_KEY_TYPE d)
{
	obj_node_t* obj_node_iter;

	if( calc_dist_obj( o, obj_v) > d)
		return false;

	obj_node_iter = obj_set_v->head->next;
	while( obj_node_iter != NULL)
	{
		if( calc_dist_obj( obj_node_iter->obj_v, obj_v) > d)
			return false;

		obj_node_iter = obj_node_iter->next;
	}

	return true;
}

/*
 *	Update the IF structure @IF_v by removing the keywords 
 *	that have been covered by an object @obj_v.
 */
bst_node_list_t* update_IF_obj( bst_t* IF_v, obj_t* obj_v)
{
	k_node_t* k_node_v;
	bst_node_t* bst_node_v;
	bst_node_list_t* bst_node_list_v, *tmp;

	bst_node_list_v = ( bst_node_list_t*)malloc( sizeof( bst_node_list_t));
	memset( bst_node_list_v, 0, sizeof( bst_node_list_t));

	/*s*/
	stat_v.memory_v += sizeof( bst_node_list_t);
	/*s*/

	k_node_v = obj_v->k_head->next;
	while( k_node_v != NULL)
	{
		bst_node_v = bst_search( IF_v, k_node_v->key);

		if( bst_node_v == NULL)
		{
			k_node_v = k_node_v->next;
			continue;
		}

		bst_delete( IF_v, bst_node_v);
		//bug.
		bst_node_v->p = NULL;
		bst_node_v->left = NULL;
		bst_node_v->right = NULL;
		
		tmp = ( bst_node_list_t*)malloc( sizeof( bst_node_list_t));
		memset( tmp, 0, sizeof( bst_node_list_t));

		/*s*/
		stat_v.memory_v += sizeof( k_node_t);
		/*s*/

		tmp->bst_node_v = bst_node_v;
		tmp->next = bst_node_list_v->next;
		bst_node_list_v->next = tmp;

		k_node_v = k_node_v->next;
	}

	/*s*/
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	return bst_node_list_v;
}

/*
 *	Release a bst_node_list_t structure.
 */
void release_bst_node_list( bst_node_list_t* bst_node_list_v)
{
	bst_node_list_t* tmp;
	
	while( bst_node_list_v != NULL)
	{
		tmp = bst_node_list_v->next;
		free( bst_node_list_v);
		bst_node_list_v = tmp;

		/*s*/
		stat_v.memory_v -= sizeof( bst_node_list_t);
		/*s*/
	}
}

/*
 *	Restore the IF_v structure @IF_v by re-including the bst_nodes of @bst_node_list_v.
 */
void restore_IF_bst_node_list( bst_t* IF_v, bst_node_list_t* bst_node_list_v)
{
	bst_node_list_t* bst_node_list_iter;

	bst_node_list_iter = bst_node_list_v->next;
	while( bst_node_list_iter != NULL)
	{
		bst_insert( IF_v, bst_node_list_iter->bst_node_v);

		bst_node_list_iter = bst_node_list_iter->next;
	}
}

/*
 *	Combine two obj_set_t structures @obj_set_v1 and @obj_set_v2.
 *		1. The combined result is stored in @obj_set_v1.
 *		2. The obj_set_v2 is released.
 *
 */
void combine_obj_set( obj_set_t* obj_set_v1, obj_set_t* obj_set_v2)
{
	obj_node_t* obj_node_v;

	//Locate the rear of obj_set_v1.
	obj_node_v = obj_set_v1->head;
	while( obj_node_v->next != NULL)
		obj_node_v = obj_node_v->next;

	obj_node_v->next = obj_set_v2->head->next;
	obj_set_v2->head->next = NULL;
	obj_set_v1->obj_n += obj_set_v2->obj_n;

	release_obj_set( obj_set_v2);
}

/*
 *	The sub-procedure of function "const_feasible_set_sub".
 */
obj_set_t* const_feasible_set_sub( bst_t* IF_v, obj_set_t* S_0, obj_t* o, B_KEY_TYPE d)
{
	obj_t* obj_v;
	obj_set_t* S;
	bst_node_t* bst_node_v;
	obj_node_t* obj_node_v;
	bst_node_list_t* bst_node_list_v;

	if( IF_v->node_n == 0)
		return alloc_obj_set( );			//An empty one.

	bst_node_v = IF_v->root;
	obj_node_v = bst_node_v->p_list_obj->head->next;
	while( obj_node_v != NULL)
	{
		//Pick an object.
		obj_v = obj_node_v->obj_v;

		//Distance constraint checking.
		if( !check_dist_constraint( S_0, obj_v, o, d))
		{
			obj_node_v = obj_node_v->next;
			continue;
		}
		
		//Update the IF_v.
		bst_node_list_v = update_IF_obj( IF_v, obj_v);

		//Update the S_0.
		//obj_v is added at the first place of S_0.
		add_obj_set_entry( obj_v, S_0);

		//Sub-procedure.
		S = const_feasible_set_sub( IF_v, S_0, o, d);

		//Restore the S_0.
		remove_obj_set_entry( S_0);

		//Restore the IF_v.
		restore_IF_bst_node_list( IF_v, bst_node_list_v);

		release_bst_node_list( bst_node_list_v);
		
		//Checking.
		if( S != NULL)
		{	
			//Include obj_v into S.
			add_obj_set_entry( obj_v, S);

			return S;
		}

		//S == NULL.

		//Try the next object candidate.
		obj_node_v = obj_node_v->next;
	}

	return NULL;
}

/*
 *	Check whether there exists a subset of @O_t such that
 *	1. it covers @psi and 2. its diameter is d.
 *
 *	Method: recursive implementation.
 *
 *	If so, return such a subset; otherwise, return NULL.
 */
obj_set_t* const_feasible_set( obj_set_t* O_t, psi_t* psi, obj_t* o, B_KEY_TYPE d)
{
	bst_t* IF_v;
	obj_set_t* S_0, *S;

	//Construct an IF structure.
	IF_v = const_IF( O_t, psi);

	//Initialize the S_0.
	S_0 = alloc_obj_set( );

	//Invoke the sub-procedure "recursively".
	S = const_feasible_set_sub( IF_v, S_0, o, d);

	//Release the resources.
	bst_release( IF_v);
	release_obj_set( S_0);

	return S;
}

/*
 *	The implementation of the "ConstructFeasibleSet-Appro" procedure in the paper.
 */
obj_set_t* ConstructFeasibleSet_Appro( obj_t* o, query_t* q)
{
	B_KEY_TYPE radius;

	obj_t* obj_v;
	obj_set_t* S;
	k_node_t* k_head1, *k_node_v;
	loc_t* loc_v, *loc_v1;
	disk_t* disk_v;
	
	S = alloc_obj_set( );

	//Location of object o.
	loc_v = get_obj_loc( o);

	//Disk(q, d(o, q)).
	loc_v1 = copy_loc( q->loc_v);
	radius = calc_dist_loc( loc_v, q->loc_v);
	disk_v = const_disk( loc_v1, radius);

	//Include object o in S.
	add_obj_set_entry( o, S);

	//Obtain the "un-covered" keywords by S.
	k_head1 = key_exclusion( q->psi_v->k_head, o->k_head);

	//Retrieve the disk-constrained NNs that contain the un-covered keywords.
	k_node_v = k_head1->next;
	while( k_node_v != NULL)
	{
		obj_v = const_NN_key( loc_v, k_node_v->key, disk_v);
		
		if( obj_v == NULL)
		{
			release_obj_set( S);
			S = NULL;
			break;
		}

		add_obj_set_entry( obj_v, S);

		k_node_v = k_node_v->next;
	}

	//Release the memory.
	release_loc( loc_v);
	release_disk( disk_v);
	release_k_list( k_head1);
	
	return S;
}

/*
 *	Compute the farthest distance between a set of objects @obj_set_v and @q.
 */
B_KEY_TYPE comp_farthest( obj_set_t* obj_set_v, query_t* q)
{
	B_KEY_TYPE far_dist, dist;
	loc_t* loc_v;
	obj_node_t* obj_node_v;

	far_dist = 0;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		
		dist = calc_dist_loc( loc_v, q->loc_v);
		if( dist > far_dist)
			far_dist = dist;

		release_loc( loc_v);
		obj_node_v = obj_node_v->next;
	}

	return far_dist;
}

/*
 *	Compute the diameter of a set of objects @obj_set_v.
 */
B_KEY_TYPE comp_diameter( obj_set_t* obj_set_v)
{
	B_KEY_TYPE dia, dist;
	obj_node_t* obj_node_v1, *obj_node_v2;

	if( obj_set_v->obj_n <= 1)
		return 0;

	//There exist at least two objects.
	dia = 0;
	obj_node_v1 = obj_set_v->head->next;
	while( obj_node_v1->next != NULL)
	{
		obj_node_v2 = obj_node_v1->next;
		while( obj_node_v2 != NULL)
		{			
			dist = calc_dist_obj( obj_node_v1->obj_v, obj_node_v2->obj_v);
			if( dist > dia)
				dia = dist;

			obj_node_v2 = obj_node_v2->next;
		}

		obj_node_v1 = obj_node_v1->next;
	}

	return dia;		
}

/*
 *	Compute the cost of set @obj_set_v wrt @q.
 */
B_KEY_TYPE comp_cost( obj_set_t* obj_set_v, query_t* q)
{
	if( cost_tag == 1) 
		return comp_farthest( obj_set_v, q) + comp_diameter( obj_set_v);
	else
		return max( comp_farthest( obj_set_v, q), comp_diameter( obj_set_v));		
}

/*
 *	Compute the lower bound @LB and the upper bound @UB of cost(S*, q).c_1,
 *	where S* is the optimal solution.
 *
 *	return NN(q, q.\psi) if successful; otherwise, return NULL.
 */
obj_set_t* comp_bounds( query_t* q, B_KEY_TYPE &LB, B_KEY_TYPE &UB)
{
	B_KEY_TYPE dist;
	k_node_t* k_node_v;
	obj_set_t* obj_set_v;
	loc_t* loc_v;
	obj_t* obj_v;

	obj_set_v = alloc_obj_set( );

	//Compute LB.
	LB = 0;
	k_node_v = q->psi_v->k_head->next;
	while( k_node_v != NULL)
	{
		obj_v = const_NN_key( q->loc_v, k_node_v->key, NULL);
		
		if( obj_v == NULL)
		{
			release_obj_set( obj_set_v);
			return NULL;
		}

		loc_v = get_obj_loc( obj_v);
		dist = calc_dist_loc( loc_v, q->loc_v);
		if( dist > LB)
			LB = dist;
		
		add_obj_set_entry( obj_v, obj_set_v);

		k_node_v = k_node_v->next;
	}

	//Compute the UB.
	UB = comp_cost( obj_set_v, q);

	return obj_set_v;
	//release_obj_set( obj_set_v);
}

/*
 *	Exclude from a set of objects @obj_set_v
 *	those objects that are located in a specified disk @disk_v.
 *
 *	Note that those objects on the boundary of the disk are kept.
 */
void obj_exclusion_disk( obj_set_t* obj_set_v, disk_t* disk_v)
{
	B_KEY_TYPE far_dist;
	obj_node_t* obj_node_v1, *obj_node_v2;

	obj_node_v1 = obj_set_v->head;
	obj_node_v2 = obj_node_v1->next;
	while( obj_node_v2 != NULL)
	{
		far_dist = calc_maxDist( obj_node_v2->obj_v->MBR, disk_v->center);
		if( far_dist < disk_v->radius)
		{
			//The object should be excluded.
			obj_node_v1->next = obj_node_v2->next;
			free( obj_node_v2);
			obj_node_v2 = obj_node_v1->next;

			obj_set_v->obj_n --;

			/*s*/
			stat_v.memory_v -= sizeof( obj_node_t);
			/*s*/

			continue;
		}

		obj_node_v1 = obj_node_v2;
		obj_node_v2 = obj_node_v2->next;
	}
}

/*
 *	Sort the objects in @obj_set_v by their distances to @q.
 *
 *	Method: heap-sort.
 *	Alternative method: based on the binary search tree.
 */
b_heap_t* heap_sort_obj_set( obj_set_t* obj_set_v, query_t* q)
{
	int cur;
	B_KEY_TYPE dist;
	b_heap_t* b_h;
	obj_node_t* obj_node_v;
	loc_t* loc_v;

	b_h = alloc_b_heap( obj_set_v->obj_n + 1);
	
	cur = 1;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		dist = calc_dist_loc( loc_v, q->loc_v);
		release_loc( loc_v);

		b_h->obj_arr[ cur].key = dist;
		b_h->obj_arr[ cur].obj_v = obj_node_v->obj_v;

		b_h_insert( b_h, cur);
		
		cur ++;

		obj_node_v = obj_node_v->next;
	}

	return b_h;
}
b_heap_t* heap_sort_obj_set_y( obj_set_t* obj_set_v)
{
	int cur;
	B_KEY_TYPE dist;
	b_heap_t* b_h;
	obj_node_t* obj_node_v;
	loc_t* loc_v;

	b_h = alloc_b_heap( obj_set_v->obj_n + 1);
	
	cur = 1;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		dist = loc_v->coord[1];
		release_loc( loc_v);

		b_h->obj_arr[ cur].key = dist;
		b_h->obj_arr[ cur].obj_v = obj_node_v->obj_v;

		b_h_insert( b_h, cur);
		
		cur ++;

		obj_node_v = obj_node_v->next;
	}

	return b_h;
}
b_heap_t* heap_sort_obj_set_x( obj_set_t* obj_set_v)
{
	int cur;
	B_KEY_TYPE dist;
	b_heap_t* b_h;
	obj_node_t* obj_node_v;
	loc_t* loc_v;

	b_h = alloc_b_heap( obj_set_v->obj_n + 1);
	
	cur = 1;
	obj_node_v = obj_set_v->head->next;
	while( obj_node_v != NULL)
	{
		loc_v = get_obj_loc( obj_node_v->obj_v);
		dist = loc_v->coord[0];
		release_loc( loc_v);

		b_h->obj_arr[ cur].key = dist;
		b_h->obj_arr[ cur].obj_v = obj_node_v->obj_v;

		b_h_insert( b_h, cur);
		
		cur ++;

		obj_node_v = obj_node_v->next;
	}

	return b_h;
}

/*
 *	The implementation of the "AchievabilityCheck" function in the paper.
 */
obj_set_t* AchievabilityCheck( tri_t triplet_v, query_t* q)
{
	//Disks.
	B_KEY_TYPE radius_1, radius_2;
	loc_t* loc_v1, *loc_v2, *loc_v3;
	disk_t* disk_v1, *disk_v2, *disk_v3, *disk_v_tmp;
	obj_set_t* O_t, *S;
	k_node_t* k_head_1, *k_head_2, *k_head_3;
	psi_t* psi_v;

	//Pre-checking.
	//Exclude the keywords covered by the triplet.
	k_head_1 = key_exclusion( q->psi_v->k_head, triplet_v.o->k_head);
	k_head_2 = key_exclusion( k_head_1, triplet_v.o_1->k_head);
	k_head_3 = key_exclusion( k_head_2, triplet_v.o_2->k_head);

	psi_v = const_psi( k_head_3);

	release_k_list( k_head_1);
	release_k_list( k_head_2);
	//Note that k_head_3 would be released when psi_v is released.

	if( psi_v->key_n == 0)
	{
		//<o, o_1, o_2> is achievable.
		//S = {o, o_1, o_2} is a verying feasible set.
		S = const_obj_set( &triplet_v);

		release_psi( psi_v);
		return S;		
	}

	//Normal-checking.
	//Construct the disks.
	//1. Disk(q, d(o, q)).
	loc_v1 = get_obj_loc( triplet_v.o);
	radius_1 = calc_dist_loc( loc_v1, q->loc_v);
	disk_v1 = const_disk( q->loc_v, radius_1);

	//2. Disk(o_1, d(o_1, o_2)) and Disk(o_2, d(o_1, o_2)).
	loc_v2 = get_obj_loc( triplet_v.o_1);
	loc_v3 = get_obj_loc( triplet_v.o_2);
	radius_2 = calc_dist_loc( loc_v2, loc_v3);
	
	disk_v2 = const_disk( loc_v2, radius_2);
	disk_v3 = const_disk( loc_v3, radius_2);

	//Decide the disk for range query.
	//Use disk_v1 at default.
	if( radius_2 < radius_1)
	{
		//Change to use disk_v2 or disk_v3.
		disk_v_tmp = disk_v1;
		disk_v1 = disk_v2;
		disk_v2 = disk_v_tmp;
	}

	//Range query on disk_v1.
	O_t = range_query( disk_v1, q);
	
	//Filter the objs that are not located in disk_v2 or disk_v3.
	obj_filter_range( O_t, disk_v2);
	obj_filter_range( O_t, disk_v3);

	release_loc( loc_v1);
	release_loc( loc_v2);
	release_loc( loc_v3);
	release_disk( disk_v1);
	release_disk( disk_v2);
	release_disk( disk_v3);
	
	//Check whether O_t covers the keywords in the query.
	if( !FeasibilityCheck( O_t, q->psi_v))
	{
		S = NULL;
		goto E;
	}

	/*s*/
	stat_v.O_t_size_sum += O_t->obj_n;
	/*s*/

	//Find a sub-set of O_t which covers psi_v.
	S = const_feasible_set( O_t, psi_v, triplet_v.o, radius_2);

	if( S == NULL)
		goto E;

	//Combine the sub-set with the critical objects.
	add_obj_set_entry( triplet_v.o, S);
	add_obj_set_entry( triplet_v.o_1, S);
	add_obj_set_entry( triplet_v.o_2, S);

	/*s*/
	stat_v.psi_n_sum += psi_v->key_n;
	/*s*/

E:	
	//Release the resources.
	release_psi( psi_v);
	release_obj_set( O_t);

	return S;
}

/*
 *	Handle one specific candidate of <o_1, o_2>.
 */
void process_obj_pair( obj_t* o_1, obj_t* o_2, B_KEY_TYPE dist, bst_t* obj_pair)
{
	bst_node_t* bst_node_v;	

	//Add the pair into the bst.
	bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
	memset( bst_node_v, 0, sizeof( bst_node_t));

	/*s*/
	stat_v.memory_v += sizeof( bst_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/

	bst_node_v->key = dist;
	bst_node_v->obj_v1 = o_1;
	bst_node_v->obj_v2 = o_2;

	bst_insert( obj_pair, bst_node_v);
}

/*
 *	Initialize the set of <o_1, o_2> pairs.
 */
bst_t* ini_obj_pair( obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q)
{
	B_KEY_TYPE d_low, d_high, dist, d_o_q;
	obj_t* o_1, *o_2;
	bst_t* obj_pair;
	obj_node_t* obj_node_v1, *obj_node_v2;
	loc_t* loc_v, *loc_v1, *loc_v2;

	obj_pair = bst_ini( );

	if( region->obj_n == 0)
		return obj_pair;

	//Compute the upper bound of d(o_1, o_2): independent on <o_1, o_2>.
	loc_v = get_obj_loc( o);
	d_o_q = calc_dist_loc( loc_v, q->loc_v);

	if( cost_tag == 1)
		d_high = d_u - d_o_q;
	else
		d_high = d_u;
	release_loc( loc_v);

	obj_node_v1 = region->head->next;
	while( obj_node_v1 != NULL)
	{
		o_1 = obj_node_v1->obj_v;

		obj_node_v2 = obj_node_v1;
		while( obj_node_v2 != NULL)
		{
			o_2 = obj_node_v2->obj_v;

			//Compute the lower bound of d(o_1, o_2): dependent on <o_1, o_2>.
			if( cost_tag == 1)
			{
				loc_v1 = get_obj_loc( o_1);
				loc_v2 = get_obj_loc( o_2);
				d_low = d_l - min( calc_dist_loc( loc_v1, q->loc_v), 
					calc_dist_loc( loc_v2, q->loc_v));	
				
				release_loc( loc_v1);
				release_loc( loc_v2);
			}
			else
			{
				d_low = d_o_q;
			}

			//process.
			dist = calc_dist_obj( o_1, o_2);
			if( cost_tag == 1)
			{
				if( !( dist >= d_high || dist < d_low))
					process_obj_pair( o_1, o_2, dist, obj_pair);
			}
			else
			{
				if( !( dist >= d_high || dist <= d_low))
					process_obj_pair( o_1, o_2, dist, obj_pair);
			}

			obj_node_v2 = obj_node_v2->next;
		}

		obj_node_v1 = obj_node_v1->next;
	}

	return obj_pair;
}

/*
 *	Update the set of <o_1, o_2> pairs.
 *		1. region should also be updated.
 */
void update_obj_pair( bst_t* obj_pair, obj_set_t* region, obj_t* o, B_KEY_TYPE d_l, B_KEY_TYPE d_u, query_t* q)
{
	B_KEY_TYPE dist, d_low, d_high, d_o_q;
	obj_node_t* obj_node_v;
	loc_t* loc_v;
	
	//Compute the upper bound of d(o_1, o_2): independent on <o_1, o_2>.
	loc_v = get_obj_loc( o);
	d_o_q = calc_dist_loc( loc_v, q->loc_v);
	if( cost_tag == 1)
		d_high = d_u - d_o_q;
	else
		d_high = d_u;
	release_loc( loc_v);

	//Remove the unsatisfied obj_pairs in the original obj_pair_v.
	if( d_high <= obj_pair->max)
		bst_trim( obj_pair, d_high);	

	//Include o into region. bug.
	add_obj_set_entry( o, region);
	
	obj_node_v = region->head->next;
	while( obj_node_v != NULL)
	{
		//Compute the lower bound of d(o_1, o_2): dependent on <o_1, o_2>.
		//The lower bound is equal to 0.

		if( cost_tag == 1)
			d_low = 0;
		else
			d_low = d_o_q;
		
		//process.
		dist = calc_dist_obj( obj_node_v->obj_v, o);
		if( cost_tag == 1)
		{
			if( !( dist >= d_high || dist < d_low))
				process_obj_pair( obj_node_v->obj_v, o, dist, obj_pair);
		}
		else
		{
			if( !( dist >= d_high || dist <= d_low))
				process_obj_pair( obj_node_v->obj_v, o, dist, obj_pair);
		}

		obj_node_v = obj_node_v->next;
	}
}

/*
 *	The implementation of the "ConstructFeasibleSet-Exact" procedure in the paper CoSK.
 */
obj_set_t* ConstructFeasibleSet_Exact( obj_t* o, query_t* q, bst_t* obj_pair)
{
	int j;
	obj_set_t* S;
	tri_t triplet_v;
	bst_node_t* cur_pair;
	obj_t* o_q;

	if( cost_tag == 2)
	{
		//Pre-checking for the Cost2 measurement.
		//Create o_q.
		o_q = ( obj_t*)malloc( sizeof( obj_t));
		alloc_obj( o_q, q->loc_v->dim);
		for( j=0; j<q->loc_v->dim; j++)
		{
			o_q->MBR[ j].min = q->loc_v->coord[ j];
			o_q->MBR[ j].max = o_q->MBR[ j].min;
		}

		S = NULL;
		triplet_v.o = o;
		triplet_v.o_1 = o;
		triplet_v.o_2 = o_q;
		S =  AchievabilityCheck( triplet_v, q);

		/*s*/
		stat_v.achi_sum ++;
		/*s*/

		if( S != NULL)
			return S;
	}		

	//Upper bound pruning: involved in CostEnum.

	S = NULL;
	triplet_v.o = o;
	cur_pair = bst_get_min( obj_pair->root);
	while( cur_pair != NULL)
	{
		triplet_v.o_1 = cur_pair->obj_v1;
		triplet_v.o_2 = cur_pair->obj_v2;

		//Lower bound pruning:
		if( calc_dist_obj( triplet_v.o_1, triplet_v.o_2) < 
			max( calc_dist_obj( triplet_v.o, triplet_v.o_1), 
			calc_dist_obj( triplet_v.o, triplet_v.o_2)))
		{
			cur_pair = bst_successor( cur_pair);
			continue;
		}

		S = AchievabilityCheck( triplet_v, q);

		/*s*/
		stat_v.achi_sum ++;
		/*s*/
		
		if( S != NULL)
			return S;

		//S == NULL.
		cur_pair = bst_successor( cur_pair);
	}

	return S;
}

/*
 *	The implementation of "CostEnum".
 *
 *	s_tag = 1: CostEnum-Exact.
 *	s_tag = 2: CostEnum-Appro.
 */
obj_set_t* CostEnum( query_t* q, int s_tag, int prune_tag)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, cost, dist;
	obj_set_t* S_a;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_t* o, *o_next;
	loc_t* loc_v;
	bst_t* obj_pair;
	b_heap_t* R_heap;

	//Compute the LB and UB of cost(S*, q).c_1;
	S_a = comp_bounds( q, LB, UB);
	cost_c = UB;
	if( S_a == NULL)
		return NULL;

	//Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);
	
	region_u = range_query( disk_u, q);
	region_l = range_query( disk_l, q);
	refine_region( region_l, disk_l);

	/*t/
	print_obj_set( region_l, stdout);
	/*t*/

	obj_exclusion_disk( region_u, disk_l);

	R = region_u;

	/*t/
	print_obj_set( R, stdout);
	/*t*/

	/*s/
	printf( "#cands: %i\n", R->obj_n);
	/*s*/

	//Pre-checking.
	if( R->obj_n == 0)
		goto E;

	//Sort the objects in R by their distances to q.
	R_heap = heap_sort_obj_set( R, q);
	top = b_h_get_top( R_heap);

	o_next = R_heap->obj_arr[ top].obj_v;

	//Construct the base set of the <o_1, o_2> pairs.
	//region_l should be strict.
	if( s_tag == 1)
	{
		obj_pair = ini_obj_pair( region_l, o_next, LB, UB, q);

		/*s/
		printf( "#pairs: %i\n", obj_pair->node_n);
		/*s*/
	}

	//Search.
	while( true)
	{
		if( o_next == NULL)
			break;

		o = o_next;

		top = b_h_get_top( R_heap);
		if( top == 0)
			o_next = NULL;
		else
			o_next = R_heap->obj_arr[ top].obj_v;

		loc_v = get_obj_loc( o);
		dist = calc_dist_loc( loc_v, q->loc_v);
		release_loc( loc_v);

		if( dist > cost_c)
			break;

		/*s*/
		stat_v.n_1_sum ++;
		/*s*/

		//The "ConstructFeasibleSet" procedure.
		if( s_tag == 1)
		{
			//Pre-checking (for the boundary case tha |S| = 1).
			S = alloc_obj_set( );
			add_obj_set_entry( o, S);
			if( !is_covered_obj_set( S, q))
			{
				release_obj_set( S);

				if( prune_tag == 1)
				{
					//Update UB by using ConstructFeasibleSet-Appro.
					S = ConstructFeasibleSet_Appro( o, q);
					if( !S)
						continue;
					
					cost = comp_cost( S, q);
					if( cost < cost_c)
					{
						release_obj_set( S_a);
						S_a = S;
						cost_c = cost;
						
						UB = cost_c;
					}
					else
						release_obj_set( S);
				}

				//Update the <o_1, o_2> pairs.
				update_obj_pair( obj_pair, region_l, o, LB, UB, q);

				/*s/
				printf( "#pairs: %i\n", obj_pair->node_n);
				/*s*/

				/*t*/
				//in_order_walk( obj_pair->root);
				/*t*/

				S = ConstructFeasibleSet_Exact( o, q, obj_pair);
			}
		}
		else	//s_tag == 2
		{
			S = ConstructFeasibleSet_Appro( o, q);

			/*t/
			printf( "base object: %i\n", o->id);
			print_obj_set( S, stdout);
			/*t*/
		}

		if( !S)
			continue;

		cost = comp_cost( S, q);
		if( cost < cost_c)
		{
			release_obj_set( S_a);
			S_a = S;
			cost_c = cost;
		}
		else
			release_obj_set( S);
	}//while
E:

	remove_identical_obj( S_a);

	//Release the resource.
	release_disk( disk_u);
	release_disk( disk_l);
	release_obj_set( R);
	release_obj_set( region_l);
	if( s_tag == 1)
		bst_release( obj_pair);
	release_b_heap( R_heap);

	return S_a;
}

/*
 *	The implementation of the "CostEnum-Exact" algorithm.
 */
obj_set_t* CostEnum_Exact( query_t* q, int prune_tag)
{
	return CostEnum( q, 1, prune_tag);
}

/*
 *	The implementation of the "CostEnum-Appro" algorithm.*	
 */
obj_set_t* CostEnum_Appro( query_t* q)
{
	return CostEnum( q, 2, 0);
}


double Adapt1( query_ssq* q, int alg_opt, int prune_opt)
{
	int i, rear, top;
	B_KEY_TYPE costV, costV_1, dist;
	b_heap_t* U;
	obj_set_t *V_1;
	void* e;
	node_t* node_v;
	obj_t* obj_v;

	loc_t* loc_v;
	query_t* q_new;
	
	costV=INT_MAX;

	//printf( "Cao-Appro2:\n");

	U = alloc_b_heap( INI_HEAP_SIZE);
	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);

	b_h_insert( U, rear++);


	/*s*/
	stat_v.n_k ++;
	/*s*/

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
		//	if(dist>costV)break;
			for( i=0; i<node_v->num; i++)
			{
				
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
			double cdist=calc_dist_loc(q->loc_v,loc_v);
		//	if(cdist>costV)break;

			if(!has_key_obj(obj_v,q->target))continue;
			//Construct a new query instance.
			//q_new = alloc_query( );
			query_t* qnew;
			qnew=alloc_query();
			qnew->loc_v=alloc_loc(2);
			qnew->loc_v = copy_loc(loc_v);
			//qnew->target=q->target;
			qnew->psi_v = alloc_psi( );
			copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);
			printf("Adapt1\n");

			V_1=alloc_obj_set( );
			if(alg_opt==9)V_1 = Cao_Appro1( qnew);
				
			if(alg_opt==10)V_1= Cao_Appro2(qnew);
				
			if(alg_opt==11)V_1=Cao_Exact(qnew);
				
			if(alg_opt==12)V_1= CostEnum_Exact(qnew,prune_opt);
				
			if(alg_opt==13)V_1=CostEnum_Appro(qnew);
				
			//add_obj_set_entry( obj_v, V_1);
			
			
		//	costV_1 = comp_cost( V_1, q);
			costV_1=calc_dist_loc(q->loc_v,loc_v)+comp_cost(V_1,qnew);
			if( costV_1 < costV)
			{
				costV = costV_1;
		
			}
			release_query( qnew);		
			release_obj_set(V_1);
			/*s*/
			stat_v.n_k ++;
			/*s*/
		}//else
	}//while
	
	release_b_heap( U);

	return costV;
}


/*
 *	The implementation of "Approx_SSQ".
 *
 */
double Approx_SSQ( query_ssq* q)
{
//	int top;
//	B_KEY_TYPE LB, UB, cost_c, cost, dist;
//	obj_set_t* S_a;				//the current best solution.
//	obj_set_t* S;				//the newly constructed feasible set.
//	obj_set_t* R;				//region R.
//	disk_t* disk_u, *disk_l;	//the outer and inner disks.
//	obj_set_t* region_u, *region_l;
	obj_t* o, *o_next;
	loc_t* loc_v;
	bst_t* obj_pair;
	double cdist;
	b_heap_t* R_heap;
	int i, rear, top;
	KEY_TYPE c_key;
	B_KEY_TYPE costV_1, dist,costV;
	b_heap_t* U;
	obj_set_t* V, *V_1, *V_tmp;
	void* e;
	node_t* node_v;
	obj_t* obj_v;
	bst_node_t* bst_node_v;
	BIT_TYPE p_list;

	query_t* q_new;
	q_new=alloc_query();
	q_new->loc_v=copy_loc(q->loc_v);
	q_new->psi_v=alloc_psi();
	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(q_new->psi_v,q->target);
	vorc::Voronoi * v_head=new vorc::Voronoi();
	vorc::Voronoi * v_last=v_head;
	for(k_node_t * psi1=q->psi_v->k_head->next;psi1->next!=NULL;psi1=psi1->next){
		for(k_node_t * psi2=psi1->next;psi2->next!=NULL;psi2=psi2->next){
		//	if(psi2==NULL)break;
			for(k_node_t * psi3=psi2->next;;psi3=psi3->next){
				if(psi3==NULL)break;
				vorc::Voronoi * v=new vorc::Voronoi();
				v->key[0]=psi1->key;
				v->key[1]=psi2->key;
				v->key[2]=psi3->key;
				v_last->next=v;
				v_last=v_last->next;
				v->root=NULL;
		}
	}
}
	for(k_node_t * psi1=q->psi_v->k_head->next;psi1->next!=NULL;psi1=psi1->next){
		for(k_node_t * psi2=psi1->next;;psi2=psi2->next){
			if(psi2==NULL)break;
			vorc::Voronoi * v=new vorc::Voronoi();
			v->key[0]=q->target;
			v->key[1]=psi1->key;
			v->key[2]=psi2->key;
			v_last->next=v;
			v_last=v_last->next;
			v->root=NULL;
		}
	}
	v_last->next=NULL;
	//Search.
	U = alloc_b_heap( INI_HEAP_SIZE);

	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);

	b_h_insert( U, rear++);
	
	costV=100000000000000000.0;

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
			if (dist>costV)break;

			if(!is_relevant_node( node_v,q_new))continue;

			for( i=0; i<node_v->num; i++)
			{
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

			cdist=calc_dist_loc(loc_v,q->loc_v);
			if(!is_relevant_obj( obj_v, q_new))continue;
			
						
			for(vorc::Voronoi * v=v_head->next;v->next!=NULL;v=v->next){
				if(!v->queue.empty()){
				while(v->queue.top()->point->r<=cdist){
					VEventc * ev;
					ev=v->queue.top();
					v->queue.pop();
					v->RemoveParabola(ev,q,&costV);
					delete ev;
					if(v->queue.empty())break;
				}
			}
			}
			
			if(cdist>costV)break;
			//Construct a new query instance.

			for(vorc::Voronoi * v=v_head->next;v->next!=NULL;v=v->next){
				for(int i=0;i<3;i++){
					if(has_key_obj(obj_v,v->key[i])){
						v->InsertParabola(new VPoint(loc_v->coord[0],loc_v->coord[1],q->loc_v->coord[0],q->loc_v->coord
							[1]),q,&costV);
						break;
					}
				}
			}

			
			/*s*/
			stat_v.n_k ++;
			/*s*/
		}//else
	}//while
	for(vorc::Voronoi * v=v_head;v->next!=NULL;){
		vorc::Voronoi * v_tmp=v;
		v=v->next;
		delete v_tmp;
	}
	release_b_heap( U);
	release_query(q_new); 
	return costV;
}


double Approx_SSQ2( query_ssq* q)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, costV, dist;
	obj_set_t* S_a;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_t* o, *o_next;
	loc_t* loc_v;
	bst_t* obj_pair;
	b_heap_t* R_heap;
	query_t* q_new;
	q_new=alloc_query();
	q_new->loc_v=copy_loc(q->loc_v);
	q_new->psi_v=alloc_psi();
	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(q_new->psi_v,q->target);
	//Compute the LB and UB of cost(S*, q).c_1;
	S_a=Cao_Appro1(q_new);
	LB=0.0;
	UB=calc_dist_loc(get_obj_loc(const_NN_key(q->loc_v ,q->target,NULL)),q->loc_v)+comp_cost(S_a,q_new);
	cost_c = UB;
	costV=UB;
	if( S_a == NULL)
		return NULL;

	//Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);
	
	region_u = range_query( disk_u, q_new);
	region_l = range_query( disk_l, q_new);
	refine_region( region_l, disk_l);

	obj_exclusion_disk( region_u, disk_l);

	R = region_u;

	vorc::Voronoi * v_head=new vorc::Voronoi();
	vorc::Voronoi * v_last=v_head;
	for(k_node_t * psi1=q->psi_v->k_head->next;psi1->next!=NULL;psi1=psi1->next){
		for(k_node_t * psi2=psi1->next;psi2->next!=NULL;psi2=psi2->next){
		//	if(psi2==NULL)break;
			for(k_node_t * psi3=psi2->next;;psi3=psi3->next){
				if(psi3==NULL)break;
				vorc::Voronoi * v=new vorc::Voronoi();
				v->key[0]=psi1->key;
				v->key[1]=psi2->key;
				v->key[2]=psi3->key;
				v_last->next=v;
				v_last=v_last->next;
				v->root=NULL;
		}
	}
}
	for(k_node_t * psi1=q->psi_v->k_head->next;psi1->next!=NULL;psi1=psi1->next){
		for(k_node_t * psi2=psi1->next;;psi2=psi2->next){
			if(psi2==NULL)break;
			vorc::Voronoi * v=new vorc::Voronoi();
			v->key[0]=q->target;
			v->key[1]=psi1->key;
			v->key[2]=psi2->key;
			v_last->next=v;
			v_last=v_last->next;
			v->root=NULL;
		}
	}
	v_last->next=NULL;
	//Sort the objects in R by their distances to q.
	R_heap = heap_sort_obj_set( R, q_new);
	top = b_h_get_top( R_heap);

	o_next = R_heap->obj_arr[ top].obj_v;

	while( true)
	{
		if( o_next == NULL)
			break;

		o = o_next;

		top = b_h_get_top( R_heap);
		if( top == 0)
			o_next = NULL;
		else
			o_next = R_heap->obj_arr[ top].obj_v;

		loc_v = get_obj_loc( o);
		dist = calc_dist_loc( loc_v, q->loc_v);
			if(!is_relevant_obj( o, q_new))continue;
			
						
			for(vorc::Voronoi * v=v_head->next;v->next!=NULL;v=v->next){
				if(!v->queue.empty()){
				while(v->queue.top()->point->r<=dist){
					VEventc * ev;
					ev=v->queue.top();
					v->queue.pop();
					v->RemoveParabola(ev,q,&costV);
					delete ev;
					if(v->queue.empty())break;
				}
			}
			}
			
			if(dist>costV)break;
			//Construct a new query instance.

			for(vorc::Voronoi * v=v_head->next;v->next!=NULL;v=v->next){
				for(int i=0;i<3;i++){
					if(has_key_obj(o,v->key[i])){
						v->InsertParabola(new VPoint(loc_v->coord[0],loc_v->coord[1],q->loc_v->coord[0],q->loc_v->coord
							[1]),q,&costV);
						break;
					}
				}
			}

			stat_v.n_k ++;
			/*s*/
	}//while
	for(vorc::Voronoi * v=v_head;v->next!=NULL;){
		vorc::Voronoi * v_tmp=v;
		v=v->next;
		delete v_tmp;
	}

	release_query(q_new); 
	release_loc( loc_v);

	//Release the resource.
	release_disk( disk_u);
	release_disk( disk_l);
	release_obj_set( R);
	release_obj_set( region_l);
	release_b_heap( R_heap);

	return costV;
}


double Adapt2_Cao_Appro1( query_ssq* q)
{
	obj_set_t* S;
	obj_t* obj_v;
	k_node_t* k_node_v;
	double costc;

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
	obj_t* obj_new=const_NN_key( q->loc_v, q->target, NULL);
	add_obj_set_entry(obj_new,S);
	query_t* qnew;
	qnew=alloc_query();
	qnew->loc_v=get_obj_loc(obj_new);
	qnew->psi_v=alloc_psi();
	copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(qnew->psi_v,q->target);
	costc=calc_dist_loc(get_obj_loc(obj_new),q->loc_v)+comp_cost(S,qnew);
	release_query(qnew);
	return costc;
}

double Adapt2_Cao_Appro2( query_ssq* q)
{
	int i, rear, top;
	KEY_TYPE c_key;
	B_KEY_TYPE costV, costV_1, dist;
	b_heap_t* U;
	obj_set_t* V;
	void* e;
	node_t* node_v;
	obj_t* obj_v;
	bst_node_t* bst_node_v;
	BIT_TYPE p_list;
	loc_t* loc_v;
	query_t* q_new;

	//V_1=alloc_obj_set( );

	//printf( "Cao-Appro2:\n");

	U = alloc_b_heap( INI_HEAP_SIZE);

	rear = 1;
	U->obj_arr[ rear].element = ( void*)IRTree_v.root;
	U->obj_arr[ rear].e_tag = 1;
	U->obj_arr[ rear].key = calc_minDist_node( IRTree_v.root, q->loc_v);

	b_h_insert( U, rear++);
	
	costV = Adapt2_Cao_Appro1( q);	

	/*s*/
	stat_v.n_k ++;
	/*s*/
	q_new=alloc_query();
	q_new->loc_v=q->loc_v;
	q_new->psi_v=alloc_psi();
	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(q_new->psi_v,q->target);

	V=Cao_Appro1(q_new);
	//Find the keyword only contained by the object that is farthest from @q in V.
	c_key = find_critical_keyword( V, q_new);
//	release_query(q_new);
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
			//q_new = alloc_query( );
			query_t* qnew;
			qnew = ( query_t*)malloc( sizeof( query_t));
			memset( qnew, 0, sizeof( query_ssq));
			qnew->loc_v=alloc_loc(2);
			qnew->loc_v = loc_v;
		//	qnew->target=q->target;
			qnew->psi_v = alloc_psi( );
			copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);
			obj_set_t * S=Cao_Appro1(qnew);
			obj_t* obj_new=const_NN_key( loc_v, q->target, NULL);
			add_obj_set_entry(obj_new,S);
			//Solve the new query.
			add_psi_entry(qnew->psi_v,q->target);
			costV_1 = calc_dist_loc(q->loc_v,get_obj_loc(obj_new))+comp_cost(S,qnew);
			//add_obj_set_entry( obj_v, V_1);

		//	costV_1 = comp_cost( V_1, q);

			if( costV_1 < costV)
			{
				costV = costV_1;
		//		V_tmp = V;
		//		V = V_1;
				
			//	release_obj_set( V_tmp);
			}
			else stat_v.n_k ++;
			/*s*/
		}//else
	}//while
	//release_obj_set( V_1);
	release_b_heap( U);
//	release_query( q_new);	
	return costV;
}




double Adapt2_CostEnum_Approx( query_ssq* q, int prune_tag)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, cost, dist;
	obj_set_t* S_a;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_t* o, *o_next;
	loc_t* loc_v;
	bst_t* obj_pair;
	b_heap_t* R_heap;
	query_t* qnew=alloc_query( );
	qnew->loc_v=q->loc_v;
	qnew->psi_v=alloc_psi();
	copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(qnew->psi_v,q->target);
	//Compute the LB and UB of cost(S*, q).c_1;
	S=Cao_Appro1(qnew);
	obj_t* obj_new=const_NN_key( q->loc_v, q->target, NULL);
	S_a = comp_bounds( qnew, LB, UB);
	qnew->loc_v=get_obj_loc(obj_new);
	UB=calc_dist_loc(q->loc_v,get_obj_loc(obj_new))+comp_cost(S,qnew);
	cost_c = UB;
	if( S_a == NULL)
		return NULL;
	qnew->loc_v=q->loc_v;
	//Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);
	
	region_u = range_query( disk_u, qnew);
	region_l = range_query( disk_l, qnew);
	refine_region( region_l, disk_l);

	/*t/
	print_obj_set( region_l, stdout);
	/*t*/

	obj_exclusion_disk( region_u, disk_l);

	R = region_u;

	/*t/
	print_obj_set( R, stdout);
	/*t*/

	/*s/
	printf( "#cands: %i\n", R->obj_n);
	/*s*/

	//Pre-checking.
	if( R->obj_n == 0)
		goto E;
	loc_v=alloc_loc(2);
	//Sort the objects in R by their distances to q.
	R_heap = heap_sort_obj_set( R, qnew);
	top = b_h_get_top( R_heap);

	o_next = R_heap->obj_arr[ top].obj_v;

	//Search.
	while( true)
	{
		if( o_next == NULL)
			break;

		o = o_next;

		top = b_h_get_top( R_heap);
		if( top == 0)
			o_next = NULL;
		else
			o_next = R_heap->obj_arr[ top].obj_v;

		loc_v = get_obj_loc( o);
		dist = calc_dist_loc( loc_v, q->loc_v);


		if( dist > cost_c)
			break;

		/*s*/
		stat_v.n_1_sum ++;
		/*s*/
		disk_t* disk_tmp;
		disk_tmp = const_disk( q->loc_v, dist);
			S = ConstructFeasibleSet_Appro( o, qnew);
			obj_new=const_NN_key( loc_v, q->target, disk_tmp);
			/*t/
			printf( "base object: %i\n", o->id);
			print_obj_set( S, stdout);
			/*t*/
		

		if( !S)
			continue;
		qnew->loc_v=get_obj_loc(obj_new);
		cost = calc_dist_loc(get_obj_loc(obj_new),q->loc_v)+comp_cost( S, qnew);
		qnew->loc_v=q->loc_v;
		if( cost < cost_c)
		{
			release_obj_set( S_a);
			S_a = S;
			cost_c = cost;
		}
		else
			release_obj_set( S);
	}//while
E:
			release_loc( loc_v);
	remove_identical_obj( S_a);

	//Release the resource.
	release_disk( disk_u);
	release_disk( disk_l);
	release_obj_set( R);
	release_obj_set( region_l);
	release_b_heap( R_heap);

	return cost_c;
}
double SSQ_A2(query_ssq* q){
	double costV=Adapt2_Cao_Appro2(q);
	disk_t* disk=alloc_disk(2);
	set_disk(disk,q->loc_v,costV);
	obj_set_t* objset=alloc_obj_set();
	query_t* qnew=alloc_query();
	qnew->loc_v=alloc_loc(2);
	qnew->loc_v=copy_loc(q->loc_v);
	qnew->psi_v=alloc_psi();
	copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(qnew->psi_v,q->target);
	objset=range_query(disk, qnew);
	b_heap_t* bheap=heap_sort_obj_set_y( objset);
	int top;
	double dist;
	obj_t* o_next;
	
	Tree tree1;
	tree1.set_show(0);
	point p;
/*		while( true)
	{

		top = b_h_get_top( bheap);
		if( top == 0)
			o_next = NULL;
		else
			o_next = bheap->obj_arr[ top].obj_v;
		if( o_next == NULL)
			break;
		p.x_coord=get_obj_loc(o_next)->coord[0];
		p.y_coord=get_obj_loc(o_next)->coord[1];
		tree1.addNode(p);
		}*/
	release_query(qnew);
	release_disk(disk);
	return costV;
}


double SKECplus( query_ssq* q, float eps)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, cost, dist;
	obj_set_t* S_a, *V_1;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_node_t* o, *o_next;
	obj_t* o_second, *obj_v;
	loc_t* loc_v;
	bst_t* obj_pair;
	b_heap_t* R_heap, A_heap;
	query_t* qnew=alloc_query( );
	qnew->loc_v=copy_loc(q->loc_v);
	qnew->psi_v=alloc_psi();
	copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(qnew->psi_v,q->target);
	//Compute the LB and UB of cost(S*, q).c_1;
	UB=Adapt1(q,9,0);
	LB=0.5*UB;
	cost_c = UB;
	cost=cost_c;
	//Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);
	
	region_u = range_query( disk_u, qnew);
	region_l = range_query( disk_l, qnew);
	refine_region( region_l, disk_l);
	obj_exclusion_disk( region_u, disk_l);

	R = region_u;
	//Pre-checking.
	if( R->obj_n == 0)
		goto E;
	loc_v=alloc_loc(2);
	//Sort the objects in R by their distances to q.

	o_next = R->head->next;

	//Search.
	while( UB-LB>eps*LB)
	{
		while(true){
			if( o_next == NULL)
				break;
			o = o_next;

			loc_v = get_obj_loc( o->obj_v);
			dist = calc_dist_loc( loc_v, q->loc_v);
			set_disk( disk_u, loc_v, UB);
			qnew->loc_v=get_obj_loc(o->obj_v);
			region_u = range_query( disk_u, qnew);

			obj_set_t* obj_set_c=alloc_obj_set();			

			obj_node* obj_c=region_u->head->next;
			while(obj_c!=NULL){
				
				float xm,ym,x1,y1,x2,y2,q,x,y;
				x1=loc_v->coord[0];
				y1=obj_c->obj_v->MBR[0].min;
				x2=loc_v->coord[1];
				y2=obj_c->obj_v->MBR[1].min;
				q=sqrt(pow(x2-x1,(float) 2.0) + pow(y2-y1,(float) 2.0));
				xm=loc_v->coord[0]+obj_c->obj_v->MBR[0].min;
				xm=xm/2;
				ym=loc_v->coord[1]+obj_c->obj_v->MBR[1].min;
				ym=ym/2;
				x = xm + sqrt(pow(float (UB/2),(float) 2.0)-pow(q/2, (float) 2.0))*(y1-y2)/q;
				y = ym + sqrt(pow(float (UB/2),(float) 2.0)-pow(q/2, (float) 2.0))*(x2-x1)/q;

				obj_t *obj_v1 = ( obj_t*)malloc( sizeof(obj_v));
				memset( obj_v1, 0, sizeof( obj_v));

				obj_v1->MBR=(range*)malloc(sizeof(range));
				memset(obj_v1->MBR, 0, sizeof(range));

				obj_v1->MBR[0].min=x;
				obj_v1->MBR[1].min=y;

				x = xm - sqrt(pow(float (UB/2),(float) 2.0)-pow(q/2, (float) 2.0))*(y1-y2)/q;
				y = ym - sqrt(pow(float (UB/2),(float) 2.0)-pow(q/2, (float) 2.0))*(x2-x1)/q;

				obj_t *obj_v2 = ( obj_t*)malloc( sizeof(obj_v));
				memset( obj_v2, 0, sizeof( obj_v));

				obj_v2->MBR=(range*)malloc(sizeof(range));
				memset(obj_v2->MBR, 0, sizeof(range));

				obj_v2->MBR[0].min=x;
				obj_v2->MBR[1].min=y;

				add_obj_set_entry(obj_v1,obj_set_c);
				add_obj_set_entry(obj_v2,obj_set_c);
				remove_obj_set_entry(region_u);
				obj_c=region_u->head->next;
			} 

			R_heap = heap_sort_obj_set( obj_set_c, qnew);
			top = b_h_get_top( R_heap);
			if(top==0)o_second=NULL;
			else o_second = R_heap->obj_arr[ top].obj_v;
			while(o_second!=NULL){
				query_t * q_new = alloc_query( );
				q_new->loc_v = copy_loc(get_obj_loc(o_second));
				q_new->psi_v = alloc_psi( );
				copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
				V_1 = Cao_Appro1( q_new);
				disk_t * disk_v=alloc_disk(2);

				set_disk(disk_v,q_new->loc_v,0.5*UB);

				obj_v = const_NN_key( q->loc_v, q->target, disk_v);

				if(obj_v!=NULL){
					q_new->loc_v=get_obj_loc(obj_v);
					
					cost = comp_cost( V_1, q_new);

					cost=cost+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

					if(cost<cost_c)cost_c=cost;

				}
				top = b_h_get_top( R_heap);
				if(top==0)o_second=NULL;
				else o_second = R_heap->obj_arr[ top].obj_v;
				release_disk( disk_v);
				release_query(q_new);
			}
			/*s*/
			stat_v.n_1_sum ++;
			release_b_heap( R_heap);
			release_obj_set(obj_set_c);
			o_next=o_next->next;
			/*s*/
		}
		if(cost>(LB+UB)/2)LB=(UB+LB)/2;
		else UB=(UB+LB)/2;
		o_next = R->head->next;
	}//while
E:
//			
//	remove_identical_obj( S_a);
	//Release the resource.
	//if(V_1!=NULL)release_obj_set( V_1);
	
	release_disk( disk_u);
	release_disk( disk_l);
	release_obj_set( R);
	release_obj_set( region_l);
//	release_loc( loc_v);
	printf("UB: %f cost: %f\n", UB, cost);
	return cost_c;
}

double SKEC( query_ssq* q)
{
	int top;
	B_KEY_TYPE LB, UB, cost_c, cost, dist;
	obj_set_t* S_a, *V_1;				//the current best solution.
	obj_set_t* S;				//the newly constructed feasible set.
	obj_set_t* R;				//region R.
	disk_t* disk_u, *disk_l;	//the outer and inner disks.
	obj_set_t* region_u, *region_l;
	obj_node_t* o, *o_next;
	obj_t* o_second, *obj_v;
	loc_t* loc_v;
	bst_t* obj_pair;
	b_heap_t* R_heap, A_heap;
	query_t* qnew=alloc_query( );
	qnew->loc_v=copy_loc(q->loc_v);
	qnew->psi_v=alloc_psi();
	copy_k_list( qnew->psi_v->k_head, q->psi_v->k_head);
	add_psi_entry(qnew->psi_v,q->target);
	//Compute the LB and UB of cost(S*, q).c_1;
	UB=Adapt1(q,9,0);
	LB=0.5*UB;
	cost_c = UB;
	cost=cost_c;
	//Initialize region R.
	//An alternative implementation is possible here.
	//Direct range query with the range be a "ring".
	disk_u = alloc_disk( IRTree_v.dim);
	disk_l = alloc_disk( IRTree_v.dim);

	set_disk( disk_u, q->loc_v, UB);
	set_disk( disk_l, q->loc_v, LB);
	
	region_u = range_query( disk_u, qnew);
	region_l = range_query( disk_l, qnew);
	refine_region( region_l, disk_l);
	obj_exclusion_disk( region_u, disk_l);

	R = region_u;
	//Pre-checking.
	if( R->obj_n == 0)
		goto E;
	loc_v=alloc_loc(2);
	//Sort the objects in R by their distances to q.

	o_next = R->head->next;

		for(obj_node_t * p1=region_u->head->next;p1->next!=NULL;p1=p1->next){
		for(obj_node_t * p2=p1->next;p2->next!=NULL;p2=p2->next){
		//	if(psi2==NULL)break;
			for(obj_node_t * p3=p2->next;;p3=p3->next){
				if(p3==NULL)break;
				float x1,y1,x2,y2,x3,y3,x,y,r;
				x1=p1->obj_v->MBR[0].min;
				y1=p1->obj_v->MBR[1].min;
				x2=p2->obj_v->MBR[0].min;
				y2=p2->obj_v->MBR[1].min;
				x3=p3->obj_v->MBR[0].min;
				y3=p3->obj_v->MBR[1].min;
				x=((pow(x1,(float)2.0)+pow(y1,(float)2.0))*(y2-y3) + (pow(x2,(float)2.0)+pow(y2,(float)2.0))*(y3-y1) + (pow(x3,(float)2.0)+pow(y3,(float)2.0))*(y1-y2))/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
				y=((pow(x1,(float)2.0)+pow(y1,(float)2.0))*(x3-x2) + (pow(x2,(float)2.0)+pow(y2,(float)2.0))*(x1-x3) + (pow(x3,(float)2.0)+pow(y3,(float)2.0))*(x2-x1))/(2*(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
				r=sqrt(pow(x,(float)2.0)+(pow(y,(float)2.0)+(pow(x1,(float)2.0)+pow(y1,(float)2.0))*(x2*y3-x3*y2) + (pow(x2,(float) 2.0)+pow(y2,(float)2.0))*(x3*y1-x1*y3) + (pow(x3,(float)2.0)+pow(y3,(float)2.0)*(x1*y2-x2*y1)))/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));

				query_t * q_new = alloc_query( );
				q_new->loc_v=alloc_loc(2);
				q_new->loc_v->coord[0] = x;
				q_new->loc_v->coord[1] = y;
				q_new->psi_v = alloc_psi( );
				copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
				V_1 = Cao_Appro1( q_new);
                if(V_1==NULL)continue;
				disk_t * disk_v=alloc_disk(2);
				set_disk(disk_v,q_new->loc_v,r);

				obj_v = const_NN_key( q->loc_v, q->target, disk_v);

				if(obj_v!=NULL){
					q_new->loc_v=get_obj_loc(obj_v);
					
					cost = comp_cost( V_1, q_new);

					cost=cost+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

					if(cost<cost_c)cost_c=cost;

					}
				release_disk( disk_v);
				release_query(q_new);
				}
	}
}
	for(obj_node_t * p1=region_u->head->next;p1->next!=NULL;p1=p1->next){
		for(obj_node_t * p2=p1->next;;p2=p2->next){
			if(p2==NULL)break;
			float x1,y1,x2,y2,x,y,r;
				x1=p1->obj_v->MBR[0].min;
				y1=p1->obj_v->MBR[1].min;
				x2=p2->obj_v->MBR[0].min;
				y2=p2->obj_v->MBR[1].min;
				
				x=(x1+x2)/2;
				y=(y1+y2)/2;
				r=sqrt(pow(x-x1,(float)2.0)+pow(y-y1,(float)2.0));

				query_t * q_new = alloc_query( );
				q_new->loc_v=alloc_loc(2);
				q_new->loc_v->coord[0] = x;
				q_new->loc_v->coord[1] = y;
				q_new->psi_v = alloc_psi( );
				copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head);
				V_1 = Cao_Appro1( q_new);
                if(V_1==NULL)continue;
				disk_t * disk_v=alloc_disk(2);
				set_disk(disk_v,q_new->loc_v,r);

				obj_v = const_NN_key( q->loc_v, q->target, disk_v);

				if(obj_v!=NULL){
					q_new->loc_v=get_obj_loc(obj_v);
					
					cost = comp_cost( V_1, q_new);

					cost=cost+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

					if(cost<cost_c)cost_c=cost;

					}
				release_disk( disk_v);
				release_query(q_new);
		}
	}

			/*s*/
	stat_v.n_1_sum ++;
E:
//			
//	remove_identical_obj( S_a);
	//Release the resource.
	//if(V_1!=NULL)release_obj_set( V_1);
	
	release_disk( disk_u);
	release_disk( disk_l);
	release_obj_set( R);
	release_obj_set( region_l);
//	release_loc( loc_v);
	printf("UB: %f cost: %f\n", UB, cost);
	return cost_c;
}
/*
 *	Generate a query based on the @key_n, @MBR and @data_v information.
 *
 *	@key_n indicates the number of keywords in the query.
 *	@MBR indicates the range of query location.
 *	@data_v contains the frequency information of keywords.
 */
query_t* gen_query( int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, j, rand_tmp, valid_key_n;
	int* rand_v;
	loc_t* loc_v;
	psi_t* psi_v;
	query_t* q;
	bst_node_t* low_n, *high_n, *bst_node_v;

	q = alloc_query( );

	//Generate the query location.
	loc_v = alloc_loc( data_v->dim);
	for( i=0; i<data_v->dim; i++)
		loc_v->coord[ i] = rand_f( MBR[ i].min, MBR[ i].max);

	q->loc_v = loc_v;

	//Generate the query keywords.
	psi_v = alloc_psi( );
	
	//valid_key_n = bst_search_range( data_v->key_freq_v, freq.min, freq.max, low_n, high_n);
	valid_key_n = bst_search_percentile_range( data_v->key_freq_v, low, high, low_n, high_n);
	if( valid_key_n < key_n)
	{
		fprintf( stderr, "The number of qualified keywords, %i, is below the requirement.\n", key_n);
		exit( 0);
	}

	//
	rand_v = ( int*)malloc( key_n * sizeof( int));
	memset( rand_v, 0, sizeof( int));

	for( i=0; i<key_n; i++)
	{
		rand_tmp = rand_i( 0, key_n-1);
		while( is_old( rand_v, i, rand_tmp))
			rand_tmp = rand_i( 0, valid_key_n-1);

		rand_v[ i] = rand_tmp;
	
		bst_node_v = low_n;
		for( j=0; j<rand_v[ i]; j++)
			bst_node_v = bst_successor( bst_node_v);

		add_psi_entry( psi_v, bst_node_v->key_id);
	}

	q->psi_v = psi_v;

	free( rand_v);

	return q;
}


/*
 *	
 */
int get_max_key( data_t* data_v)
{
	int i, k_max;
	k_node_t* k_node_iter;

	k_max = 0;	
	for( i=0; i<data_v->obj_n; i++)
	{
		k_node_iter = data_v->obj_v[ i].k_head->next;
		while( k_node_iter)
		{
			if( k_node_iter->key > k_max)
				k_max = ( int)k_node_iter->key;

			k_node_iter = k_node_iter->next;
		}
	}	

	return k_max;
}

/*
 *
 */
void collect_key_freq( data_t* data_v, int* freq_v, int size)
{
	int i;
	k_node_t* k_node_iter;

	for( i=0; i<data_v->obj_n; i++)
	{
		k_node_iter = data_v->obj_v[ i].k_head->next;
		while( k_node_iter)
		{
			/*t*/
			if( k_node_iter->key >= size)
			{
				printf( "Keywod out of range!\n");
				k_node_iter = k_node_iter->next;

				continue;
			}
			/*t*/

			freq_v[ ( int)k_node_iter->key] ++;
			
			k_node_iter = k_node_iter->next;
		}
	}
}

/*
 *
 */
int compare_int (const void * a, const void * b)
{
	if( *( ( int*)a) - *( ( int*)b) > 0)
		return 1;
	else if( *( ( int*)a) - *( ( int*)b) < 0)
		return -1;
	else
		return 0;
}

/*
 *	
 */
int compare_key_freq (const void * a, const void * b)
{
	if( ( ( key_freq_pair_t*)a)->freq > ( ( key_freq_pair_t*)b)->freq)
		return 1;
	else if( ( ( key_freq_pair_t*)a)->freq < ( ( key_freq_pair_t*)b)->freq)
		return -1;
	else
		return 0;
}

/*
 *	An alternative implementation of "generating an query"
 */
query_t* gen_query2( int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, k_max, k_sta, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v;
	int* rand_v;
	loc_t* loc_v;
	psi_t* psi_v;
	query_t* q;

	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
		return NULL;
	
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);

	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));

	collect_key_freq( data_v, freq_v, k_max + 1);
	
	//Sort the keywords based on their frequency information.
	qsort( freq_v, k_max + 1, sizeof( int), compare_int);

	//Locate the keyword with the smallest freequency.
	for( i=1; i<=k_max; i++)
	{
		if( freq_v[ i] != 0)
			break;
	}

	k_sta = i;
	k_num = k_max - k_sta + 1;

	//Locate the keywords within the percentile range.
	p_low = low / 100.0;
	p_high = high /100.0;
	
	valid_k_sta = k_sta;
	while( ( valid_k_sta - k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_max)
		valid_k_sta ++;

	if( valid_k_sta == k_max)
	{
		printf( "Gen_query2 bug.\n");
		return NULL;
	}

	valid_k_end = valid_k_sta;
	while( ( valid_k_end - k_sta + 1) / ( float)k_num < p_high && valid_k_end < k_max)
		valid_k_end ++;

	valid_k_num = valid_k_end - valid_k_sta + 1;

	if( valid_k_num < key_n)
	{
		printf( "Valid keywords are not enough.\n");
		return NULL;
	}

	//Generate the query.
	q = alloc_query( );

	//query location.
	loc_v = alloc_loc( data_v->dim);
	for( i=0; i<data_v->dim; i++)
		loc_v->coord[ i] = rand_f( MBR[ i].min, MBR[ i].max);

	q->loc_v = loc_v;

	//query keywords.
	psi_v = alloc_psi( );

	//
	rand_v = ( int*)malloc( key_n * sizeof( int));
	memset( rand_v, 0, sizeof( int));

	for( i=0; i<key_n; i++)
	{
		rand_tmp = rand_i( 0, valid_k_num-1);
		while( is_old( rand_v, i, rand_tmp))
			rand_tmp = rand_i( 0, valid_k_num-1);

		rand_v[ i] = rand_tmp;

		add_psi_entry( psi_v, valid_k_sta + rand_v[ i]);
	}

	q->psi_v = psi_v;

	free( rand_v);
	free( freq_v);	

	return q;
}

/*
 *	Generate a set queries.
 */
query_t** gen_query_set( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i;
	query_t** q_set;

	q_set = ( query_t**)malloc( sizeof( query_t*) * query_n);
	memset( q_set, 0, sizeof( query_t*) * query_n);

	for( i=0; i<query_n; i++)
	{
		/*t/
		printf( "%i\n", i+1);
		/*t*/

		q_set[ i] = gen_query2( key_n, MBR, data_v, low, high);
	}
		
	return q_set;
}

/*
 *	Collect the number of keywords.
 */
int	get_key_num( int* freq_v, int size)
{
	int i, cnt;

	cnt = 0;
	for( i=0; i<size; i++)
	{
		if( freq_v[ i] != 0)
			cnt++;
	}

	return cnt;
}

/*
 *	The implementation of gen_query_set is problematic.
 *
 */
query_t** gen_query_set2( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, j, k, cnt, k_max, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v, *rand_v;
	loc_t* loc_v;
	psi_t* psi_v;
	query_t** q_set;
	query_ssq** q_set_ssq;
	key_freq_pair_t* key_freq_pair_v;

	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
		return NULL;
	
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);
	
	/*t/
	printf( "Maximum key id: %i\n", k_max);
	/*t*/

	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));

	collect_key_freq( data_v, freq_v, k_max + 1);

	//Restore the frequency information.
	//1. Collect the number of keys.
	k_num = get_key_num( freq_v, k_max + 1);

	if( k_num < key_n)
	{
		printf( "Insufficient keywords!\n");
		exit( 0);
	}

	//2. Restore.
	key_freq_pair_v = ( key_freq_pair_t*)malloc( k_num * sizeof( key_freq_pair_t));
	memset( key_freq_pair_v, 0, k_num * sizeof( key_freq_pair_t));

	cnt = 0;
	for( i=0; i<=k_max; i++)
	{
		if( freq_v[ i] != 0)
		{
			key_freq_pair_v[ cnt].key = i;
			key_freq_pair_v[ cnt++].freq = freq_v[ i];
		}
	}

	/*t*/
	if( cnt != k_num)
	{
		printf( "k_num inconsistency!\n");
		exit( 0);
	}
	/*t*/
	
	//3. Sort the keywords id based on the freq.
	qsort( key_freq_pair_v, k_num, sizeof( key_freq_pair_t), compare_key_freq);
	//for(i=0;i<=k_num;i++)printf("i:%d freq:%d ",i,freq_v[i]);
	//Sort the keywords based on their frequency information.
	//qsort( freq_v, k_max + 1, sizeof( int), compare_int);

	//Locate the keywords within the percentile range.
	p_low = ( 100 - high) / 100.0;
	p_high = ( 100 - low) / 100.0;
	
	valid_k_sta = 0;
	while( ( valid_k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_num)
		valid_k_sta ++;

	if( valid_k_sta == k_num)
	{
		printf( "Gen_query2 bug.\n");
		return NULL;
	}

	valid_k_end = valid_k_sta;
	while( ( valid_k_end + 1) / ( float)k_num < p_high && valid_k_end < k_num)
		valid_k_end ++;

	valid_k_num = valid_k_end - valid_k_sta + 1;
	if( valid_k_num < key_n)
	{
		printf( "Insufficient valid keywords!\n");
		exit( 0);
	}

	//Generate the queries.
	q_set = ( query_t**)malloc( sizeof( query_t*) * query_n);
	memset( q_set, 0, sizeof( query_t*) * query_n);

	for( i=0; i<query_n; i++)
	{
		q_set[ i] = alloc_query( );
		
		//query keywords.
		psi_v = alloc_psi( );
		
		//
		rand_v = ( int*)malloc( key_n * sizeof( int));
		memset( rand_v, 0, sizeof( int));
		
		for( k=0; k<key_n; k++)
		{
			rand_tmp = rand_i( 0, valid_k_num-1);
			while( is_old( rand_v, k, rand_tmp))
				rand_tmp = rand_i( 0, valid_k_num-1);
			
			rand_v[ k] = rand_tmp;
			
			add_psi_entry( psi_v, key_freq_pair_v[ valid_k_sta + rand_v[ k]].key);
		}
		
		q_set[ i]->psi_v = psi_v;

		


		//query location.
		loc_v = alloc_loc( data_v->dim);
		for( j=0; j<data_v->dim; j++)
			loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
		
		q_set[ i]->loc_v = loc_v;

		free( rand_v);
	}

	free( freq_v);
	free( key_freq_pair_v);

	return q_set;
}
query_ssq** gen_query_set2_ssq( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, j, k, cnt, k_max, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v, *rand_v;
	loc_t* loc_v;
	psi_t* psi_v;
	query_t** q_set;
	query_ssq** q_set_ssq;
	key_freq_pair_t* key_freq_pair_v;

	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
		return NULL;
	
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);
	
	/*t/
	printf( "Maximum key id: %i\n", k_max);
	/*t*/

	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));

	collect_key_freq( data_v, freq_v, k_max + 1);

	//Restore the frequency information.
	//1. Collect the number of keys.
	k_num = get_key_num( freq_v, k_max + 1);

	if( k_num < key_n)
	{
		printf( "Insufficient keywords!\n");
		exit( 0);
	}

	//2. Restore.
	key_freq_pair_v = ( key_freq_pair_t*)malloc( k_num * sizeof( key_freq_pair_t));
	memset( key_freq_pair_v, 0, k_num * sizeof( key_freq_pair_t));

	cnt = 0;
	for( i=0; i<=k_max; i++)
	{
		if( freq_v[ i] != 0)
		{
			key_freq_pair_v[ cnt].key = i;
			key_freq_pair_v[ cnt++].freq = freq_v[ i];
		}
	}

	/*t*/
	if( cnt != k_num)
	{
		printf( "k_num inconsistency!\n");
		exit( 0);
	}
	/*t*/

	//3. Sort the keywords id based on the freq.
	qsort( key_freq_pair_v, k_num, sizeof( key_freq_pair_t), compare_key_freq);
	//Sort the keywords based on their frequency information.
	//qsort( freq_v, k_max + 1, sizeof( int), compare_int);

	//Locate the keywords within the percentile range.
	p_low = ( 100 - high) / 100.0;
	p_high = ( 100 - low) / 100.0;
	
	valid_k_sta = 0;
	while( ( valid_k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_num)
		valid_k_sta ++;

	if( valid_k_sta == k_num)
	{
		printf( "Gen_query2 bug.\n");
		return NULL;
	}

	valid_k_end = valid_k_sta;
	while( ( valid_k_end + 1) / ( float)k_num < p_high && valid_k_end < k_num)
		valid_k_end ++;

	valid_k_num = valid_k_end - valid_k_sta + 1;
	if( valid_k_num < key_n)
	{
		printf( "Insufficient valid keywords!\n");
		exit( 0);
	}

	//Generate the queries.
	q_set = ( query_t**)malloc( sizeof( query_t*) * query_n);
	memset( q_set, 0, sizeof( query_t*) * query_n);
	q_set_ssq = ( query_ssq**)malloc( sizeof( query_ssq*) * query_n);
	memset( q_set_ssq, 0, sizeof( query_ssq*) * query_n);
	for( i=0; i<query_n; i++)
	{
		q_set[ i] = alloc_query( );
		
		//query keywords.
		psi_v = alloc_psi( );
		
		//
		rand_v = ( int*)malloc( key_n * sizeof( int));
		memset( rand_v, 0, sizeof( int));
		
		for( k=0; k<key_n; k++)
		{
			rand_tmp = rand_i( 0, valid_k_num-1);
			while( is_old( rand_v, k, rand_tmp))
				rand_tmp = rand_i( 0, valid_k_num-1);
			
			rand_v[ k] = rand_tmp;
			
			add_psi_entry( psi_v, key_freq_pair_v[ valid_k_sta + rand_v[ k]].key);
		}
		q_set_ssq[i] = ( query_ssq*)malloc( sizeof( query_ssq));
		memset( q_set_ssq[i], 0, sizeof( query_ssq));
		q_set_ssq[ i]->psi_v = psi_v;
		//while( is_old( rand_v, k, rand_tmp))
				rand_tmp = rand_i( 0, valid_k_num-1);

		q_set_ssq[i]->target=key_freq_pair_v[ valid_k_sta+rand_tmp].key;
		//query location.
		loc_v = alloc_loc( data_v->dim);
		for( j=0; j<data_v->dim; j++)
			loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
		
		q_set_ssq[ i]->loc_v = loc_v;

		
		free( rand_v);
	}

	free( freq_v);
	free( key_freq_pair_v);

	return q_set_ssq;
}
/*
 *	Generate a location based on two objects @obj_v1, @obj_v2.
 */
loc_t* gen_loc( obj_t* obj_v1, obj_t* obj_v2, KEY_TYPE key1, KEY_TYPE key2)
{
	int i;
	loc_t* loc_v, *loc_v1, *loc_v2;
	disk_t* disk_v;
	query_t* q;
	obj_set_t* obj_set_v;
	k_node_t*  k_node_v;

	loc_v1 = get_obj_loc( obj_v1);
	loc_v2 = get_obj_loc( obj_v2);
	
	//Check the feasibility.
	//range query.
	disk_v = alloc_disk( IRTree_v.dim);
	for( i=0; i<IRTree_v.dim; i++)
	{
		disk_v->center->coord[ i] = ( loc_v1->coord[ i] + loc_v2->coord[ i]) / 2;
	}
	disk_v->radius = calc_dist_loc( loc_v1, loc_v2) / 2;

	q = alloc_query( );
	q->loc_v = alloc_loc( IRTree_v.dim);
	q->psi_v = alloc_psi( );
	k_node_v = q->psi_v->k_head;
	add_keyword_entry( k_node_v, key1);
	add_keyword_entry( k_node_v, key2);

	obj_set_v = range_query( disk_v, q);

	if( obj_set_v->obj_n <= 2)
		loc_v = copy_loc( disk_v->center);
	else
		loc_v = NULL;

	release_disk( disk_v);
	release_query( q);
	release_obj_set( obj_set_v);
	release_loc( loc_v1);
	release_loc( loc_v2);
	
	return loc_v;
}

/*
 *	Retrieve all the objects that contain keyword @key.
 */
obj_set_t* retrieve_obj_key( data_t* data_v, KEY_TYPE key)
{
	int i;
	obj_set_t* obj_set_v;
	
	obj_set_v = alloc_obj_set( );
	
	for( i=0; i<data_v->obj_n; i++)
	{
		if( has_key_obj( data_v->obj_v + i, key))
			add_obj_set_entry( data_v->obj_v + i, obj_set_v);
	}

	return obj_set_v;
}

/*
 *	Generate a location based on two keywords @key1, @key2.
 *
 *	A certain of randomness is injected when selecting objects based on the keys.
 */
loc_t* gen_loc( KEY_TYPE key1, KEY_TYPE key2, data_t* data_v)
{
	int i, j, obj_n1, obj_n2, rand_v1, rand_v2, cnt1, cnt2;
	bool* tag1, *tag2;
	obj_node_t* obj_node_v1, *obj_node_v2;
	obj_t* obj_v1, *obj_v2;
	loc_t* loc_v, *loc_v1;
	obj_set_t* obj_set_v1, *obj_set_v2;

	//Retrieve all the objects that contain @key1.
	obj_set_v1 = retrieve_obj_key( data_v, key1);
	obj_n1 = obj_set_v1->obj_n;
	if( obj_n1 == 0)
	{
		release_obj_set( obj_set_v1);
		return NULL;
	}

	tag1 = ( bool*)malloc( obj_n1 * sizeof( bool));
	memset( tag1, 0, obj_n1 * sizeof( bool));

	cnt1 = 0;
	while( cnt1 < obj_n1)
	{
		//Select one object that contains @key randomly.
		rand_v1 = rand_i( 0, obj_n1-1);
		while( tag1[ rand_v1])
			rand_v1 = rand_i( 0, obj_n1-1);

		tag1[ rand_v1] = true;
		cnt1 ++;

		obj_node_v1 = obj_set_v1->head->next;
		for( i=0; i<rand_v1; i++)
			obj_node_v1 = obj_node_v1->next;

		obj_v1 = obj_node_v1->obj_v;

		if( has_key_obj( obj_v1, key2))
			continue;

		loc_v1 = get_obj_loc( obj_v1);

		//Retrieve the object that contains @key2.
		obj_set_v2 = const_k_NN_key( loc_v1, key2, NULL, K_NN_PARAMETER);

		if( obj_set_v2 == NULL)
			continue;
		if( obj_set_v2->obj_n == 0)
		{
			release_obj_set( obj_set_v2);
			continue;
		}

		obj_n2 = obj_set_v2->obj_n;
		tag2 = ( bool*)malloc( obj_n2 * sizeof( bool));
		memset( tag2, 0, obj_n2 * sizeof( bool));

		//Try all possible pairs of two objects.
		cnt2 = 0;
		while( cnt2 < obj_n2)
		{
			rand_v2 = rand_i( 0, obj_n2-1);
			while( tag2[ rand_v2])
				rand_v2 = rand_i( 0, obj_n2-1);
			
			tag2[ rand_v2] = true;
			cnt2 ++;

			obj_node_v2 = obj_set_v2->head->next;
			for( j=0; j<rand_v2; j++)
				obj_node_v2 = obj_node_v2->next;

			obj_v2 = obj_node_v2->obj_v;
			if( has_key_obj( obj_v2, key1))
				continue;
	
			if( ( loc_v = gen_loc( obj_v1, obj_v2, key1, key2)) != NULL)
				return loc_v;
		}//while( cnt2)

		free( tag2);
		release_obj_set( obj_set_v2);
		release_loc( loc_v1);
	}//while( cnt1)

	free( tag1);
	release_obj_set( obj_set_v1);

	return NULL;
}

/*
 *	Generate a location based on the keywords in query @q.
 *
 */
void gen_loc( query_t* q, data_t* data_v, range* MBR)
{
	int i;
	loc_t* loc_v;
	k_node_t* k_node_v1, *k_node_v2;

	//Try all combinations of two keywords in q.
	if( q->psi_v->key_n <= 0)
	{
		q->loc_v = alloc_loc( data_v->dim);
		return;
	}

	k_node_v1 = q->psi_v->k_head->next;
	while( k_node_v1->next != NULL)
	{
		k_node_v2 = k_node_v1->next;
		while( k_node_v2 != NULL)
		{
			loc_v = gen_loc( k_node_v1->key, k_node_v2->key, data_v);

			if( loc_v != NULL)
			{
				q->loc_v = loc_v;
				return;
			}

			k_node_v2 = k_node_v2->next;
		}

		k_node_v1 = k_node_v1->next;
	}

	/*note*/
	printf( "Query with Random location.\n");
	q->loc_v = alloc_loc( data_v->dim);
	for( i=0; i<data_v->dim; i++)
		q->loc_v->coord[ i] = rand_f( MBR[ i].min, MBR[ i].max);
	/*note*/

	return;
}

/*
 *	The implementation of gen_query_set is problematic.
 *
 */
query_t** gen_query_set3( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, k, cnt, k_max, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v, *rand_v;
	psi_t* psi_v;
	query_t** q_set;
	key_freq_pair_t* key_freq_pair_v;

	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
		return NULL;
	
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);
	
	/*t*/
	printf( "Maximum key id: %i\n", k_max);
	/*t*/

	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));

	collect_key_freq( data_v, freq_v, k_max + 1);

	//Restore the frequency information.
	//1. Collect the number of keys.
	k_num = get_key_num( freq_v, k_max + 1);

	if( k_num < key_n)
	{
		printf( "Insufficient keywords!\n");
		exit( 0);
	}

	//2. Restore.
	key_freq_pair_v = ( key_freq_pair_t*)malloc( k_num * sizeof( key_freq_pair_t));
	memset( key_freq_pair_v, 0, k_num * sizeof( key_freq_pair_t));

	cnt = 0;
	for( i=0; i<=k_max; i++)
	{
		if( freq_v[ i] != 0)
		{
			key_freq_pair_v[ cnt].key = i;
			key_freq_pair_v[ cnt++].freq = freq_v[ i];
		}
	}

	/*t*/
	if( cnt != k_num)
	{
		printf( "k_num inconsistency!\n");
		exit( 0);
	}
	/*t*/

	//3. Sort the keywords id based on the freq.
	qsort( key_freq_pair_v, k_num, sizeof( key_freq_pair_t), compare_key_freq);
	
	//Sort the keywords based on their frequency information.
	//qsort( freq_v, k_max + 1, sizeof( int), compare_int);

	//Locate the keywords within the percentile range.
	p_low = low / 100.0;
	p_high = high /100.0;
	
	valid_k_sta = 0;
	while( ( valid_k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_num)
		valid_k_sta ++;

	if( valid_k_sta == k_num)
	{
		printf( "Gen_query2 bug.\n");
		return NULL;
	}

	valid_k_end = valid_k_sta;
	while( ( valid_k_end + 1) / ( float)k_num < p_high && valid_k_end < k_num)
		valid_k_end ++;

	valid_k_num = valid_k_end - valid_k_sta + 1;
	if( valid_k_num < key_n)
	{
		printf( "Insufficient valid keywords!\n");
		exit( 0);
	}

	//Generate the queries.
	q_set = ( query_t**)malloc( sizeof( query_t*) * query_n);
	memset( q_set, 0, sizeof( query_t*) * query_n);

	for( i=0; i<query_n; i++)
	{
		q_set[ i] = alloc_query( );
		
		//query keywords.
		psi_v = alloc_psi( );
		
		//
		rand_v = ( int*)malloc( key_n * sizeof( int));
		memset( rand_v, 0, sizeof( int));
		
		for( k=0; k<key_n; k++)
		{
			rand_tmp = rand_i( 0, valid_k_num-1);
			while( is_old( rand_v, k, rand_tmp))
				rand_tmp = rand_i( 0, valid_k_num-1);
			
			rand_v[ k] = rand_tmp;
			
			add_psi_entry( psi_v, key_freq_pair_v[ valid_k_sta + rand_v[ k]].key);
		}
		
		q_set[ i]->psi_v = psi_v;


		//query location.
		/*
		loc_v = alloc_loc( data_v->dim);
		for( j=0; j<data_v->dim; j++)
			loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
		
		q_set[ i]->loc_v = loc_v;
		*/
		/*t/
		if( i == 29)
			printf( "");
		/*t*/
		gen_loc( q_set[i], data_v, MBR);

		free( rand_v);
	}

	free( freq_v);
	free( key_freq_pair_v);

	return q_set;
}
query_ssq** gen_query_set3_ssq( int query_n, int key_n, range* MBR, data_t* data_v, int low, int high)
{
	int i, k, cnt, k_max, k_num, rand_tmp, valid_k_sta, valid_k_end, valid_k_num;
	float p_low, p_high;
	int* freq_v, *rand_v;
	psi_t* psi_v;
	query_ssq** q_set;
	key_freq_pair_t* key_freq_pair_v;

	if( low > high || low < 0 || low > 100 || high < 0 || high > 100)
		return NULL;
	
	//Collect the "keyword range" information.
	k_max = get_max_key( data_v);
	
	/*t*/
	printf( "Maximum key id: %i\n", k_max);
	/*t*/

	//Collect the "frequency" information.
	freq_v = ( int*)malloc( ( k_max + 1) * sizeof( int));
	memset( freq_v, 0, ( k_max + 1) * sizeof( int));

	collect_key_freq( data_v, freq_v, k_max + 1);

	//Restore the frequency information.
	//1. Collect the number of keys.
	k_num = get_key_num( freq_v, k_max + 1);

	if( k_num < key_n)
	{
		printf( "Insufficient keywords!\n");
		exit( 0);
	}

	//2. Restore.
	key_freq_pair_v = ( key_freq_pair_t*)malloc( k_num * sizeof( key_freq_pair_t));
	memset( key_freq_pair_v, 0, k_num * sizeof( key_freq_pair_t));

	cnt = 0;
	for( i=0; i<=k_max; i++)
	{
		if( freq_v[ i] != 0)
		{
			key_freq_pair_v[ cnt].key = i;
			key_freq_pair_v[ cnt++].freq = freq_v[ i];
		}
	}

	/*t*/
	if( cnt != k_num)
	{
		printf( "k_num inconsistency!\n");
		exit( 0);
	}
	/*t*/

	//3. Sort the keywords id based on the freq.
	qsort( key_freq_pair_v, k_num, sizeof( key_freq_pair_t), compare_key_freq);
	
	//Sort the keywords based on their frequency information.
	//qsort( freq_v, k_max + 1, sizeof( int), compare_int);

	//Locate the keywords within the percentile range.
	p_low = low / 100.0;
	p_high = high /100.0;
	
	valid_k_sta = 0;
	while( ( valid_k_sta + 1) / ( float)k_num < p_low && valid_k_sta < k_num)
		valid_k_sta ++;

	if( valid_k_sta == k_num)
	{
		printf( "Gen_query2 bug.\n");
		return NULL;
	}

	valid_k_end = valid_k_sta;
	while( ( valid_k_end + 1) / ( float)k_num < p_high && valid_k_end < k_num)
		valid_k_end ++;

	valid_k_num = valid_k_end - valid_k_sta + 1;
	if( valid_k_num < key_n)
	{
		printf( "Insufficient valid keywords!\n");
		exit( 0);
	}

	//Generate the queries.
	q_set = ( query_ssq**)malloc( sizeof( query_ssq*) * query_n);
	memset( q_set, 0, sizeof( query_ssq*) * query_n);

	for( i=0; i<query_n; i++)
	{
		q_set[ i] = ( query_ssq*)malloc( sizeof( query_ssq));
		memset( q_set[i], 0, sizeof( query_ssq));
		//query keywords.
		psi_v = alloc_psi( );
		
		//
		rand_v = ( int*)malloc( key_n * sizeof( int));
		memset( rand_v, 0, sizeof( int));
		
		for( k=0; k<key_n; k++)
		{
			rand_tmp = rand_i( 0, valid_k_num-1);
			while( is_old( rand_v, k, rand_tmp))
				rand_tmp = rand_i( 0, valid_k_num-1);
			
			rand_v[ k] = rand_tmp;
			
			add_psi_entry( psi_v, key_freq_pair_v[ valid_k_sta + rand_v[ k]].key);
		}
		
		q_set[ i]->psi_v = psi_v;

		while( is_old( rand_v, k, rand_tmp))
		rand_tmp = rand_i( 0, valid_k_num-1);
		q_set[i]->target=key_freq_pair_v[ valid_k_end].key;
		//query location.
		/*
		loc_v = alloc_loc( data_v->dim);
		for( j=0; j<data_v->dim; j++)
			loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
		
		q_set[ i]->loc_v = loc_v;
		*/
		/*t/
		if( i == 29)
			printf( "");
		/*t*/
		loc_t * loc_v = alloc_loc( data_v->dim);
		for( int j=0; j<data_v->dim; j++)
			loc_v->coord[ j] = rand_f( MBR[ j].min, MBR[ j].max);
		
		q_set[ i]->loc_v = loc_v;

		free( rand_v);
	}

	free( freq_v);
	free( key_freq_pair_v);

	return q_set;
}
