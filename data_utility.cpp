/*
 *	Author: Victor Junqiu WEI
 *	Author: wjqjsnj@gmail.com
 */

#include "data_utility.h"
#include "bst.h"

#ifndef WIN32

/*
* GetCurTime is used to get the current running time in the current process.
*
* @Param curTime gets back the time information.
*
* @Return void.
*/
void GetCurTime( rusage* curTime)
{
	int ret = getrusage( RUSAGE_SELF, curTime);
	if( ret != 0)
	{
		fprintf( stderr, "The running time info couldn't be collected successfully.\n");
		//FreeData( 2);
		exit( 0);
	}
}

/*
* GetTime is used to get the 'float' format time from the start and end rusage structure.
* 
* @Param timeStart, timeEnd indicate the two time points.
* @Param userTime, sysTime get back the time information.
* 
* @Return void.
*/
void GetTime( struct rusage* timeStart, struct rusage* timeEnd, float* userTime, float* sysTime)
{
	(*userTime) = ((float)(timeEnd->ru_utime.tv_sec - timeStart->ru_utime.tv_sec)) + 
		((float)(timeEnd->ru_utime.tv_usec - timeStart->ru_utime.tv_usec)) * 1e-6;
	(*sysTime) = ((float)(timeEnd->ru_stime.tv_sec - timeStart->ru_stime.tv_sec)) +
		((float)(timeEnd->ru_stime.tv_usec - timeStart->ru_stime.tv_usec)) * 1e-6;
}

#endif


/*
*	IRTree_read_config reads the configuration info fro constructing the IRTree.
*/
IRTree_config_t* read_config_irtree( )
{
	//char des[MAX_DESCRIPTION_LENG];
	FILE* c_fp;

	IRTree_config_t* cfg = ( IRTree_config_t*)malloc( sizeof( IRTree_config_t));

	if( ( c_fp = fopen( CONFIG_FILE, "r")) == NULL)
	{
		fprintf( stderr, "The config file cannot be opened.\n");
		exit(0);
	}

	//reads the configuration info.
	fscanf( c_fp, "%s%s%s", cfg->loc_file, cfg->doc_file, cfg->tree_file);
	fscanf( c_fp, "%i%i%i", &cfg->obj_n, &cfg->key_n, &cfg->dim);		//data related.
	//fscanf( c_fp, "%i%i%i", &(cfg->M), &(cfg->m), &( cfg->split_opt));//R-tree related.
	fscanf( c_fp, "%i", &( cfg->split_opt));							//R-tree related.
	//fscanf( c_fp, "%i", &cfg->key_n);									//IF related.

	fclose( c_fp);

	return cfg;
}

/*
 *	Read the configuration for the CoSKQ problem.
 */
coskq_config_t* read_config_coskq( )
{
	coskq_config_t* cfg;
	FILE* c_fp;

	cfg = ( coskq_config_t*)malloc( sizeof( coskq_config_t));
	memset( cfg, 0, sizeof( coskq_config_t));

	if( ( c_fp = fopen( COSKQ_CONFIG_FILE, "r")) == NULL)
	{
		fprintf( stderr, "The coskq_config file cannot be opened.\n");
		exit(0);
	}

	//cost function option.
	fscanf( c_fp, "%i", &cfg->cost_measure);

	//algorithm option.
	fscanf( c_fp, "%i", &cfg->alg_opt);

	//data.
	fscanf( c_fp, "%i%i%s", &cfg->dim, &cfg->obj_n, cfg->loc_file);	//data related.

	fscanf( c_fp, "%i%s", &cfg->key_n, cfg->doc_file);

	fscanf( c_fp, "%i%s", &cfg->tree_tag, cfg->tree_file);
	
	//IR-tree.
	//fscanf( c_fp, "%i", &cfg->split_opt);

	//cost measurement.
	//fscanf( c_fp, "%i", &cfg->cost_measure);

	//CoSKQ.
	fscanf( c_fp, "%i%i", &cfg->q_key_n, &cfg->q_set_size);
	//fscanf( c_fp, "%f%f", &cfg->key_freq_range.min, &cfg->key_freq_range.max);
	fscanf( c_fp, "%i%i", &cfg->low, &cfg->high);
	//fscanf( c_fp, "%i",	&cfg->prune_opt);

	//fscanf( c_fp, "%f", &cfg->diag);

	//ratio.
	//fscanf( c_fp, "%f", &cfg->ratio_thr);

	fclose( c_fp);

	return cfg;
}

/*
 *	Add a key to the keyword list.
 */
void add_keyword_entry( k_node_t* &k_node_v, KEY_TYPE key)
{
	k_node_v->next = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( k_node_v->next, 0, sizeof( k_node_t));
	k_node_v->next->key = key;
	k_node_v = k_node_v->next;

	/*s*/
	stat_v.memory_v += sizeof( k_node_t);
	if( stat_v.memory_v > stat_v.memory_max)
		stat_v.memory_max = stat_v.memory_v;
	/*s*/
}

/*
 *	Copy the k_list info of @k_head2 to @k_head1.
 */
void copy_k_list( k_node_t* k_head1, k_node_t* k_head2)
{
	k_node_t* k_node_iter1, *k_node_iter2;

	k_node_iter1 = k_head1;
	k_node_iter2 = k_head2->next;
	while( k_node_iter2 != NULL)
	{
		add_keyword_entry( k_node_iter1, k_node_iter2->key);

		k_node_iter2 = k_node_iter2->next;
	}
}

/*
 *	Print the a list of keywords @k_head in @o_fp.
 */
void print_k_list( k_node_t* k_head, FILE* o_fp)
{
	k_node_t* k_node_iter;

	k_node_iter = k_head->next;
	while( k_node_iter != NULL)
	{
		fprintf( o_fp, "%.0lf  ", k_node_iter->key);

		k_node_iter = k_node_iter->next;
	}

	fprintf( o_fp, "\n");
}

/*
 *	Print the statistics maintained in stat_v.
 */
void print_coskq_stat( coskq_config_t* cfg, int cnt)
{
	FILE* s_fp;

	if( !( s_fp = fopen( COSKQ_STAT_FILE, "a")))
	{
		fprintf( stderr, "Cannot open the coskq_stat file.\n");
		exit( 0);
	}
    fprintf(s_fp,"keywords:%i algorithm option:%d\n",cfg->q_key_n,cfg->alg_opt);
	//average cost.
	fprintf( s_fp, "%lf\n\n", stat_v.aver_cost);

	//time.
	fprintf( s_fp, "%f\n%f\n\n", stat_v.irtree_build_time, stat_v.q_time);
	
	//memory.
	fprintf( s_fp, "%f\n", stat_v.memory_max / ( 1024 * 1024));

	//IR-tree memory.
	fprintf( s_fp, "%f\n", stat_v.tree_memory_max / ( 1024 * 1024));

	//data memory.
	//fprintf( s_fp, "%f\n\n", stat_v.data_memory_max / ( 1024 * 1024));
	
	//Other measurements.
	if( stat_v.achi_sum != 0)
	{
		stat_v.O_t_size_aver = stat_v.O_t_size_sum / stat_v.achi_sum;
		stat_v.psi_n_aver = stat_v.psi_n_sum / stat_v.achi_sum;
	}

	fprintf( s_fp, "%f\n%f\n%f\n%f\n\n", stat_v.n_1_sum / cnt, stat_v.achi_sum / cnt, 
		stat_v.O_t_size_aver, stat_v.psi_n_aver);

	fprintf( s_fp, "%f\n\n", stat_v.O_simp_size);

	//Cao_Appro2.
	if( cnt != 0)
		stat_v.n_k /= cnt;
	else
		stat_v.n_k = 0;
	//fprintf( s_fp, "%f\n\n", stat_v.n_k);

	//Ratio related.
	fprintf( s_fp, "%f\n%f\n%f\n%f\n\n", stat_v.ratio_min, stat_v.ratio_max, 
										stat_v.ratio_aver, stat_v.ratio_dev);

	fclose( s_fp);
}

/*
 *	Allocate the memory for an object.
 */
void alloc_obj( obj_t* obj_v, int dim)
{
	//obj_v = ( obj_t*)malloc( sizeof( obj_t));
	//memset( obj_v, 0, sizeof( obj_t));

	obj_v->k_head = ( k_node_t*)malloc( sizeof( k_node_t));
	memset( obj_v->k_head, 0, sizeof( k_node_t));

	obj_v->MBR = ( range*)malloc( dim * sizeof( range));
	memset( obj_v->MBR, 0, dim * sizeof( range));
}

/*
 *	Read the data based on the IRTree_config_t info.
 */
data_t*	read_data_irtree( IRTree_config_t* cfg)
{
	int i, j;
	KEY_TYPE key;
	char des;
	char keys[ TEXT_COL_MAX];
	char* tok;

	k_node_t* k_node_v;
	FILE* i_fp;

	data_t* data_v = ( data_t*)malloc( sizeof( data_t));
	memset( data_v, 0, sizeof( data_t));

	data_v->dim = cfg->dim;
	data_v->obj_n = cfg->obj_n;
	data_v->key_n = cfg->key_n;

	data_v->obj_v = ( obj_t*)malloc( sizeof( obj_t) * data_v->obj_n);
	memset( data_v->obj_v, 0, sizeof( obj_t) * data_v->obj_n);

	//data_v->key_freq_v = bst_ini( );

	//Read the loc info.
	if( ( i_fp = fopen( cfg->loc_file, "r")) == NULL)
	{
		fprintf( stderr, "Cannot open the loc file.\n");
		exit( 0);
	}

	for( i=0; i<data_v->obj_n; i++)
	{
		alloc_obj( data_v->obj_v + i, data_v->dim);

		fscanf( i_fp, "%i", &( data_v->obj_v[i ].id));

		for( j=0; j<data_v->dim; j++)
		{
			fscanf( i_fp, "%c%f", &des, &( data_v->obj_v[i].MBR[j].min));
			data_v->obj_v[i].MBR[j].max = data_v->obj_v[i].MBR[j].min;
		}
	}

	fclose( i_fp);

	//Read the keywords info.
	if( ( i_fp = fopen( cfg->doc_file, "r")) == NULL)
	{
		fprintf( stderr, "Cannot open the doc file.\n");
		exit( 0);
	}

	for( i=0; i<data_v->obj_n; i++)
	{
		k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
		memset( k_node_v, 0, sizeof( k_node_t));

		data_v->obj_v[ i].k_head = k_node_v;	

		fgets( keys, TEXT_COL_MAX, i_fp);

		tok = strtok( keys, " ,");
		while( ( tok = strtok( NULL, " ,")))
		{
			key = atoi( tok);

			add_keyword_entry( k_node_v, key);
/*
			//Update the frequency info of the keyword.
			bst_node_v = bst_search( data_v->key_freq_v, key);
			if( bst_node_v == NULL)
			{
				bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
				memset( bst_node_v, 0, sizeof( bst_node_t));

				bst_node_v->key = key;
				bst_node_v->freq = 1.0 / data_v->obj_n;
				
				bst_insert( data_v->key_freq_v, bst_node_v);
			}
			else
				bst_node_v->freq += 1.0 / data_v->obj_n;
*/
		}
	}

	fclose( i_fp);

	return data_v;
}

/*
 *
 */
data_t* alloc_data( int num)
{
	data_t* data_v;

	data_v = ( data_t*)malloc( sizeof( data_t));
	memset( data_v, 0, sizeof( data_t));

	data_v->obj_n = num;

	data_v->obj_v = ( obj_t*)malloc( sizeof( obj_t) * data_v->obj_n);
	memset( data_v->obj_v, 0, sizeof( obj_t) * data_v->obj_n);

	return data_v;
}

/*
 *	Read the data based on the coskq_config_t info.
 */
data_t*	read_data_coskq( coskq_config_t* cfg)
{
	int i, j;
	KEY_TYPE key;
	char des;
	char keys[ TEXT_COL_MAX];
	char* tok;
	data_t* data_v;

	k_node_t* k_node_v;
	FILE* i_fp;

	data_v = alloc_data( cfg->obj_n);
	
	data_v->dim = cfg->dim;
	data_v->key_n = cfg->key_n;
	//data_v->obj_n = cfg->obj_n;

	//data_v->key_freq_v = bst_ini( );

	//Read the loc info.
	if( ( i_fp = fopen( cfg->loc_file, "r")) == NULL)
	{
		fprintf( stderr, "Cannot open the loc file.\n");
		exit( 0);
	}

	for( i=0; i<data_v->obj_n; i++)
	{
		alloc_obj( data_v->obj_v + i, data_v->dim);

		fscanf( i_fp, "%i", &( data_v->obj_v[i ].id));

		for( j=0; j<data_v->dim; j++)
		{
			fscanf( i_fp, "%c%f", &des, &( data_v->obj_v[i].MBR[j].min));
			data_v->obj_v[i].MBR[j].max = data_v->obj_v[i].MBR[j].min;
		}
	}

	fclose( i_fp);

	//Read the keywords info.
	if( ( i_fp = fopen( cfg->doc_file, "r")) == NULL)
	{
		fprintf( stderr, "Cannot open the doc file.\n");
		exit( 0);
	}

	for( i=0; i<data_v->obj_n; i++)
	{
		k_node_v = ( k_node_t*)malloc( sizeof( k_node_t));
		memset( k_node_v, 0, sizeof( k_node_t));

		/*t/
		if( strcmp( cfg->loc_file, "web-loc") == 0)
		{
			if( ( i+1) % 100 == 0)
				printf( "%i\n", i / 100 + 1);
		}
		/*t*/

		data_v->obj_v[ i].k_head = k_node_v;	

		fgets( keys, TEXT_COL_MAX, i_fp);

		tok = strtok( keys, " ,");
		while( ( tok = strtok( NULL, " ,")))
		{
			key = atoi( tok);
			
			add_keyword_entry( k_node_v, key);
/*
			//Update the frequency info of the keyword.
			bst_node_v = bst_locate( data_v->key_freq_v, key);
			if( bst_node_v == NULL)
			{
				bst_node_v = ( bst_node_t*)malloc( sizeof( bst_node_t));
				memset( bst_node_v, 0, sizeof( bst_node_t));

				bst_node_v->key_id = key;
				bst_node_v->key = 1.0 / data_v->obj_n;
				
				bst_insert( data_v->key_freq_v, bst_node_v);
			}
			else
			{
				bst_node_v->key += 1.0 / data_v->obj_n;
				bst_update( data_v->key_freq_v, bst_node_v);
			}
/*
/*t/
			in_order_walk_non_recur( data_v->key_freq_v->root);
/*t*/
		}
	}

	fclose( i_fp);

	return data_v;
}

/*
 *	Release a list.
 */
void release_k_list( k_node_t* k_node_v)
{
	k_node_t* k_node_v1;
	while( k_node_v->next != NULL)
	{
		k_node_v1 = k_node_v;
		k_node_v = k_node_v->next;
		free( k_node_v1);
	}
	free( k_node_v);
}

/*
*	IRTree_free_data release the data read by 'IRTree_read_data'.
*/
void release_data( data_t* data_v)
{
	int i;

	for( i=0; i<data_v->obj_n; i++)
	{
		free( data_v->obj_v[ i].MBR);
		release_k_list( data_v->obj_v[ i].k_head);
	}

	free( data_v->obj_v);
	free( data_v);
}


/*
 *	Generates a number between [min, max] randomly.
 */
int rand_i( int min, int max)
{
	int i_rand;
	float ratio;

	if( min > max)
		return 0;

	if( max - min > RAND_MAX)
	{
		//printf( "rand_i: [%i, %i] out of range.\n", min, max);
		//exit( 0);

		ratio = rand_f( 0, 1);
		i_rand = min + ( int)( ratio * ( max - min));

		return i_rand;
	}

	//srand( time( NULL));

	if( min == max)
		i_rand = 0;
	else
		i_rand = rand( ) % ( max - min + 1);

	return i_rand + min;
}

/*
 *	Generate a random float value in [min_v, max_v].
 */
float rand_f( float min_v, float max_v)
{
	float rand_v;

	if( min_v > max_v)
		return 0;

	//srand( time(NULL));

	rand_v = float( rand( )) / RAND_MAX; //rand_v falls in [0,1] randomly.

	return ( max_v - min_v) * rand_v + min_v;
}


/*
*	Generate a float value that satisfy the normal/gaussian distribution N[mean, s_var^2].
*	
*	Polar form method, modified from the Internet.
*/
float gaussian_f( float mean, float s_var)
{
	float x1, x2, w, y1, y2;

	do
	{
		x1 = 2 * rand_f( 0, 1) - 1;
		x2 = 2 * rand_f( 0, 1) - 1;
		w = x1 * x1 + x2 * x2;
	} while( w >= 1);

	w = float(sqrt( float(-2.0 * log(w)) / w));
	y1 = x1 * w;					//y1 follows N[0, 1^2].

	y2 = y1 * s_var + mean;			//y2 follows N[mean, s_var^2].

	return y2;
}

/*
 *	Check whether an variable has been generated before.
 *
 */
bool is_old( int* rand_v, int cur, int var)
{
	int i;
	
	for( i=0; i<cur; i++)
	{
		if( rand_v[ i] == var)
			return true;
	}

	return false;
}

/*
 *
 */
range* collect_data_range( data_t* data_v)
{
	int i, j;
	range* MBR;

	MBR = ( range*)malloc( data_v->dim * sizeof( range));
	memset( MBR, 0, data_v->dim * sizeof( range));

	for( i = 0; i<data_v->dim; i++)
	{
		MBR[ i].min = FLT_MAX;
		MBR[ i].max = -FLT_MAX;

		for( j=0; j<data_v->obj_n; j++)
		{
			if( data_v->obj_v[ j].MBR[ i].min < MBR[ i].min)
				MBR[ i].min = data_v->obj_v[ j].MBR[ i].min;

			if( data_v->obj_v[ j].MBR[ i].max > MBR[ i].max)
				MBR[ i].max = data_v->obj_v[ j].MBR[ i].max;
		}
	}

	return MBR;
}

/*
 *
 */
void print_data( data_t* data_v, syn_config_t* s_cfg)
{
	int i, j;
	k_node_t* k_node_iter;

	FILE* loc_fp, *doc_fp;

	if( ( loc_fp = fopen( s_cfg->new_loc_file, "w")) == NULL ||
		( doc_fp = fopen( s_cfg->new_doc_file, "w")) == NULL)
	{
		printf( "Cannot open loc/doc files.\n");
		exit( 0);
	}

	//Print the loc info.
	for( i=0; i<data_v->obj_n; i++)
	{
		fprintf( loc_fp, "%i", data_v->obj_v[ i].id);
		for( j=0; j<data_v->dim; j++)
			fprintf( loc_fp, ",%f", data_v->obj_v[ i].MBR[ j].min);

		fprintf( loc_fp, "\n");
	}

	//Print the doc info.
	for( i=0; i<data_v->obj_n; i++)
	{
		fprintf( doc_fp, "%i", data_v->obj_v[ i].id);
		
		k_node_iter = data_v->obj_v[ i].k_head->next;
		while( k_node_iter)
		{
			fprintf( doc_fp, ",%i", ( int)k_node_iter->key);

			k_node_iter = k_node_iter->next;
		}
		fprintf( doc_fp, "\n");
	}

	fclose( loc_fp);
	fclose( doc_fp);
}

/*
 *
 */
void gen_syn_data( syn_config_t* s_cfg)
{
	int i, j, rand_v;
	float dev;
	data_t* base_data_v, *new_data_v;
	coskq_config_t* cfg;
	range* MBR;
	
	new_data_v = alloc_data( s_cfg->new_obj_n);
	new_data_v->dim = s_cfg->dim;
	new_data_v->key_n = s_cfg->key_n;
	//new_data_v->obj_n = s_cfg->new_obj_n;
	
	//Read the base data.
	cfg = ( coskq_config_t*)malloc( sizeof( coskq_config_t));
	memset( cfg, 0, sizeof( coskq_config_t));
	
	strcpy( cfg->loc_file, s_cfg->loc_file);
	strcpy( cfg->doc_file, s_cfg->doc_file);
	cfg->obj_n = s_cfg->obj_n;
	cfg->dim = s_cfg->dim;

	base_data_v = read_data_coskq( cfg);

	/*t/
	if( s_cfg->obj_n > RAND_MAX)
	{
		printf( "Insufficient random number range.\n");
		exit( 0);
	}
	/*t*/

	//Generate the synthetic data.
	MBR = collect_data_range( base_data_v);

	for( i=0; i<s_cfg->new_obj_n; i++)
	{
		alloc_obj( new_data_v->obj_v + i, new_data_v->dim);

		new_data_v->obj_v[ i].id = i + 1;

		//Generate the loc info.
		//Choose a base obj randomly.
		rand_v = rand_i( 0,  s_cfg->obj_n-1);

		//Use the base obj for generating the loc info.
		for( j=0; j<new_data_v->dim; j++)
		{
			dev = s_cfg->dev * ( MBR[ j].max - MBR[ j].min);
			new_data_v->obj_v[ i].MBR[ j].min = gaussian_f( base_data_v->obj_v[ rand_v].MBR[ j].min, dev);
			new_data_v->obj_v[ i].MBR[ j].max = new_data_v->obj_v[ i].MBR[ j].min;
		}
		printf("finish %dth object\n",i);
		//Generate the doc info.
		rand_v = rand_i( 0, s_cfg->obj_n-1);

		//Use the keywords info of the base obj.
		copy_k_list( new_data_v->obj_v[ i].k_head, base_data_v->obj_v[ rand_v].k_head);
	}

	//Print the synthetic data.
 	print_data( new_data_v, s_cfg);

	release_data( base_data_v);
	release_data( new_data_v);
	free( cfg);
	free( MBR);
}

/*
 *
 */
void batch_gen_syn_data( )
{
	int ins_n;
	syn_config_t s_cfg;
	FILE* c_fp;

	if( ( c_fp = fopen( SYN_CONFIG_FILE, "r")) == NULL)
	{
		printf( "Cannot open the syn_config file.\n");
		exit( 0);
	}

	ins_n = 1;
	while( 	fscanf( c_fp, "%s%s", s_cfg.loc_file, s_cfg.doc_file) != EOF)
	{
		printf( "#Instance %i ...\n", ins_n++);

		//read the config.
		fscanf( c_fp, "%i%i%i", &s_cfg.dim, &s_cfg.obj_n, &s_cfg.key_n);
		fscanf( c_fp, "%s%s", s_cfg.new_loc_file, s_cfg.new_doc_file);
		fscanf( c_fp, "%i%f", &s_cfg.new_obj_n, &s_cfg.dev);

		//generate the synthetic data.
		gen_syn_data_keyword( &s_cfg);
		printf("finish instance %i\n",ins_n);
	}

}



void gen_syn_data_keyword( syn_config_t* s_cfg)
{
	int i, j, rand_v;
	float dev;
	data_t* base_data_v, *new_data_v;
	coskq_config_t* cfg;
	range* MBR;
	
	new_data_v = alloc_data( s_cfg->obj_n);
	new_data_v->dim = s_cfg->dim;
	new_data_v->key_n = s_cfg->key_n;
	//new_data_v->obj_n = s_cfg->new_obj_n;
	
	//Read the base data.
	cfg = ( coskq_config_t*)malloc( sizeof( coskq_config_t));
	memset( cfg, 0, sizeof( coskq_config_t));
	
	strcpy( cfg->loc_file, s_cfg->loc_file);
	strcpy( cfg->doc_file, s_cfg->doc_file);
	cfg->obj_n = s_cfg->obj_n;
	cfg->dim = s_cfg->dim;

	base_data_v = read_data_coskq( cfg);

	/*t/
	if( s_cfg->obj_n > RAND_MAX)
	{
		printf( "Insufficient random number range.\n");
		exit( 0);
	}
	/*t*/

	//Generate the synthetic data.
	MBR = collect_data_range( base_data_v);

	for( i=0; i<s_cfg->obj_n; i++)
	{
		alloc_obj( new_data_v->obj_v + i, new_data_v->dim);

		new_data_v->obj_v[ i].id = i + 1;

		//Generate the loc info.
		//Choose a base obj randomly.
		rand_v = rand_i( 0,  s_cfg->obj_n-1);

		//Use the base obj for generating the loc info.
		for( j=0; j<new_data_v->dim; j++)
		{
			new_data_v->obj_v[ i].MBR[ j].min = base_data_v->obj_v[ rand_v].MBR[ j].min;
			new_data_v->obj_v[ i].MBR[ j].max = new_data_v->obj_v[ i].MBR[ j].min;
		}
		printf("finish %dth object\n",i);
		//Generate the doc info.
		copy_k_list( new_data_v->obj_v[ i].k_head, base_data_v->obj_v[ rand_v].k_head);
		//Use the keywords info of the base obj.
		for(j=0;j<s_cfg->new_obj_n-1;j++){
			rand_v = rand_i( 0,  s_cfg->obj_n-1);
			k_node_t* k_node_iter1, *k_node_iter2;

			k_node_iter1 = new_data_v->obj_v[ i].k_head;
			while(k_node_iter1->next!=NULL)k_node_iter1=k_node_iter1->next;
			k_node_iter2 = base_data_v->obj_v[ rand_v].k_head->next;
			while( k_node_iter2 != NULL)
			{
				add_keyword_entry( k_node_iter1, k_node_iter2->key);

				k_node_iter2 = k_node_iter2->next;
			}
		}
	}

	//Print the synthetic data.
 	print_data( new_data_v, s_cfg);

	release_data( base_data_v);
	release_data( new_data_v);
	free( cfg);
	free( MBR);
}
