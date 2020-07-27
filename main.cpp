// LSSK.cpp : Defines the entry point for the console application.
//



/*
 *	Author: Victor Junqiu WEI
 *	Author: wjqjsnj@gmail.com
 */
#include "irtree.h"
#include "costenum.h"
#include "cao_alg.h"
#include "string.h"


#include <float.h>
//#include<time.h>  


#define QUERY_SIZE		50
//#define	RATIO_THRESHOLD	1.3

#define QUERY_LOG		"query_log.txt"


IRTree_t		IRTree_v;
coskq_stat_t	stat_v;

int				cost_tag;



void nskq( int ratio_tag);


//void collect_ratio_stat( );

int main()
{
	//Checking the running environment.
	//printf( "RAND_MAX: %i\n", RAND_MAX);

	//printf( "%I64d\n", i);
	//printf( "INT_MAX: %i\n", INT_MAX);


//	batch_gen_syn_data( );

	//gen_irtree( );

	nskq(0);
	
	//read_tree( "hotel-tree");

	//print_and_check_tree( 1, "hotel-tree");

/*t*/
	//printf( "max unsigned long int: %u\n", ULONG_MAX);
	//printf( "max BIT_TYPE: %u\n", UINT_MAX);
/*t*/

/*t/
	int o_tag;

	//scanf( "%i", &o_tag);
	o_tag = 2;
	test_IRTree( o_tag);
/*t*/

/*Test the bit operators.*
	int	k;
	BIT_TYPE v;
	
	printf( "Input the BIT_TYPE:");
	scanf( "%u", &v);
	while( 1)
	{
		printf( "insert:\n");
		printf( "Input the location:");
		scanf( "%i", &k);

		insert_k_bit( v, k);
		printf( "%u\n\n", v);

		printf( "get:\n");
		printf( "Input the location:");
		scanf( "%i", &k);
		printf( "%i^th bit is %i\n\n", k, get_k_bit( v, k));

		printf( "delete:\n");
		printf( "Input the location:");
		scanf( "%i", &k);

		delete_k_bit( v, k);
		printf( "%u\n\n", v);		
	}
/*t*/	

	//test_bst( );

	return 0;
}


/*
 *	The interface for the NSKQ query.
 *
 */
void nskq( int ratio_tag)
{
	int i, j;
	B_KEY_TYPE cost, cost_opt, cost_sum, ratio_sum;
	B_KEY_TYPE* ratio;
	coskq_config_t* cfg;
	//coskq_stat_t* sta_v;
	data_t* data_v;
	query_ssq** q_set_ssq;
	query_t** q_set;
	range* MBR;
	obj_set_t* S, *S_opt;
	FILE *r_fp;
	memset( &stat_v, 0, sizeof( coskq_stat_t));
	
	//Read the cofig.
	cfg = read_config_coskq( );

	//
	cfg->prune_opt = 0;

	cost_tag = cfg->cost_measure;
	if( cost_tag == 1)
		printf( "The MaxSum measurement:\n");
	else
		printf( "The Diameter measurement:\n");

	//Read the data.
	printf( "Reading data ...\n");
	data_v = read_data_coskq( cfg);

#ifndef WIN32
	struct rusage IR_tree_sta, IR_tree_end;
	float sys_t, usr_t, usr_t_sum = 0;

	GetCurTime( &IR_tree_sta);
#endif

	//Option 1: Build the tree from scratch.
	//Build the IR-tree.
	if( cfg->tree_tag == 0)
	{
		printf( "Building IR-tree ...\n");
		build_IRTree( data_v);

		print_and_check_tree( 1, cfg->tree_file);
		//check_IF( );
	}
	else
	{
		//Option 2: Read an existing tree.
		printf( "Reading IR-Tree ...\n");
		read_tree( cfg->tree_file);
	}	

#ifndef WIN32
	GetCurTime( &IR_tree_end);
	GetTime( &IR_tree_sta, &IR_tree_end, &stat_v.irtree_build_time, &sys_t);
#endif

	//Get the whole range.
	MBR = get_MBR_node( IRTree_v.root, IRTree_v.dim);
    for(cfg->alg_opt=18;cfg->alg_opt<=18;cfg->alg_opt++){
        for(ratio_tag=0;ratio_tag<=1;ratio_tag++){
   //         if(ratio_tag==1){
     //           if(cfg->alg_opt<=6||cfg->alg_opt==8||cfg->alg_opt==11||cfg->alg_opt==12)break;
       //     }
            for(cfg->q_key_n=3;cfg->q_key_n<=3;cfg->q_key_n+=3){
	//Generate the set of querys.
	printf( "Generating queries ...\n");
	q_set = gen_query_set2( cfg->q_set_size, cfg->q_key_n, MBR, data_v, cfg->low, cfg->high);
	q_set_ssq = gen_query_set2_ssq( cfg->q_set_size, cfg->q_key_n, MBR, data_v, cfg->low, cfg->high);
	if( q_set_ssq == NULL)
	{
		printf( "Query generation failed!\n");
		exit( 0);
	}
	stat_v.ratio_min = 3;
	stat_v.ratio_max = 1;
	cost_opt=10000.0;
    usr_t_sum = 0;
	//Query.
	printf( "Performing Queries ...\n");
	if( cfg->alg_opt == 1)
		printf( "MaxSum(Dia)-Exact:\n");
	else if( cfg->alg_opt == 2)
		printf( "MaxSum(Dia)-Appro:\n");
	else if( cfg->alg_opt == 3)
		printf( "Cao-Exact:\n");
	else if( cfg->alg_opt == 4)
		printf( "Cao-Appro1:\n");
	else if(cfg->alg_opt == 5)
		printf( "Cao-Appro2:\n");
	else if(cfg->alg_opt==6)
		printf("MaxSum(Dia)-SSQ");
	else if(cfg->alg_opt==7)
		printf("NSKQ-Approx1.");
	else if(cfg->alg_opt==8)
		printf("SSQ-Exact");
	else if(cfg->alg_opt==9)
		printf("Adapt1_Cao1");
	else if(cfg->alg_opt==10)
		printf("Adapt1_Cao2");
	else if(cfg->alg_opt==11)
		printf("Adapt1_Exact_Cao");
	else if(cfg->alg_opt==12)
		printf("Adapt1_CoSKQ_Long_Exact");
	else if(cfg->alg_opt==13)
		printf("Adapt1_CoSKQ_Long_Appro");
	else if(cfg->alg_opt==14)
		printf("Adapt2_Cao1");
	else if(cfg->alg_opt==15)
		printf("Adapt2_Cao2");
	else if(cfg->alg_opt==16)
		printf("Adapt2_Long");
	else if(cfg->alg_opt==17)
		printf("NSKQ-Approx2");

	else if(cfg->alg_opt==18)
		printf("SKECplus");
	else if(cfg->alg_opt==19)
		printf("SKEC");
	ratio = ( double*)malloc( cfg->q_set_size * sizeof( double));
	memset( ratio, 0, cfg->q_set_size * sizeof( double));
	
	ratio_sum = 0;
	cost_sum = 0;
	for( i=0; i<cfg->q_set_size; i++)
	{
		printf( "Query #%i ...\n", i+1);

#ifndef WIN32
	struct rusage query_sta, query_end;

	GetCurTime( &query_sta);
#endif
		S_opt = NULL;
		if( cfg->alg_opt == 1)
		{
			//printf( "MaxSum(Dia)-Exact:\n");
			S = CostEnum_Exact( q_set[ i], cfg->prune_opt);
		}
		else if( cfg->alg_opt == 2)
		{
			//printf( "MaxSum(Dia)-Appro:\n");
			S = CostEnum_Appro( q_set[ i]);
			if( ratio_tag == 1)
				S_opt = Cao_Exact( q_set[ i]);
		}
		else if( cfg->alg_opt == 3)
		{
			//printf( "Cao-Exact:\n");
			S = Cao_Exact( q_set[ i]);
		}
		else if( cfg->alg_opt == 4)
		{
			//printf( "Cao-Appro1:\n");
			S = Cao_Appro1( q_set[ i]);
			if( ratio_tag == 1)
				S_opt = Cao_Exact( q_set[ i]);
		}
		else if(cfg->alg_opt == 5)
		{
			//printf( "Cao-Appro2:\n");
			S = Cao_Appro2( q_set[ i]);
			if( ratio_tag == 1)
				S_opt = Cao_Exact( q_set[ i]);
		}
		else if(cfg->alg_opt == 7)
		{
			cost=Approx_SSQ( q_set_ssq[i]);
			if( ratio_tag == 1)
				
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if(cfg->alg_opt==8){
	//		cost = SSQ_Exact( q_set_ssq[i], 1, cfg->prune_opt);
	//		cost_opt=cost;
		}
		else if( cfg->alg_opt==9 || cfg->alg_opt==10 || cfg->alg_opt==11 || cfg->alg_opt==12 || cfg->alg_opt==13){
			cost=Adapt1( q_set_ssq[ i],cfg->alg_opt,cfg->prune_opt);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==14){
			cost=Adapt2_Cao_Appro1( q_set_ssq[i]);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==15){
			cost=Adapt2_Cao_Appro2( q_set_ssq[i]);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==16){
			cost=Adapt2_CostEnum_Approx( q_set_ssq[i], cfg->prune_opt);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==17){
			cost=SSQ_A2( q_set_ssq[i]);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==18){
			cost=SKECplus( q_set_ssq[i],0.63);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
		else if( cfg->alg_opt==19){
			cost=SKEC( q_set_ssq[i]);
			if( ratio_tag == 1)
				cost_opt = Adapt1( q_set_ssq[ i],12,cfg->prune_opt);
		}
#ifndef WIN32 
	GetCurTime( &query_end);

	GetTime( &query_sta, &query_end, &usr_t, &sys_t);
	usr_t_sum += usr_t;
#endif
	if(cfg->alg_opt<=5){
		cost = comp_cost( S, q_set[ i]);
	}
		cost_sum += cost;
		
		if( S_opt != NULL)
		{
			//Ratio processing.
			if(cfg->alg_opt<=5)
			cost_opt = comp_cost( S_opt, q_set[ i]);
		}
			ratio[ i] = cost / cost_opt;
			if( ratio[ i] < stat_v.ratio_min)
				stat_v.ratio_min = ( float)ratio[ i];
			if( ratio[ i] > stat_v.ratio_max)
				stat_v.ratio_max = ( float)ratio[ i];

			ratio_sum += ( float)ratio[ i];
			stat_v.ratio_aver = ( float)ratio_sum / ( i+1);

			//Calculate the deviation.
			for( j=0; j<=i; j++)
				stat_v.ratio_dev += ( float)pow( ( ratio[ j] - stat_v.ratio_aver), 2);
			stat_v.ratio_dev /= i + 1;
			stat_v.ratio_dev = ( float)sqrt( stat_v.ratio_dev);
		

		//Print the accumulated query results.
		if( i == 0)
		{
			if( ( r_fp = fopen( COSKQ_RESULT_FILE, "w")) == NULL)
			{
				fprintf( stderr, "Cannot open the nskq_result file.\n");
				exit( 0);
			}
		}
		else
		{
			if( ( r_fp = fopen( COSKQ_RESULT_FILE, "a")) == NULL)
			{
				fprintf( stderr, "Cannot open the nskq_result file.\n");
				exit( 0);
			}
		}

		fprintf( r_fp, "Query #%i:\n", i+1);
		fprintf( r_fp, "Keywords: ");
		print_k_list( q_set[ i]->psi_v->k_head, r_fp);

		//Print the query result.
		if(cfg->alg_opt<=5)print_obj_set( S, r_fp);
		fprintf( r_fp, "Cost: %0.4lf\n\n", cost);
		printf( "Cost: %0.4lf\n", cost);
		fprintf( r_fp, "Time: %lf\n\n", usr_t);
		fprintf( r_fp, "Ratio: %lf\n\n", ratio[ i]);
		printf( "Time: %f\n", usr_t);
		printf( "Ratio: %lf\n\n", ratio[ i]);


		fclose( r_fp);

		//Print the statistics.
		stat_v.aver_cost = cost_sum / ( i+1);

#ifndef WIN32
		stat_v.q_time = usr_t_sum / ( i+1);
#endif
		
		if( cfg->alg_opt == 3)
		{
		//	printf( "feasible_set_n: %f\n", stat_v.feasible_set_n / ( i+1));
		//	printf( "node_set_n: %f\n", stat_v.node_set_n / ( i+1));
		}

		//
		if(cfg->alg_opt<=5)release_obj_set( S);
		if( S_opt != NULL)
			release_obj_set( S_opt);
	}

	//Collect the statistics.
	/*
	stat_v.q_time /= cfg->q_set_size;
	stat_v.aver_cost = cost / cfg->q_set_size;
	
	stat_v.memory_max /= cfg->q_set_size;
	stat_v.n_1_sum /= cfg->q_set_size;
	*/
	FILE* s_fp;
	if( !( s_fp = fopen( COSKQ_STAT_FILE, "a")))
	{
		fprintf( stderr, "Cannot open the nskq_stat file.\n");
		exit( 0);
	}
    fprintf(s_fp,"gnscal----- ratio_rag:%i ",ratio_tag);

	fclose( s_fp);
		print_coskq_stat( cfg, i+1);
	for( i=0; i<cfg->q_set_size; i++)
		release_query( q_set[ i]);
	for( i=0; i<cfg->q_set_size; i++){
			release_loc( q_set_ssq[i]->loc_v);
			release_psi( q_set_ssq[i]->psi_v);
			free( q_set_ssq[i]);
	}
	//Release the resources.
	free( q_set);
	free(q_set_ssq);
            }
        }
    }
	//printf("Memory balance: %f\n", stat_v.memory_v / cfg->q_set_size);
	free( cfg);
	release_data( data_v);
	free( MBR);
}

