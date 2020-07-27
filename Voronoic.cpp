#include <iostream>
#include <algorithm>
#include <cmath>
#include "Voronoic.h"

using namespace vorc;
double mod(double x){
	if(x<0.0)x=x+360.0;
	return x;
}
double GetAngle(double x,double y,double tx,double ty)
	{
		double angle;
		if((x-tx)==0 && (y-ty)>=0)return 90.0;
		if((x-tx)==0 && (y-ty)<0)return 270.0;
		
		angle=atan((y-ty)/(x-tx)) *180.0/3.141592653;
		if((x-tx)<0)angle=angle+180.0;
		else{if(angle<0)angle=360.0+angle;}
		
		return angle;
	}

Voronoi::Voronoi()
{
	edges = 0;
}
/*
Edges * Voronoi::GetEdges(Vertices * v, int w, int h)
{
	places = v;
	width = w;
	height = h;
	root = 0;
	std::cout<<"getedge";
	if(!edges) edges = new Edges();
	else 
	{
		for(Vertices::iterator	i = points.begin(); i != points.end(); ++i) delete (*i);
		for(Edges::iterator		i = edges->begin(); i != edges->end(); ++i) delete (*i);
		points.clear();
		edges->clear();
	}
	
	for(Vertices::iterator i = places->begin(); i!=places->end(); ++i)
	{
		queue.push(new VEventc( *i, true));
	}

	VEventc * e;
	while(!queue.empty())
	{
		e = queue.top();
		queue.pop();
		rr = e->point->r;
		if(deleted.find(e) != deleted.end()) { delete(e); deleted.erase(e); continue;}
		if(e->pe) InsertParabola(e->point);
		else {RemoveParabola(e);}
		delete(e);
	}
	
	FinishEdge(root);

	for(Edges::iterator i = edges->begin(); i != edges->end(); ++i)
	{
		if( (*i)->neighbour) 
		{
			(*i)->start = (*i)->neighbour->end;
			delete (*i)->neighbour;
		}
	}
	
	return edges;
}
*/
void	Voronoi::InsertParabola(VPoint * p,query_ssq* q,double * costV)
{
	if(!root){root = new VEllipse(p); return;}
	double a0=361.0,a1=0.0;
    double r;
	double angle,sqa,sqb,A,B,C,D,E,F,x1,x2,y1,y2,a,b,c,f,g,disc;
	VPoint * pc=0;


VEllipse * lm=root;
	while(!lm->isLeaf){
		lm=lm->Left();
	}
	VEllipse * rm=root;
	while(!rm->isLeaf){
		rm=rm->Right();
	}

	if(!root->isLeaf){
	
	
	if(rm->site==lm->site){
		a0=GetAngleOfEdge(lm->parent,p->r);
		a1=GetAngleOfEdge(rm->parent,p->r);
	}
	else{
	double al=GetAngleOfEdge(lm->parent,p->r);
	double ar=GetAngleOfEdge(rm->parent,p->r);
	
	double f1=(lm->site->x - rm->site->x) / (rm->site->y - lm->site->y);
	double g1=(lm->site->y+rm->site->y)/2.0-(lm->site->x+rm->site->x)*f1/2.0;

	pc=new VPoint((lm->site->x+p->tx)/2,(lm->site->y+p->ty)/2,q->loc_v->coord[0],q->loc_v->coord[1]);
	angle=lm->site->angle*3.141592653/180.0;
	sqa=p->r*p->r/4.0;
	sqb=sqa-lm->site->r*lm->site->r/4.0;

	 A=sqa*sin(angle)*sin(angle)+sqb*cos(angle)*cos(angle);
	 B=2.0*(sqb-sqa)*sin(angle)*cos(angle);
	 C=sqa*cos(angle)*cos(angle)+sqb*sin(angle)*sin(angle);
	 D=-2.0*A*pc->x-B*pc->y;
	 E=-1.0*B*pc->x - 2.0*C*pc->y;
	 F=A*pc->x*pc->x+B*pc->x*pc->y+C*pc->y*pc->y-sqa*sqb;
	// f=(p->y-p->ty)/(p->x-p->tx);
	// g=p->y-f*p->x;
	 a=A+B*f1+C*f1*f1;
	 b=B*g1+2.0*C*f1*g1+D+E*f1;
	 c=C*g1*g1+E*g1+F;

	 disc = b*b - 4.0 * a * c;
	 x1 = (-b + std::sqrt(disc)) / (2.0*a);
	 x2 = (-b - std::sqrt(disc)) / (2.0*a);
	 y1=f1*x1+g1;
	 y2=f1*x2+g1;

	 if( (mod(int (GetAngle(x1,y1,q->loc_v->coord[0],q->loc_v->coord[1])-ar)%360))> mod(int (al-ar)%360)) a0=GetAngle(x1,y1,q->loc_v->coord[0],q->loc_v->coord[1]);
	else a0=GetAngle(x2,y2,q->loc_v->coord[0],q->loc_v->coord[1]);}
	}
	/*
	if(root->isLeaf && ((p->r-root->site->r)<1)) // degenerovaný pøípad - obì spodní místa ve stejné výšce
	{
		VPoint * fp = root->site;
		root->isLeaf = false;
		root->SetLeft( new VEllipse(fp) );
		root->SetRight(new VEllipse(p)  );
		VPoint * s = new VPoint(0, 0); // zaèátek hrany uprostøed míst
		points.push_back(s);
		if(p->angle < fp->angle) root->edge = new VEdge(s, fp, p); // rozhodnu, který vlevo, který vpravo
		else root->edge = new VEdge(s, p, fp);
		edges->push_back(root->edge);
	//	std::cout << root->edge->start->x, root->edge->start->y, root->edge->end->x, root->edge->end->y;
		return;
	}
	*/
	VEllipse * par;
	if(a0==361.0)par=root;
	else{
	if(mod(p->angle-a0)<mod(a1-a0))
		{
			par=lm;
			//printf("%d,%d,par=lm\n",int(p->angle-a0)%360,int(a1-a0)%360);
	}
	else par = GetParabolaByAngle(p->angle,a0,a1,p->r);
	}
	//printf("insert:%f,%f,%f\n",p->angle,a0,a1);
	if(par->cEventc)
	{
		deleted.insert(par->cEventc);
		par->cEventc = 0;
		if(par==lm){
			deleted.insert(rm->cEventc);
			rm->cEventc = 0;
		}
		if(par==rm){
			deleted.insert(lm->cEventc);
			lm->cEventc = 0;
		}
	}
	pc=new VPoint((par->site->x+p->tx)/2,(par->site->y+p->ty)/2,q->loc_v->coord[0],q->loc_v->coord[1]);
	
	 angle=par->site->angle*3.141592653/180.0;
	
	 sqa=p->r*p->r/4.0;
	 sqb=sqa-par->site->r*par->site->r/4.0;

	 A=sqa*sin(angle)*sin(angle)+sqb*cos(angle)*cos(angle);
	 B=2.0*(sqb-sqa)*sin(angle)*cos(angle);
	 C=sqa*cos(angle)*cos(angle)+sqb*sin(angle)*sin(angle);
	 D=-2.0*A*pc->x-B*pc->y;
	 E=-1.0*B*pc->x - 2.0*C*pc->y;
	 F=A*pc->x*pc->x+B*pc->x*pc->y+C*pc->y*pc->y-sqa*sqb;
	 f=(p->y-p->ty)/(p->x-p->tx);
	 g=p->y-f*p->x;
	 a=A+B*f+C*f*f;
	 b=B*g+2.0*C*f*g+D+E*f;
	 c=C*g*g+E*g+F;

	 disc = b*b - 4.0 * a * c;
	 x1 = (-b + std::sqrt(disc)) / (2.0*a);
	 x2 = (-b - std::sqrt(disc)) / (2.0*a);
	 y1=f*x1+g;
	 y2=f*x2+g;

	double rx,ry;

	if(abs(GetAngle(x1,y1,q->loc_v->coord[0],q->loc_v->coord[1])-p->angle)<10){
		rx=x1;
		ry=y1;
	}
	else{
		ry=y2;
		rx=x2;
	}
	VPoint * start = new VPoint(rx,ry,q->loc_v->coord[0],q->loc_v->coord[1]);
	points.push_back(start);
	
	VEdge * el = new VEdge(start, p, par->site);
	VEdge * er = new VEdge(start, par->site,p);
//	std::cout << el->start->x, el->start->y, el->end->x, el->end->y;
	el->neighbour = er;
//	edges->push_back(el);
	par->edge = er;
	par->isLeaf = false;

	VEllipse * p0 = new VEllipse(par->site);
	VEllipse * p1 = new VEllipse(p);
	VEllipse * p2 = new VEllipse(par->site);

	par->SetRight(p2);
	par->SetLeft(new VEllipse());
	par->Left()->edge = el;

	par->Left()->SetLeft(p0);
	par->Left()->SetRight(p1);


	B_KEY_TYPE costV_1;
		loc_t* loc_v=alloc_loc(2);
		//loc_v->dim=2;
		loc_v->coord[0]=(p->x+par->site->x)/2.0;
		loc_v->coord[1]=(p->y+par->site->y)/2.0;
		obj_set_t * V_1,* V;
		obj_t * obj_v;
//		alloc_obj(obj_v,2);
	query_t * q_new = alloc_query( );
	q_new->loc_v = loc_v;

			q_new->psi_v = alloc_psi( );
			copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
		//	V_1 = Cao_Appro1( q_new);
			//add_obj_set_entry( obj_v, V_1);

		//	if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}
			

			disk_t * disk_v=alloc_disk(2);
            r=std::sqrt(pow(p->x-par->site->x,2.0)+pow(p->y-par->site->y,2.0))/2.0+1.0;
			set_disk(disk_v,loc_v,r);

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);

			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
			
			//costV_1 = comp_cost( V_1, q_new);
         
			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);
			
			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
//				V_tmp = V;
//				V = V_1;
				
//				release_obj_set( V_tmp);
			}
			
				
			}
			release_query( q_new);		
	//		release_obj_set( V_1);
		//	release_obj(obj_v);

	// pøestavuju strom .. vkládám novou parabolu

//	release_loc(loc_v);
	CheckCircle(p0,p->r,q);
	CheckCircle(p2,p->r,q);
}

void	Voronoi::RemoveParabola(VEventc * e,query_ssq* q,double * costV)
{
	VEllipse * p1 = e->arch;
    double r;
	if(!p1->cEventc)return;
	if(!e)return;
	if(p1->cEventc!=e)return;
	VEllipse * xl = VEllipse::GetLeftParent(p1);
	VEllipse * xr = VEllipse::GetRightParent(p1);

	VEllipse * p0 =  VEllipse::GetLeftChild(xl);
	VEllipse * p2 =  VEllipse::GetRightChild(xr);

	VEllipse * lm=root;
	while(!lm->isLeaf){
		lm=lm->Left();
	}

	VEllipse * rm=root;
	while(!rm->isLeaf){
		rm=rm->Right();
	}

	B_KEY_TYPE costV_1;
	obj_set_t * V_1;
	obj_t * obj_v;
	query_t * q_new;
	disk_t * disk_v;

	if(!p0){
		if(p2->cEventc){ deleted.insert(p2->cEventc); p2->cEventc = 0; }
		if(rm->cEventc){deleted.insert(rm->cEventc); rm->cEventc = 0;}
		xl=rm->parent;
		p0 =  VEllipse::GetLeftChild(xl);
		if(p0->cEventc){ deleted.insert(p0->cEventc); p0->cEventc = 0; }
		VPoint * p=GetEdgeIntersection(xl->edge, xr->edge,q);
	//	if(p==NULL)return;
	points.push_back(p);

	xl->edge->end = p;
	xr->edge->end = p;



	//	B_KEY_TYPE costV_1;
		loc_t* loc_v=alloc_loc(2);
		//loc_v->dim=2;
		loc_v->coord[0]=p->x;
		loc_v->coord[1]=p->y;
		
	q_new = alloc_query( );
	q_new->loc_v = loc_v;

			q_new->psi_v = alloc_psi( );
			copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			//V_1 = Cao_Appro1( q_new);
//			add_obj_set_entry( obj_v, V_1);

		//	if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}

			disk_t * disk_v=alloc_disk(2);
            r=std::sqrt(pow(p->x-p2->site->x,2.0)+pow(p->y-p2->site->y,2.0))+1.0;
			set_disk(disk_v,loc_v,r);

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
			//costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
			//	V = V_1;
				
			//	release_obj_set( V_tmp);
			}

				
			}
			loc_v->coord[0]=(p0->site->x+p2->site->x)/2.0;
			loc_v->coord[1]=(p0->site->y+p2->site->y)/2.0;
		
	//q_new = alloc_query( );
	q_new->loc_v = loc_v;

		//	q_new->psi_v = alloc_psi( );
		//	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			//V_1 = Cao_Appro1( q_new);
			//add_obj_set_entry( obj_v, V_1);

			//if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}


            r=std::sqrt(pow(p0->site->x-p2->site->x,2.0)+pow(p0->site->y-p2->site->y,2.0))/2.0+1.0;
			set_disk(disk_v,loc_v,r);

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
			//costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
			//	V = V_1;
				
			//	release_obj_set( V_tmp);
			}

				
			}
			release_query( q_new);		
		//	release_obj_set( V_1);



	VEllipse * higher=new VEllipse();
	//VEllipse * par = p1;
	//while(par != root)
	//{
	//	par = par->parent;
	//	if(par == xl) higher = xl;
	//	if(par == xr) higher = xr;
	//}
	higher=xl;
	higher->edge = new VEdge(p, p2->site, p0->site);
	//higher->SetLeft(p0);
	VEllipse * p21=new VEllipse(p2->site);
	higher->SetRight(p21);
//	edges->push_back(higher->edge);
	if(xr==root){
		xr->Right()->parent=0;
		root=xr->Right();}
	else{
//	xr->Right()->parent=xr->parent;
	xr->parent->SetLeft(xr->Right());
	}
//	xl->edge=new VEdge(p, p2->site, p0->site);
//	VEllipse * gparent = p1->parent->parent->Left();
//	xr->isLeaf=true;
//	xr->site=p2->site;
//	rm->site=p2->site;
//	rm->isLeaf=true;


	//VEllipse * gp=xl->parent;
//	gp->SetRight(higher);
//	delete p2;
//	delete p1->parent;
	delete p1;
//	delete rm->parent;
//	delete rm;
	delete rm;
	delete xr;

	CheckCircle(p0,e->point->r,q);
	CheckCircle(p2,e->point->r,q);
	return;
	}


		if(!p2){
		if(p0->cEventc){ deleted.insert(p0->cEventc); p0->cEventc = 0; }
		if(lm->cEventc){deleted.insert(lm->cEventc); lm->cEventc = 0;}
			xr=VEllipse::GetRightParent(lm);
			p2 =  VEllipse::GetRightChild(xr);
			if(p2->cEventc){ deleted.insert(p2->cEventc); p2->cEventc = 0; }
		VPoint * p=GetEdgeIntersection(xl->edge, xr->edge,q);

	points.push_back(p);

	xl->edge->end = p;
	xr->edge->end = p;

	VEllipse * higher=new VEllipse();
	higher=xr;
	higher->edge = new VEdge(p, p2->site, p0->site);
	higher->SetLeft(new VEllipse(p0->site));
	if(xl==root){
		xl->Left()->parent=0;
		root=xl->Left();}
	else{
//	xl->Left()->parent=xl->parent;
	xl->parent->SetRight(xl->Left());
	
	}
//	VEllipse * gp=xr->parent;
	//gp->SetLeft(higher);
	

//			B_KEY_TYPE costV_1;
		loc_t* loc_v;
		loc_v=alloc_loc(2);
		loc_v->coord[0]=p->x;
		loc_v->coord[1]=p->y;
//		obj_set_t * V_1,* V,* V_tmp;
//		obj_t * obj_v;
	q_new = alloc_query( );
	q_new->loc_v = loc_v;

			q_new->psi_v = alloc_psi( );
			copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			//V_1 = Cao_Appro1( q_new);
			//add_obj_set_entry( obj_v, V_1);

		//	if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}

			disk_v=alloc_disk(2);
            r=std::sqrt(pow(p->x-p2->site->x,2.0)+pow(p->y-p2->site->y,2.0))+1.0;
			set_disk(disk_v,loc_v,r);
			
			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
			//costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
			//	V = V_1;
				
				//release_obj_set( V_tmp);
			}
			}

				loc_v->coord[0]=(p0->site->x+p2->site->x)/2.0;
			loc_v->coord[1]=(p0->site->y+p2->site->y)/2.0;
		
	//q_new = alloc_query( );
	q_new->loc_v = loc_v;

		//	q_new->psi_v = alloc_psi( );
		//	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			//V_1 = Cao_Appro1( q_new);
		//	add_obj_set_entry( obj_v, V_1);

		//	if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}

			disk_t * disk_v=alloc_disk(2);
            r=std::sqrt(pow(p0->site->x-p2->site->x,2.0)+pow(p0->site->y-p2->site->y,2.0))/2.0+1.0;
			set_disk(disk_v,loc_v,r);

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
			//costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
			//	V = V_1;
				
			//	release_obj_set( V_tmp);
			}

				
			}

			release_query( q_new);		
		//	release_obj_set( V_1);
//	delete p1->parent;
	delete p1;
//	delete lm->parent;
	delete lm;
	delete xl;
	

	CheckCircle(p2,e->point->r,q);
	CheckCircle(p0,e->point->r,q);
	return;
	}



	if(p0 == p2) std::cout << "chyba - pravá a levá parabola má stejné ohnisko!\n";

	if(p0->cEventc){ deleted.insert(p0->cEventc); p0->cEventc = 0; }
	if(p2->cEventc){ deleted.insert(p2->cEventc); p2->cEventc = 0; }

	VPoint * p=GetEdgeIntersection(xl->edge, xr->edge,q);

	points.push_back(p);

	xl->edge->end = p;
	xr->edge->end = p;

	VEllipse * higher;
	VEllipse * par = p1;
	while(par != root)
	{
		par = par->parent;
		if(par == xl) higher = xl;
		if(par == xr) higher = xr;
	}
	higher->edge = new VEdge(p, p2->site, p0->site);
//	edges->push_back(higher->edge);
	   
	VEllipse * gparent = p1->parent->parent;
	if(p1->parent->Left() == p1)
	{
		if(gparent->Left()  == p1->parent) gparent->SetLeft ( p1->parent->Right() );
		if(gparent->Right() == p1->parent) gparent->SetRight( p1->parent->Right() );
	}
	else
	{
		if(gparent->Left()  == p1->parent) gparent->SetLeft ( p1->parent->Left()  );
		if(gparent->Right() == p1->parent) gparent->SetRight( p1->parent->Left()  );
	}

	delete p1->parent;
	delete p1;

		loc_t* loc_v=alloc_loc(2);
		
		loc_v->coord[0]=p->x;
		loc_v->coord[1]=p->y;

	q_new = alloc_query( );
	q_new->loc_v = loc_v;

			q_new->psi_v = alloc_psi( );
			copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
			//V_1 = Cao_Appro1( q_new);
			//add_obj_set_entry( obj_v, V_1);

			//if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}

			disk_v=alloc_disk(2);
			set_disk(disk_v,loc_v,r);
            r=std::sqrt(pow(p->x-p2->site->x,2.0)+pow(p->y-p2->site->y,2.0))/4.0;

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
//			costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
				//V = V_1;
				
			//	release_obj_set( V_tmp);
			}
			}
						loc_v->coord[0]=(p0->site->x+p2->site->x)/2.0;
			loc_v->coord[1]=(p0->site->y+p2->site->y)/2.0;
		
	//q_new = alloc_query( );
	q_new->loc_v = loc_v;

		//	q_new->psi_v = alloc_psi( );
		//	copy_k_list( q_new->psi_v->k_head, q->psi_v->k_head); //obj_v->k_head);

			//Solve the new query.
		//	V_1 = Cao_Appro1( q_new);
		//	add_obj_set_entry( obj_v, V_1);

		//	if( V_1 == NULL)
		//	{
		//		release_query( q_new);
			//	continue;
		//	}

			set_disk(disk_v,loc_v,r);
            r=std::sqrt(pow(p0->site->x-p2->site->x,2.0)+pow(p0->site->y-p2->site->y,2.0))/2.0+1.0;

			obj_v = const_NN_key( q->loc_v, q->target, disk_v);
			if(obj_v!=NULL){
			q_new->loc_v=get_obj_loc(obj_v);
//			costV_1 = comp_cost( V_1, q_new);

			costV_1=2*r+calc_dist_loc(get_obj_loc(obj_v),q->loc_v);

			
			if( costV_1 < *costV)
			{
				*costV = costV_1;
			//	V_tmp = V;
			//	V = V_1;
				
			//	release_obj_set( V_tmp);
			}

				
			}
			release_query( q_new);		
		//	release_obj_set( V_1);
	CheckCircle(p0,e->point->r,q);
	CheckCircle(p2,e->point->r,q);
}

void	Voronoi::FinishEdge(VEllipse * n)
{
	if(n->isLeaf) {delete n; return;}
	double mx;
	if(n->edge->direction->x > 0.0)	mx = max(width,	n->edge->start->x + 10 );
	else							mx = min(-width,		n->edge->start->x - 10);
	
	VPoint * end = new VPoint(mx, mx * n->edge->f + n->edge->g); 
	n->edge->end = end;
	points.push_back(end);
			
	FinishEdge(n->Left() );
	FinishEdge(n->Right());
	delete n;
}

double	Voronoi::GetAngleOfEdge(VEllipse * par, double rr)
{
	VPoint * p= par->edge->left;
	VPoint * r=par->edge->right;
	double f=par->edge->f;
	double g=par->edge->g;

	VPoint * pc=new VPoint((p->x+p->tx)/2.0,(p->y+p->ty)/2.0,p->tx,p->ty);

	double angle=p->angle*3.141592653/180.0;
	double sqa=rr*rr;
	sqa=sqa/4.0;
	double sqb=sqa-p->r*p->r/4.0;

	double A=sqa*pow(sin(angle),2.0)+sqb*pow(cos(angle),2.0);
	double B=2*(sqb-sqa)*sin(angle)*cos(angle);
	double C=sqa*pow(cos(angle),2.0)+sqb*pow(sin(angle),2.0);
	double D=-2.0*A*pc->x-B*pc->y;
	double E=-1.0*B*pc->x - 2.0*C*pc->y;
	double F=A*pc->x*pc->x+B*pc->x*pc->y+C*pc->y*pc->y-sqa*sqb;

	double a=A+B*f+C*f*f;
	double b=B*g+2.0*C*f*g+D+E*f;
	double c=C*g*g+E*g+F;

	double disc = b*b - 4.0 * a * c;
	double x1 = (-b + std::sqrt(disc)) / (2.0*a);
	double x2 = (-b - std::sqrt(disc)) / (2.0*a);
	double y1=f*x1+g;
	double y2=f*x2+g;
	
	//double ry;
	if((x1 - par->edge->start->x)/par->edge->direction->x < 0) return GetAngle(x2,y2,p->tx,p->ty);
	if((y1 - par->edge->start->y)/par->edge->direction->y < 0) return GetAngle(x2,y2,p->tx,p->ty);
		
	
	//if(p->r>r->r)ry=std::max(GetAngle(x1,y1),GetAngle(x2,y2));
	//else ry=std::min(GetAngle(x1,y1),GetAngle(x2,y2));

	return GetAngle(x1,y1,p->tx,p->ty);
}

VEllipse * Voronoi::GetParabolaByAngle(double angle,double an0, double an1, double rr)
{
	VEllipse * par = root;
	double x = 0.0;
	if(an0==361.0){
		return root;
	}
	else{
		
	while(!par->isLeaf) // projdu stromem dokud nenarazím na vhodný list
	{
		x = GetAngleOfEdge(par, rr);
	//	printf("angles:%f,%f,%f,%f\n",x, angle,an0,an1);
	//	if(int(x-an0)%360<int(an1-an0)%360)return 0;
	//	if(x==an0)x=x-1;
		//if(x>270)x=x-360.0;
		if(mod(x-an1)<mod(angle-an1)) par = par->Left();
		else par = par->Right();				
	}}
	
	return par;
}


void	Voronoi::CheckCircle(VEllipse * b,double rr,query_ssq * q)
{
	VEllipse * lp = VEllipse::GetLeftParent (b);
	VEllipse * rp = VEllipse::GetRightParent(b);

	VEllipse * a  = VEllipse::GetLeftChild (lp);
	VEllipse * c  = VEllipse::GetRightChild(rp);

	VEllipse * lm=root;
	while(!lm->isLeaf){
		lm=lm->Left();
	}
	VEllipse * rm=root;
	while(!rm->isLeaf){
		rm=rm->Right();
	}

	if(!a){
		lp=VEllipse::GetLeftParent(rm);
		a  = VEllipse::GetLeftChild (lp);
	}

	if(!c){
		rp=VEllipse::GetRightParent(lm);
		c  = VEllipse::GetRightChild(rp);
	}
		



	if(a && c && a->site == c->site) return;

	VPoint * s = 0;
	s = GetEdgeIntersection(lp->edge, rp->edge,q);
	if(s == 0) return;

	double dx = lp->edge->left->x - s->x;
	double dy = lp->edge->left->y - s->y;

	double d = std::sqrt( (dx * dx) + (dy * dy) );

	if(s->r + d < rr) { return;}
	double dp=std::sqrt(pow(s->x-s->tx,2)+pow(s->y-s->ty,2));
	double ex=s->x+d*(s->x-s->tx)/dp;
	double ey=s->y+d*(s->y-s->ty)/dp;
	VEventc * e = new VEventc(new VPoint(ex,ey,q->loc_v->coord[0],q->loc_v->coord[1]), false);
	points.push_back(e->point);
	b->cEventc = e;
	e->arch = b;
	queue.push(e);
}

VPoint * Voronoi::GetEdgeIntersection(VEdge * a, VEdge * b, query_ssq * q)
{
	double x = (b->g-a->g) / (a->f - b->f);
	double y = a->f * x + a->g;

	if((x - a->start->x)/a->direction->x < 0) return 0;
	if((y - a->start->y)/a->direction->y < 0) return 0;
		
	if((x - b->start->x)/b->direction->x < 0) return 0;
	if((y - b->start->y)/b->direction->y < 0) return 0;	

	VPoint * p = new VPoint(x, y,q->loc_v->coord[0],q->loc_v->coord[1]);		
	points.push_back(p);
	return p;
}


	
