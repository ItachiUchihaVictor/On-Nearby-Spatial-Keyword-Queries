#include "VEllipse.h"
#include "VPoint.h"

/*
	Constructors
*/

VEllipse::VEllipse()
{
	site	= 0;
	isLeaf	= false;
	cEventc	= 0;
	edge	= 0;
	parent	= 0;
	_left=0;
	_right=0;
}

VEllipse::VEllipse(VPoint * s)
{
	site	= s; 
	isLeaf	= true;
	cEventc	= 0;
	edge	= 0;
	parent	= 0;
	_left=0;
	_right=0;
}

/*
	Tree operations (described in the header file)
*/

VEllipse * VEllipse::GetLeft			(VEllipse * p)
{
	return GetLeftChild(GetLeftParent(p));
}


VEllipse * VEllipse::GetRight			(VEllipse * p)
{
	return GetRightChild(GetRightParent(p));
}

VEllipse * VEllipse::GetLeftParent	(VEllipse * p)
{
	VEllipse * par		= p->parent;
	VEllipse * pLast	= p;
	while(par->Left() == pLast) 
	{ 
		if(!par->parent) return 0;
		pLast = par; 
		par = par->parent; 
	}
	return par;
}

VEllipse * VEllipse::GetRightParent	(VEllipse * p)
{
	VEllipse * par		= p->parent;
	VEllipse * pLast	= p;
	while(par->Right() == pLast) 
	{ 
		if(!par->parent) return 0;
		pLast = par; par = par->parent; 
	}
	return par;
}

VEllipse * VEllipse::GetLeftChild		(VEllipse * p)
{
	if(!p) return 0;
	VEllipse * par = p->Left();
	while(!par->isLeaf) par = par->Right();
	return par;
}

VEllipse * VEllipse::GetRightChild	(VEllipse * p)
{
	if(!p) return 0;
	VEllipse * par = p->Right();
	while(!par->isLeaf) par = par->Left();
	return par;
}
