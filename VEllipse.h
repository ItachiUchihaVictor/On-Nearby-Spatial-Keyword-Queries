#ifndef VEllipse_h
#define VEllipse_h

#include "VPoint.h"
#include "VEdge.h"
//#include "VEvent.h"
//#include "VEvent_c.h"
/*
	Parabolas and events have pointers to each other, so we declare class VEvent, which will be defined later.
*/

class VEventc;

/*
	A class that stores information about an item in beachline sequence (see Fortune's algorithm). 
	It can represent an arch of parabola or an intersection between two archs (which defines an edge).
	In my implementation I build a binary tree with them (internal nodes are edges, leaves are archs).
*/

class VEllipse
{
public:

	/*
		isLeaf		: flag whether the node is Leaf or internal node
		site		: pointer to the focus point of parabola (when it is parabola)
		edge		: pointer to the edge (when it is an edge)
		cEvent		: pointer to the event, when the arch disappears (circle event)
		parent		: pointer to the parent node in tree
	*/

	bool		isLeaf;
	VPoint *	site;
	VEdge *		edge;
	VEventc *	cEventc;
	VEllipse * parent;

	/*
		Constructors of the class (empty for edge, with focus parameter for an arch).
	*/
	
	VEllipse	();
	VEllipse	(VPoint * s);

	/*
		Access to the children (in tree).
	*/

	void		SetLeft (VEllipse * p) {_left  = p; p->parent = this;}
	void		SetRight(VEllipse * p) {
		_right = p; 
		p->parent = this;}

	VEllipse *	Left () { return _left;  }
	VEllipse * Right() { return _right; }

	/*
		Some useful tree operations

		GetLeft			: returns the closest left leave of the tree
		GetRight		: returns the closest right leafe of the tree
		GetLeftParent	: returns the closest parent which is on the left
		GetLeftParent	: returns the closest parent which is on the right
		GetLeftChild	: returns the closest leave which is on the left of current node
		GetRightChild	: returns the closest leave which is on the right of current node
	*/

	static VEllipse * GetLeft			(VEllipse * p);
	static VEllipse * GetRight			(VEllipse * p);
	static VEllipse * GetLeftParent	(VEllipse * p);
	static VEllipse * GetRightParent	(VEllipse * p);
	static VEllipse * GetLeftChild		(VEllipse * p);
	static VEllipse * GetRightChild	(VEllipse * p);
	
private:
	VEllipse * _left;
	VEllipse * _right;

};

#endif