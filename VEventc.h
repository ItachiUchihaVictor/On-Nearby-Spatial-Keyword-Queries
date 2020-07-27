#ifndef VEventc_h
#define VEventc_h

#include <iostream>
#include "VPoint.h"
#include "VEllipse.h"

/*
	The class for storing place / circle event in event queue.

	point		: the point at which current event occurs (top circle point for circle event, focus point for place event)
	pe			: whether it is a place event or not
	y			: y coordinate of "point", events are sorted by this "y"
	arch		: if "pe", it is an arch above which the event occurs
*/

class VEventc
{
public:
	VPoint *	point;
	bool		pe;
	double		y;
	double		x;
	VEllipse * arch;

	/*
		Constructor for the class
		
		p		: point, at which the event occurs
		pev		: whether it is a place event or not
	*/
	
	VEventc(VPoint * p, bool pev)
	{
		point	= p;
		pe		= pev;
		y		= p->y;
		x		= p->x;
		arch	= 0;
	}

	/*
		function for comparing two events (by "y" property)
	*/

	struct CompareEvent : public std::binary_function<VEventc*, VEventc*, bool>
	{
		bool operator()(const VEventc* l, const VEventc* r) const { return (l->point->r > r->point->r); }
	};
};

#endif
