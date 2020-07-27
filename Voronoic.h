#ifndef Voronoic_h
#define Voronoic_h

#include <list>
#include <queue>
#include <set>
#include <functional>
#include "VPoint.h"
#include "VEdge.h"
#include "VEllipse.h"
#include "VEventc.h"
#include "cao_alg.h"

namespace vorc
{
	/*
		Useful data containers for Vertices (places) and Edges of Voronoi diagram
	*/

	typedef std::list<VPoint *>		Vertices	;
	typedef std::list<VEdge *>		Edges		;
	
	/*
		Class for generating the Voronoi diagram
	*/

	class Voronoi
	{
	public:
		double key[3];
		Voronoi * next;
		/*
			Constructor - without any parameters
		*/

		Voronoi();

		/*
			The only public function for generating a diagram
			
			input:
				v		: Vertices - places for drawing a diagram
				w		: width  of the result (top left corner is (0, 0))
				h		: height of the result
			output:
				pointer to list of edges
			
			All the data structures are managed by this class
		*/
		
	//	Edges *			GetEdges(Vertices * v, int w, int h);


//	private:

		/*
						places		: container of places with which we work
						edges		: container of edges which will be teh result
						width		: width of the diagram
						height		: height of the diagram
						root		: the root of the tree, that represents a beachline sequence
						ly			: current "y" position of the line (see Fortune's algorithm)
		*/

		Vertices *		places;
		Edges *			edges;
		double			width, height;
		VEllipse *		root;
		double			rr;
		std::priority_queue<VEventc *, std::vector<VEventc *>, VEventc::CompareEvent> queue;
		/*
						deleted		: set  of deleted (false) Events (since we can not delete from PriorityQueue
						points		: list of all new points that were created during the algorithm
						queue		: priority queue with events to process
		*/

		std::set<VEventc *>	deleted;
		std::list<VPoint *> points;


		/*
						InsertParabola		: processing the place event
						RemoveParabola		: processing the circle event
						FinishEdge			: recursively finishes all infinite edges in the tree
						GetXOfEdge			: returns the current x position of an intersection point of left and right parabolas
						GetParabolaByX		: returns the Parabola that is under this "x" position in the current beachline
						CheckCircle			: checks the circle event (disappearing) of this parabola
						GetEdgeInterse
		*/

		void		InsertParabola(VPoint * p,query_ssq* q,double * costV);
		void		RemoveParabola(VEventc * e,query_ssq* q,double * costV);
		void		FinishEdge(VEllipse * n);
		double		GetAngleOfEdge(VEllipse * par, double r);
		VEllipse * GetParabolaByAngle(double angle,double an0,double an1,double rr);
	//	double		GetY(VPoint * p, double x);
		void		CheckCircle(VEllipse * b,double rr,query_ssq * q);
		VPoint *	GetEdgeIntersection(VEdge * a, VEdge * b,query_ssq * q);
	};
}

#endif
