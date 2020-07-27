/* $Id$
 * 
 * Copyright (C) 2007 Jose Alberto Cisneros Perez. All rights reserved.
 *
 * This file may be used under the terms of the GNU General Public License
 * versions 3.0 as published by the Free Software Foundation and
 * appearing in the COPYING file included in the packaging of this file.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 *  WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *  This file implements a concatenable queue where the information is 
 *  stored in the leaves of a 2-3 Tree.
*/

#ifndef __CQUEUE_H__
#define __CQUEUE_H__


#include <iostream>



enum node_type_d
{
	LEAF_NODE = 0 ,
	TWO_NODE,
	THREE_NODE,
};

struct node_d 
{
	node_d() : type(LEAF_NODE),
					data1(0.0),
					data2(0.0),
					x_coord(0.0),
					left(NULL),
					middle(NULL),
					right(NULL),
					parent(NULL),
					biggest_y(NULL),
					smallest_y(NULL)
	{
	}

	node_type_d type; 
	double data1, data2, x_coord; // on leaf nodes, data2 is the y-coordinate of the point 
	node_d *left, *middle, *right, *parent, *biggest_y, *smallest_y;//
}; 

class ConcatenableQueue 
{
public:
	ConcatenableQueue();
	ConcatenableQueue(const double &x_coord, const double &y_coord);
	void addNode(const double &x_coord, const double &y_coord);
	void deleteNode(const double &value);
	bool isEmpty() const;
	void print();
	void search(const double &value);
	int height() const;
	void addChildNode(node_d *child, node_d *aux, const double &value);
	node_d *rightmostNodeAtLevel(const int &level, node_d*);
	node_d *leftmostNodeAtLevel(const int &level, node_d*);	
	node_d *root();
	void set_root(node_d *new_root);
	ConcatenableQueue clone();

	//Procedures concatenate and split concatenable queues
	static ConcatenableQueue concatenate(ConcatenableQueue CQ1, ConcatenableQueue CQ2);
	//split_left, CQ1 keeps the node with value val 
	static void split_left (ConcatenableQueue CQ, ConcatenableQueue *CQ1, ConcatenableQueue *CQ2, const double &split_value);

	//split_right, CQ2 keeps the node with value val 
	static void split_right (ConcatenableQueue CQ, ConcatenableQueue *CQ1, ConcatenableQueue *CQ2, const double &split_value);

private:
	node_d *createLeafNode(const double &x_coord, const double &y_coord);
	node_d *create2Node(const double &l_value, const double &m_value, node_d *left_child, node_d *right_child);
	node_d *searchForInsert(node_d *new_node, node_d *current_node);
	node_d *searchLeafNode(const double &value, node_d *current);
	void deleteNode(node_d *node_to_delete);
	void updateValuesInsert(node_d *current_node, double value);
	void updateValuesDelete(node_d *current_node, const double &current_value, const double &new_value);
	void printValues(node_d *current_node);
	void addQueuePointsToQueue(ConcatenableQueue *CQ, node_d *current_node);

	node_d *root_;
};

#endif /*__CQUEUE_H__*/
